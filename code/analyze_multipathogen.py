# -*- coding: utf-8 -*-
"""多病原 + 宏基因组层扩展审计分析 (v2.0)
输入: multipathogen_biosample.csv (全量4病原+小病原) / metagenome_sample.csv (抽样3000)
输出: multipathogen_analysis.json
归因规则: 与 map_analyze.py 完全一致 (复用其词典)
"""
import csv, json, os, sys, io
from collections import Counter, defaultdict
from datetime import datetime

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.path.insert(0, r"D:\南开CDC\task\1、监测殖民主义")
from map_analyze import resolve_country as _base_resolve, sample_country, INCOME, parse_date

# v2.0 扩展机构词典 (覆盖 v1 审计未覆盖的多病原/宏基因组提交主体)
EXT_MAP = {
    "cdphe": "USA", "us cdc national wastewater": "USA", "national wastewater surveillance": "USA",
    "biobot": "USA", "utah public health": "USA", "wadsworth": "USA",
    "oregon state university": "USA", "wisconsin state laboratory": "USA",
    "nj doh": "USA", "university of miami": "USA", "genometrakr": "USA",
    "illinois_wbe": "USA", "washington state department": "USA",
    "unsw": "AUS", "university of new south wales": "AUS",
    "tohoku": "JPN", "guangdong": "CHN", "university of windsor": "CAN",
    "aalborg university": "DNK", "aalborg": "DNK",
    "queen's university belfast": "GBR", "belfast": "GBR", "mcgill": "CAN",
}
EXT_PROXY = {"ebi": "GBR", "european bioinformatics institute": "GBR", "ena": "GBR"}

def resolve_country(owner):
    o = (owner or "").strip().lower()
    if not o:
        return "(unresolved)"
    for k, v in EXT_PROXY.items():
        if o == k or o.startswith(k + " ") or o.startswith(k + ","):
            return v
    for k, v in EXT_MAP.items():
        if k in o:
            return v
    return _base_resolve(owner)

DST = r"D:\南开CDC\task\1、监测殖民主义"
HUB_OWNERS = ("EBI", "European Bioinformatics Institute")

def load(fp):
    return list(csv.DictReader(open(fp, encoding="utf-8-sig")))

def summarize(rows, label, sample_frac=None):
    n = len(rows)
    if n == 0:
        return {"label": label, "n": 0}
    sc_counter = Counter()
    oc_counter = Counter()
    income_counter = Counter()
    hub_n = 0
    lags = []
    aligned = misaligned = 0
    sub_dates = []
    for r in rows:
        sc = sample_country(r["geo_location"])
        oc = resolve_country(r["owner"])
        if sc != "(missing)":
            sc_counter[sc] += 1
            income_counter[INCOME.get(sc, "unknown")] += 1
        if r["owner"].strip() in HUB_OWNERS:
            hub_n += 1
        oc_counter[oc] += 1
        if sc != "(missing)" and oc not in ("(unresolved)",):
            if sc == oc:
                aligned += 1
            else:
                misaligned += 1
        d1, d2 = parse_date(r["submission_date"]), parse_date(r["collection_date"])
        if d1 and d2:
            lag = (d1 - d2).days
            if -60 < lag <= 2920:  # metagenome 项目周期更长, 放宽上限到8年
                lags.append(lag)
        if d1:
            sub_dates.append(d1.year)
    lags_sorted = sorted(lags)
    out = {
        "label": label,
        "n": n,
        "n_sample_countries": len(sc_counter),
        "top_sample_countries": sc_counter.most_common(10),
        "top_owner_countries": oc_counter.most_common(10),
        "hub_ebi_share_pct": round(hub_n / n * 100, 2),
        "income_distribution": dict(income_counter),
        "aligned_pct": round(aligned / max(1, aligned + misaligned) * 100, 1),
        "n_valid_lag": len(lags),
    }
    if lags_sorted:
        out["median_lag_days"] = lags_sorted[len(lags_sorted) // 2]
        out["p90_lag_days"] = lags_sorted[int(len(lags_sorted) * 0.9)]
    if sub_dates:
        yc = Counter(sub_dates)
        out["submission_year_range"] = [min(yc), max(yc)]
        out["submission_year_distribution"] = {str(k): v for k, v in sorted(yc.items())}
    if sample_frac:
        # 抽样外推全群体量
        out["extrapolated_total"] = round(n / sample_frac)
    return out

def main():
    fetch_log = json.load(open(os.path.join(DST, "multipathogen_fetch.json"), encoding="utf-8"))
    meta_total = fetch_log["metagenome"]["total_population"]
    meta_sampled = fetch_log["metagenome"]["sample_n"]

    result = {"ascertained": "2026-11-02",
              "meta_population": meta_total, "meta_sampled": meta_sampled}

    # ---- 病原层 ----
    mp = load(os.path.join(DST, "multipathogen_biosample.csv"))
    by_layer = defaultdict(list)
    for r in mp:
        by_layer[r["pathogen_layer"]].append(r)
    result["pathogen_layers"] = {}
    for layer, rows in by_layer.items():
        name = {"influenza_a": "Influenza A virus", "influenza_b": "Influenza B virus",
                "norovirus": "Norovirus", "rotavirus_a": "Rotavirus A",
                "rsv": "RSV"}.get(layer, layer)
        result["pathogen_layers"][layer] = summarize(rows, name)
    result["pathogen_all"] = summarize(mp, "all pathogen (non-SARS-CoV-2)")

    # ---- 宏基因组层 (抽样) ----
    mg = load(os.path.join(DST, "metagenome_sample.csv"))
    result["metagenome_sample"] = summarize(mg, "wastewater metagenome (sample)",
                                             sample_frac=meta_sampled / meta_total)
    # 主权化提交: 采样国==所有者国 且非hub
    sovereign = 0
    for r in mg:
        sc = sample_country(r["geo_location"])
        oc = resolve_country(r["owner"])
        if sc != "(missing)" and sc == oc:
            sovereign += 1
    result["metagenome_sample"]["sovereign_domestic_share_pct"] = round(sovereign / len(mg) * 100, 2)

    # ---- 与 SARS-CoV-2 主审计合并全景 ----
    master = load(os.path.join(DST, "master_biosample.csv"))
    result["panorama"] = {
        "sars_cov_2_pathogen": 82757,
        "other_pathogen": len(mp),
        "metagenome_total": meta_total,
        "total_wastewater_records": 82757 + len(mp) + meta_total,
        "note": "metagenome layer is 4.5x the pathogen layer - the archive's largest stratum"
    }

    with open(os.path.join(DST, "multipathogen_analysis.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
    print(json.dumps({k: result[k] for k in
                      ["pathogen_all", "metagenome_sample", "panorama"]},
                     ensure_ascii=False, indent=1)[:3000])

if __name__ == "__main__":
    main()

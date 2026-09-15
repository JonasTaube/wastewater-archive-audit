# -*- coding: utf-8 -*-
"""错位矩阵分析: 采样国 x 提交机构国 -> 数据提取不对称性审计
输入: master_biosample.csv (fetch_full.py 产出)
输出: analysis_summary.json / misalignment_matrix.csv / country_stats.csv
"""
import csv, json, os, re, sys
from collections import Counter, defaultdict
from datetime import datetime

DST = r"D:\南开CDC\task\1、监测殖民主义"
MASTER = os.path.join(DST, "master_biosample.csv")

# ---------- 1. 机构名 -> 国别 ----------
# 数据库代提交机构
PROXY_DB = {
    "ebi": "GBR", "embl-ebi": "GBR", "european bioinformatics institute": "GBR",
    "ena": "GBR", "ncbi": "USA", "ddbj": "JPN", "dna data bank of japan": "JPN",
}
# 已知机构映射 (精确/子串匹配, 小写)
ORG_MAP = {
    "broad institute": "USA", "harvard": "USA", "mit": "USA",
    "wellcome sanger": "GBR", "sanger institute": "GBR",
    "public health england": "GBR", "uk health security": "GBR", "ukhsa": "GBR",
    "nibsc": "GBR", "university of cambridge": "GBR", "imperial college": "GBR",
    "university college london": "GBR", "london school of hygiene": "GBR",
    "glassgow": "GBR", "edinburgh": "GBR", "oxford": "GBR",
    "centers for disease control": "USA", "cdc atlanta": "USA",
    "national institutes of health": "USA", "nih ": "USA",
    "washington university": "USA", "university of washington": "USA",
    "ucsd": "USA", "ucsf": "USA", "ucla": "USA", "uc davis": "USA",
    "stanford": "USA", "yale": "USA", "cornell": "USA", "duke": "USA",
    "johns hopkins": "USA", "new york university": "USA", "nyu": "USA",
    "mount sinai": "USA", "scripps": "USA", "illumina": "USA",
    "abbott": "USA", "thermo fisher": "USA", "helix": "USA",
    "university of michigan": "USA", "university of wisconsin": "USA",
    "university of california": "USA", "colorado state": "USA",
    "battelle": "USA", "verily": "USA",
    "institut pasteur": "FRA", "pasteur institute": "FRA",
    "sorbonne": "FRA", "universite de paris": "FRA",
    "charite": "DEU", "tu muenchen": "DEU", "technical university of munich": "DEU",
    "helmholtz": "DEU", "robert koch": "DEU", "university of hamburg": "DEU",
    "delft": "NLD", "amsterdam": "NLD", "rotterdam": "NLD",
    "erasmus mc": "NLD", "wageningen": "NLD",
    "karolinska": "SWE", "sciencelife": "SWE", "lund university": "SWE",
    "copenhagen": "DNK", "dtu": "DNK", "statens serum": "DNK",
    "university of zurich": "CHE", "eth zurich": "CHE", "eawag": "CHE",
    "university of basel": "CHE", "epfl": "CHE",
    "university of milan": "ITA", "university of rome": "ITA", "padua": "ITA",
    "university of barcelona": "ESP", "csic": "ESP",
    "university of lisbon": "PRT", "porto": "PRT",
    "university of melbourne": "AUS", "university of sydney": "AUS",
    "csiro": "AUS", "university of queensland": "AUS", "monash": "AUS",
    "university of toronto": "CAN", "public health agency of canada": "CAN",
    "university of british columbia": "CAN", "waterloo": "CAN",
    "university of tokyo": "JPN", "nagoya university": "JPN",
    "osaka university": "JPN", "university of tohoku": "JPN",
    "seoul national university": "KOR", "korea university": "KOR",
    "yonsei": "KOR", "samsung": "KOR",
    "national university of singapore": "SGP", "ntu singapore": "SGP",
    "chinese academy of sciences": "CHN", "china cdc": "CHN",
    "chinese center for disease": "CHN", "peking university": "CHN",
    "tsinghua": "CHN", "fudan": "CHN", "wuhan university": "CHN",
    "bgi": "CHN", "shanghai": "CHN", "beijing": "CHN",
    "zhejiang university": "CHN", "sun yat-sen": "CHN", "guangzhou": "CHN",
    "hong kong": "HKG", "shenzhen": "CHN",
    "university of hong kong": "HKG", "chinese university of hong kong": "HKG",
    "taiwan cdc": "TWN", "national taiwan university": "TWN",
    "csir": "IND", "indian institute": "IND", "iit ": "IND",
    "national institute of virology": "IND", "pune": "IND", "delhi": "IND",
    "university of the witwatersrand": "ZAF", "wits": "ZAF",
    "national health laboratory service": "ZAF", "capetown": "ZAF",
    "UCT": "ZAF", "krisp": "ZAF", " Stellenbosch": "ZAF",
    "kehma": "ZAF", "nicd": "ZAF",
    "facultad de ciencias, udelar": "URY", "universidad de la republica": "URY",
    "tata institute": "IND", "christian medical college": "IND",
    "university of louisville": "USA", "cnrs": "FRA",
    "arizona state": "USA", "university of illinois": "USA",
    "central michigan": "USA", "university of manitoba": "CAN",
    "ginkgo": "USA", "vall d'hebron": "ESP", "university of arizona": "USA",
    "ethiopian public health institute": "ETH",
    "universidade de sao paulo": "BRA", "fiocruz": "BRA",
    "university of sao paulo": "BRA", "sao paulo": "BRA",
    "universidad de chile": "CHL", "ubuenosaires": "ARG",
    "universidad nacional autonoma de mexico": "MEX", "unam": "MEX",
    "university of lagos": "NGA", "nigeria cdc": "NGA", "iroko": "NGA",
    "university of nairobi": "KEN", "kemu": "KEN",
    "university of ghana": "GHA", "makerere": "UGA",
    "cairo university": "EGY", "ain shams": "EGY",
    "pasteur institute of iran": "IRN", "tehran": "IRN",
    "agakhan": "PAK", "karachi": "PAK",
    "icddr": "BGD", "dhaka": "BGD",
    "university of colombo": "LKA",
    "chulalongkorn": "THA", "mahidol": "THA",
    "vietnam academy": "VNM", "hanoi": "VNM",
    "universidad de buenos aires": "ARG", "cordoba": "ARG",
    "universidad de antioquia": "COL", "bogota": "COL",
    "quito": "ECU", "lima": "PER",
}
# 城市->国家 兜底词典 (机构名含这些城市名时)
CITY_MAP = {
    "london": "GBR", "cambridge, ma": "USA", "boston": "USA",
    "new york": "USA", "san diego": "USA", "san francisco": "USA",
    "seattle": "USA", "chicago": "USA", "atlanta": "USA", "houston": "USA",
    "paris": "FRA", "lyon": "FRA", "berlin": "DEU", "munich": "DEU",
    "hamburg": "DEU", "heidelberg": "DEU", "zurich": "CHE", "geneva": "CHE",
    "basel": "CHE", "milan": "ITA", "rome": "ITA", "madrid": "ESP",
    "barcelona": "ESP", "lisbon": "PRT", "stockholm": "SWE", "gothenburg": "SWE",
    "oslo": "NOR", "helsinki": "FIN", "copenhagen": "DNK", "dublin": "IRL",
    "amsterdam": "NLD", "utrecht": "NLD", "brussels": "BEL", "antwerp": "BEL",
    "vienna": "AUT", "prague": "CZE", "warsaw": "POL", "budapest": "HUN",
    "athens": "GRC", "istanbul": "TUR", "moscow": "RUS", "doha": "QAT",
    "dubai": "ARE", "riyadh": "SAU", "tel aviv": "ISR", "jerusalem": "ISR",
    "tokyo": "JPN", "osaka": "JPN", "kyoto": "JPN", "seoul": "KOR",
    "taipei": "TWN", "singapore": "SGP", "bangkok": "THA",
    "kuala lumpur": "MYS", "jakarta": "IDN", "manila": "PHL",
    "sydney": "AUS", "melbourne": "AUS", "brisbane": "AUS", "perth": "AUS",
    "auckland": "NZL", "toronto": "CAN", "vancouver": "CAN", "montreal": "CAN",
    "mexico city": "MEX", "bogota": "COL", "santiago": "CHL",
    "buenos aires": "ARG", "montevideo": "URY", "lima": "PER", "quito": "ECU",
    "johannesburg": "ZAF", "cape town": "ZAF", "pretoria": "ZAF",
    "lagos": "NGA", "accra": "GHA", "nairobi": "KEN", "kampala": "UGA",
    "addis ababa": "ETH", "cairo": "EGY", "casablanca": "MAR",
    "mumbai": "IND", "new delhi": "IND", "pune": "IND", "bangalore": "IND",
    "hyderabad": "IND", "karachi": "PAK", "lahore": "PAK", "dhaka": "BGD",
    "colombo": "LKA", "kathmandu": "NPL", "beijing": "CHN", "shanghai": "CHN",
    "wuhan": "CHN", "guangzhou": "CHN", "shenzhen": "CHN", "hangzhou": "CHN",
    "chengdu": "CHN", "nanjing": "CHN", "xiamen": "CHN", "harbin": "CHN",
}
UNKNOWN = "(unresolved)"

def resolve_country(owner):
    o = (owner or "").strip().lower()
    if not o:
        return UNKNOWN
    for k, v in PROXY_DB.items():
        if o == k or o.startswith(k + " ") or o.startswith(k + ","):
            return v
    for k, v in ORG_MAP.items():
        if k.lower() in o:
            return v
    for k, v in CITY_MAP.items():
        if k in o:
            return v
    return UNKNOWN

# ---------- 2. 采样国标准化 ----------
GEO_FIX = {
    "usa": "USA", "united states": "USA", "u.s.a": "USA", "us": "USA",
    "united kingdom": "GBR", "uk": "GBR", "england": "GBR", "scotland": "GBR",
    "wales": "GBR", "northern ireland": "GBR",
    "south korea": "KOR", "korea": "KOR", "republic of korea": "KOR",
    "hong kong": "HKG", "russia": "RUS", "russian federation": "RUS",
    "czech republic": "CZE", "vietnam": "VNM", "iran": "IRN",
    "the netherlands": "NLD", "netherlands": "NLD",
    "turkey": "TUR", "cote d'ivoire": "CIV", "ivory coast": "CIV",
    "chile": "CHL", "peru": "PER", "india": "IND", "china": "CHN",
    "taiwan": "TWN", "bolivia": "BOL", "venezuela": "VEN",
    "czechia": "CZE", "saudi arabia": "SAU", "south africa": "ZAF",
    "new zealand": "NZL", "puerto rico": "PRI", "guam": "GUM",
    "denmark": "DNK", "austria": "AUT", "switzerland": "CHE",
    "slovenia": "SVN", "norway": "NOR", "france": "FRA", "italy": "ITA",
    "germany": "DEU", "spain": "ESP", "canada": "CAN", "ethiopia": "ETH",
    "uruguay": "URY", "liechtenstein": "LIE", "india": "IND",
    "chile": "CHL", "peru": "PER",
}

def sample_country(geo):
    if not geo:
        return "(missing)"
    c = geo.split(":")[0].strip()
    return GEO_FIX.get(c.lower(), c)

# 收入组 (World Bank, 主要采样国; 未列出=unknown)
INCOME = {
    "USA": "HIC", "GBR": "HIC", "DNK": "HIC", "SWE": "HIC", "CHE": "HIC",
    "NLD": "HIC", "DEU": "HIC", "FRA": "HIC", "ITA": "HIC", "ESP": "HIC",
    "CAN": "HIC", "AUS": "HIC", "JPN": "HIC", "KOR": "HIC", "SGP": "HIC",
    "NZL": "HIC", "IRL": "HIC", "NOR": "HIC", "FIN": "HIC", "AUT": "HIC",
    "BEL": "HIC", "ISR": "HIC", "ARE": "HIC", "SAU": "HIC", "QAT": "HIC",
    "HKG": "HIC", "TWN": "HIC", "PRI": "HIC",
    "CHN": "UMIC", "BRA": "UMIC", "MEX": "UMIC", "ZAF": "UMIC",
    "ARG": "UMIC", "CHL": "UMIC", "URY": "UMIC", "COL": "UMIC",
    "PER": "UMIC", "ECU": "UMIC", "THA": "UMIC", "TUR": "UMIC",
    "RUS": "UMIC", "MAR": "UMIC", "BOL": "UMIC", "VEN": "UMIC",
    "NPL": "LMIC", "EGY": "LMIC", "MAR2": "LMIC", "LKA": "LMIC",
    "IND": "LMIC", "PAK": "LMIC", "BGD": "LMIC", "KEN": "LMIC",
    "NGA": "LMIC", "GHA": "LMIC", "CIV": "LMIC", "VNM": "LMIC",
    "IDN": "UMIC", "PHL": "LMIC", "MYS": "UMIC", "UKR": "LMIC",
    "ETH": "LIC", "UGA": "LIC", "MOZ": "LIC", "MWI": "LIC",
    "AFG": "LIC", "MLI": "LIC", "BFA": "LIC",
}

def parse_date(s):
    s = (s or "").strip()
    for fmt in ("%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%d", "%Y-%m", "%Y/%m/%d"):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            continue
    return None

def main():
    rows = list(csv.DictReader(open(MASTER, encoding="utf-8-sig")))
    print(f"Loaded {len(rows)} rows", flush=True)

    n = len(rows)
    aligned = misaligned = unresolved_owner = missing_geo = proxy_sub = 0
    pair_counter = Counter()
    country_rows = defaultdict(list)

    for r in rows:
        sc = sample_country(r["geo_location"])
        oc = resolve_country(r["owner"])
        r["_sample_country"] = sc
        r["_owner_country"] = oc
        if oc in ("EBI", "GBR") and r["owner"].strip().lower().startswith(("ebi", "ena")):
            proxy_sub += 1
        if r["owner"].strip().lower() in PROXY_DB:
            proxy_sub += 1
        if sc == "(missing)":
            missing_geo += 1
        elif oc == UNKNOWN:
            unresolved_owner += 1
        elif sc == oc:
            aligned += 1
        else:
            misaligned += 1
            pair_counter[(sc, oc)] += 1
        if sc != "(missing)":
            country_rows[sc].append(r)

    # 延迟计算
    lag_days = defaultdict(list)
    for r in rows:
        cd = parse_date(r["collection_date"])
        sd = parse_date(r["submission_date"])
        if cd and sd:
            lag_days[r["_sample_country"]].append((sd - cd).days)

    # --- v2: 集中度 + 修复的 proxy 计数 + 收入组 ---
    def gini(counter):
        vals = sorted(counter.values(), reverse=True)
        tot = sum(vals)
        if tot == 0:
            return 0.0
        cum = 0.0
        for i, v in enumerate(vals, 1):
            cum += i * v
        return (2 * cum) / (tot * len(vals)) - (len(vals) + 1) / len(vals)

    geo_counter = Counter(r["_sample_country"] for r in rows if r["_sample_country"] != "(missing)")
    hhi = sum((v / sum(geo_counter.values())) ** 2 for v in geo_counter.values())
    top3_share = sum(v for _, v in geo_counter.most_common(3)) / sum(geo_counter.values()) * 100

    PROXY_NAMES = {"ebi", "european bioinformatics institute", "ncbi", "ena", "ddbj"}
    proxy_fixed = sum(1 for r in rows if r["owner"].strip().lower() in PROXY_NAMES)
    ebi_n = sum(1 for r in rows if r["owner"].strip().lower() in
                {"ebi", "european bioinformatics institute"})

    income_dist = Counter(INCOME.get(k, "unknown") for k in geo_counter.elements())
    lic_umic = sum(v for k, v in geo_counter.items()
                   if INCOME.get(k, "unknown") in ("LMIC", "LIC", "UMIC"))

    stats_countries = []
    for sc, rs in sorted(country_rows.items(), key=lambda kv: -len(kv[1])):
        own_c = Counter(x["_owner_country"] for x in rs)
        mis = sum(v for k, v in own_c.items() if k not in (sc, UNKNOWN))
        lags = sorted(l for l in lag_days.get(sc, []) if -60 <= l <= 1460)
        stats_countries.append({
            "sample_country": sc,
            "income_group": INCOME.get(sc, "unknown"),
            "n_samples": len(rs),
            "n_distinct_owner_countries": len([k for k in own_c if k != UNKNOWN]),
            "misaligned_n": mis,
            "misalignment_rate_pct": round(mis / len(rs) * 100, 1),
            "foreign_submission_rate_pct": round(
                sum(v for k, v in own_c.items() if k not in (sc, UNKNOWN)) / len(rs) * 100, 1),
            "top3_owner_countries": own_c.most_common(3),
            "median_lag_days": (lags[len(lags)//2] if lags else None),
        })

    summary = {
        "total_samples": n,
        "aligned_same_country_pct": round(aligned / n * 100, 1),
        "misaligned_cross_country_pct": round(misaligned / n * 100, 1),
        "owner_unresolved_pct": round(unresolved_owner / n * 100, 1),
        "geo_missing_pct": round(missing_geo / n * 100, 1),
        "proxy_hub_submission_n": proxy_fixed,
        "proxy_hub_submission_pct": round(proxy_fixed / n * 100, 1),
        "ebi_single_hub_pct": round(ebi_n / n * 100, 1),
        "gini_country_concentration": round(gini(geo_counter), 3),
        "hhi": round(hhi, 3),
        "top3_country_share_pct": round(top3_share, 1),
        "n_sample_countries": len([k for k in geo_counter if k != "(missing)"]),
        "income_group_distribution": dict(income_dist),
        "lmic_umic_lic_total_samples": lic_umic,
        "lmic_umic_lic_pct": round(lic_umic / n * 100, 2),
        "top30_misaligned_pairs": [
            {"sampled": k[0], "submitted_by": k[1], "n": v}
            for k, v in pair_counter.most_common(30)
        ],
        "country_stats_top40": stats_countries[:40],
    }
    with open(os.path.join(DST, "analysis_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    with open(os.path.join(DST, "misalignment_matrix.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["sampled_country", "submitted_by_country", "n_samples"])
        w.writerows([[k[0], k[1], v] for k, v in pair_counter.most_common()])

    with open(os.path.join(DST, "country_stats.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(stats_countries[0].keys()) if stats_countries else ["x"])
        w.writeheader()
        w.writerows(stats_countries)

    print(json.dumps({k: summary[k] for k in [
        "total_samples", "aligned_same_country_pct", "misaligned_cross_country_pct",
        "owner_unresolved_pct", "n_sample_countries", "proxy_hub_submission_pct",
        "ebi_single_hub_pct", "gini_country_concentration", "top3_country_share_pct",
        "income_group_distribution", "lmic_umic_lic_pct"]},
        ensure_ascii=False, indent=2), flush=True)
    print("TOP10 misaligned pairs:", flush=True)
    for p in summary["top30_misaligned_pairs"][:10]:
        print(f"  {p['sampled']} -> {p['submitted_by']}: {p['n']}", flush=True)

if __name__ == "__main__":
    main()

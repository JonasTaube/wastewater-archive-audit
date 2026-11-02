# -*- coding: utf-8 -*-
"""多病原 + 宏基因组层扩展审计数据抓取 (v2.0)
Part A: 全量抓取 Influenza A/B, Norovirus, Rotavirus A, RSV 污水记录 (~1.4k条)
Part B: wastewater metagenome (373k条) 随机抽样 3000 条
输出: multipathogen_biosample.csv / metagenome_sample.csv / multipathogen_fetch.json
"""
import json, time, csv, os, re, random
import xml.etree.ElementTree as ET
import requests

DST = r"D:\南开CDC\task\1、监测殖民主义"
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
P = {"tool": "wbe_equity_audit_v2", "email": "apoc1920@gmail.com"}
FIELDS = ["accession", "owner", "owner_country", "submission_date",
          "geo_location", "collection_date", "organism", "isolation_source",
          "project", "sample_title", "pathogen_layer"]

TERMS = [
    ("influenza_a", "Influenza A virus[Organism] AND wastewater[All Fields]"),
    ("influenza_b", "Influenza B virus[Organism] AND wastewater[All Fields]"),
    ("norovirus", "Norovirus[Organism] AND wastewater[All Fields]"),
    ("rotavirus_a", "Rotavirus A[Organism] AND wastewater[All Fields]"),
    ("rsv", "Respiratory syncytial virus[Organism] AND wastewater[All Fields]"),
]
META_TERM = "wastewater metagenome[Organism]"
META_N_SAMPLE = 3000

def _retry(fn, *a, **kw):
    for attempt in range(5):
        try:
            return fn(*a, **kw)
        except (requests.RequestException, OSError) as e:
            wait = 2 ** (attempt + 1)
            print(f"  [retry {attempt+1}/5] {type(e).__name__}: {e}; sleep {wait}s", flush=True)
            time.sleep(wait)
    raise RuntimeError("retries exhausted")

def esearch_all(term, retmax=1000):
    """返回全部 id 列表(分页)。"""
    ids, start = [], 0
    while True:
        def _do():
            r = requests.get(f"{BASE}/esearch.fcgi",
                             params={**P, "db": "biosample", "term": term,
                                     "retmax": retmax, "retstart": start, "retmode": "json"},
                             timeout=30)
            r.raise_for_status()
            return r.json()["esearchresult"]
        j = _retry(_do)
        count, batch = int(j["count"]), j["idlist"]
        if not batch:
            break
        ids.extend(batch)
        start += len(batch)
        if start >= count:
            break
        time.sleep(0.4)
    return count, ids

def efetch(ids):
    def _do():
        r = requests.post(
            f"{BASE}/efetch.fcgi?db=biosample&retmode=xml&tool=wbe_equity_audit_v2&email=apoc1920@gmail.com",
            data={"id": ",".join(ids)}, timeout=90)
        r.raise_for_status()
        return r.text
    return _retry(_do)

def norm(s):
    return re.sub(r"[\s\-]+", "_", (s or "").strip().lower())

def parse(el, layer):
    d = dict.fromkeys(FIELDS, "")
    d["accession"] = el.get("accession", "")
    d["submission_date"] = el.get("submission_date", "")
    d["pathogen_layer"] = layer
    own = el.find("Owner")
    if own is not None:
        nm = own.find("Name")
        d["owner"] = (nm.text or "").strip() if nm is not None else (own.text or "").strip()
    og = el.find(".//Organism")
    if og is not None:
        d["organism"] = og.get("taxonomy_name", "")
    for attr in el.iter("Attribute"):
        hn = norm(attr.get("harmonized_name") or attr.get("attribute_name"))
        val = (attr.text or "").strip()
        if hn == "geo_loc_name":
            d["geo_location"] = val
        elif hn == "collection_date":
            d["collection_date"] = val
        elif hn == "isolation_source":
            d["isolation_source"] = val
    for link in el.iter("Link"):
        if link.get("target") == "bioproject":
            d["project"] = link.get("label", "")
    ttl = el.find(".//Title")
    if ttl is not None:
        d["sample_title"] = (ttl.text or "").strip()[:200]
    return d

def fetch_records(ids, layer, csv_path, mode="a"):
    got = 0
    for i in range(0, len(ids), 200):
        batch = ids[i:i+200]
        root = ET.fromstring(efetch(batch))
        rows = [parse(bs, layer) for bs in root.iter("BioSample")]
        with open(csv_path, mode, newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS)
            if mode == "a" and i == 0 and not os.path.getsize(csv_path):
                pass
            w.writerows(rows)
        got += len(rows)
        print(f"    [{layer}] {got}/{len(ids)}", flush=True)
        time.sleep(0.4)
    return got

def main():
    mp_csv = os.path.join(DST, "multipathogen_biosample.csv")
    mg_csv = os.path.join(DST, "metagenome_sample.csv")
    log = {"pathogens": {}, "metagenome": {}, "ascertained": "2026-11-02"}
    # ---------- Part A: 病原层全量 ----------
    with open(mp_csv, "w", newline="", encoding="utf-8-sig") as f:
        csv.DictWriter(f, fieldnames=FIELDS).writeheader()
    for layer, term in TERMS:
        count, ids = esearch_all(term)
        print(f"  {layer}: {count} records", flush=True)
        n = fetch_records(ids, layer, mp_csv, mode="a") if ids else 0
        log["pathogens"][layer] = {"count": count, "fetched": n}
    # ---------- Part B: 宏基因组层抽样 ----------
    count, ids = esearch_all(META_TERM)
    print(f"  wastewater metagenome: {count} total", flush=True)
    random.seed(42)
    sample_ids = random.sample(ids, min(META_N_SAMPLE, len(ids)))
    with open(mg_csv, "w", newline="", encoding="utf-8-sig") as f:
        csv.DictWriter(f, fieldnames=FIELDS).writeheader()
    n = fetch_records(sample_ids, "wastewater_metagenome", mg_csv, mode="a")
    log["metagenome"] = {"total_population": count, "sample_n": len(sample_ids), "fetched": n}
    with open(os.path.join(DST, "multipathogen_fetch.json"), "w", encoding="utf-8") as f:
        json.dump(log, f, ensure_ascii=False, indent=2)
    print("DONE", flush=True)
    print(json.dumps(log, indent=2), flush=True)

if __name__ == "__main__":
    main()

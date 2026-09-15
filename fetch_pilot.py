# -*- coding: utf-8 -*-
"""污水测序档案数据提取审计 - 可行性验证 v2: 修复字段解析"""
import json, time, csv, os, sys, re
import xml.etree.ElementTree as ET

try:
    import requests
except ImportError:
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "--quiet", "install", "requests"])
    import requests

OUT_DIR = r"C:\Users\lenovo\WorkBuddy\2026-09-15-12-56-54\outputs\wastewater_extraction"
os.makedirs(OUT_DIR, exist_ok=True)
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
P = {"tool": "wbe_equity_audit", "email": "apoc1920@gmail.com"}

def norm(s):
    return re.sub(r"[\s\-]+", "_", (s or "").strip().lower())

def esearch(term, retmax=1000):
    r = requests.get(f"{BASE}/esearch.fcgi",
                     params={**P, "db": "biosample", "term": term,
                             "retmax": retmax, "retstart": 0, "retmode": "json"}, timeout=30)
    r.raise_for_status()
    j = r.json()["esearchresult"]
    return int(j["count"]), j["idlist"]

def efetch(ids):
    r = requests.post(f"{BASE}/efetch.fcgi?db=biosample&retmode=xml&tool=wbe_equity_audit&email=apoc1920@gmail.com",
                      data={"id": ",".join(ids)}, timeout=60)
    r.raise_for_status()
    return r.text

def parse(el):
    d = {"accession": el.get("accession", ""),
         "owner": "", "submission_date": el.get("submission_date", ""),
         "geo_location": "", "collection_date": "", "organism": "",
         "isolation_source": "", "project": ""}
    own = el.find("Owner")
    if own is not None:
        nm = own.find("Name")
        if nm is not None and (nm.text or "").strip():
            d["owner"] = nm.text.strip()
        else:
            d["owner"] = (own.text or "").strip()
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
    return d

def country_of(geo):
    if not geo:
        return "(missing)"
    return geo.split(":")[0].strip()

def main():
    term = "SARS-CoV-2[Organism] AND wastewater[All Fields]"
    count, ids = esearch(term, retmax=1000)
    print(f"TOTAL HITS: {count}; parsing first {len(ids)}", flush=True)

    rows = []
    for i in range(0, len(ids), 200):
        root = ET.fromstring(efetch(ids[i:i+200]))
        for bs in root.iter("BioSample"):
            rows.append(parse(bs))
        time.sleep(0.4)
    print(f"Parsed: {len(rows)}", flush=True)

    csv_path = os.path.join(OUT_DIR, "biosample_window1k.csv")
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    from collections import Counter
    geo_c = Counter(country_of(r["geo_location"]) for r in rows)
    n = len(rows)
    summary = {
        "total_hits_all_time": count,
        "window_parsed": n,
        "n_distinct_countries": len([k for k in geo_c if k != "(missing)"]),
        "geo_missing_pct": round(geo_c.get("(missing)", 0) / n * 100, 1),
        "field_coverage_pct": {
            "collection_date": round(sum(1 for r in rows if r["collection_date"]) / n * 100, 1),
            "submission_date": round(sum(1 for r in rows if r["submission_date"]) / n * 100, 1),
            "owner_institution": round(sum(1 for r in rows if r["owner"]) / n * 100, 1),
        },
        "geo_country_top20": geo_c.most_common(20),
        "top10_owners": Counter((r["owner"][:55] if r["owner"] else "(missing)") for r in rows).most_common(10),
    }
    with open(os.path.join(OUT_DIR, "summary_stats.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)

if __name__ == "__main__":
    main()

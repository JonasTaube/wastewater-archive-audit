# -*- coding: utf-8 -*-
"""全量抓取 NCBI BioSample: SARS-CoV-2 x wastewater (82,757 条)
- 断点续传: progress.json 记录 next_start, 重跑自动续
- 增量落盘: 每批 200 条立即追加 master CSV
- 限速: 无 API key, 3 req/s (每请求 sleep 0.4s)
"""
import json, time, csv, os, re, sys
import xml.etree.ElementTree as ET
import requests

DST = r"D:\南开CDC\task\1、监测殖民主义"
MASTER = os.path.join(DST, "master_biosample.csv")
PROGRESS = os.path.join(DST, "progress.json")
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
P = {"tool": "wbe_equity_audit", "email": "apoc1920@gmail.com"}
FIELDS = ["accession", "owner", "owner_country", "submission_date",
          "geo_location", "collection_date", "organism", "isolation_source",
          "project", "sample_title"]

def norm(s):
    return re.sub(r"[\s\-]+", "_", (s or "").strip().lower())

def _retry(fn, *a, **kw):
    """带指数退避的重试: 最多5次, 间隔 2/4/8/16/32s"""
    for attempt in range(5):
        try:
            return fn(*a, **kw)
        except (requests.RequestException, OSError) as e:
            wait = 2 ** (attempt + 1)
            print(f"  [retry {attempt+1}/5] {type(e).__name__}: {e}; sleep {wait}s", flush=True)
            time.sleep(wait)
    raise RuntimeError("retries exhausted")

def esearch_page(start, retmax=1000):
    def _do():
        r = requests.get(f"{BASE}/esearch.fcgi",
                         params={**P, "db": "biosample",
                                 "term": "SARS-CoV-2[Organism] AND wastewater[All Fields]",
                                 "retmax": retmax, "retstart": start, "retmode": "json"}, timeout=30)
        r.raise_for_status()
        return r.json()["esearchresult"]
    j = _retry(_do)
    return int(j["count"]), j["idlist"]

def efetch(ids):
    def _do():
        r = requests.post(
            f"{BASE}/efetch.fcgi?db=biosample&retmode=xml&tool=wbe_equity_audit&email=apoc1920@gmail.com",
            data={"id": ",".join(ids)}, timeout=60)
        r.raise_for_status()
        return r.text
    return _retry(_do)

def parse(el):
    d = dict.fromkeys(FIELDS, "")
    d["accession"] = el.get("accession", "")
    d["submission_date"] = el.get("submission_date", "")
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

def load_progress():
    if os.path.exists(PROGRESS):
        with open(PROGRESS, encoding="utf-8") as f:
            return json.load(f)
    return {"next_start": 0, "total": None}

def save_progress(p):
    with open(PROGRESS, "w", encoding="utf-8") as f:
        json.dump(p, f)

def main():
    prog = load_progress()
    start = prog["next_start"]
    if not os.path.exists(MASTER):
        with open(MASTER, "w", newline="", encoding="utf-8-sig") as f:
            csv.DictWriter(f, fieldnames=FIELDS).writeheader()

    while True:
        count, ids = esearch_page(start, 1000)
        if prog["total"] is None:
            prog["total"] = count
            print(f"TOTAL: {count}", flush=True)
        if not ids:
            break
        for i in range(0, len(ids), 200):
            batch = ids[i:i+200]
            root = ET.fromstring(efetch(batch))
            rows = [parse(bs) for bs in root.iter("BioSample")]
            with open(MASTER, "a", newline="", encoding="utf-8-sig") as f:
                w = csv.DictWriter(f, fieldnames=FIELDS)
                w.writerows(rows)
            start += len(batch)
            prog["next_start"] = start
            save_progress(prog)
            print(f"progress: {start}/{prog['total']}", flush=True)
            time.sleep(0.4)
        if len(ids) < 1000 or start >= count:
            break
    print("DONE", flush=True)

if __name__ == "__main__":
    main()

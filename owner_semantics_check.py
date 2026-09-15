# -*- coding: utf-8 -*-
"""Owner-semantics verification: sample 300 EBI-owner (SAMEA) records, fetch ENA
center_name for each, classify whether the record's registration centre is EBI
itself or the true national submitting centre."""
import csv, json, io, sys, time, random, re
from collections import Counter

import requests, warnings
warnings.filterwarnings("ignore")

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
DST = r"D:\南开CDC\task\1、监测殖民主义"

rows = list(csv.DictReader(open(DST + r"\master_biosample.csv", encoding="utf-8-sig")))
ebi = [r for r in rows if r["owner"] in ("EBI", "European Bioinformatics Institute")]
random.seed(42)
sample = random.sample(ebi, 300)

s = requests.Session(); s.trust_env = False
out, fails = [], 0
for i, r in enumerate(sample):
    acc = r["accession"]
    try:
        resp = s.get(f"https://www.ebi.ac.uk/ena/browser/api/xml/{acc}", timeout=30)
        m = re.search(r'center_name="([^"]*)"', resp.text)
        cn = m.group(1) if m else "NOT_FOUND"
    except Exception:
        cn = "FETCH_FAIL"; fails += 1
    out.append(dict(acc=acc, country=(r["geo_location"] or "").split(":")[0], center_name=cn))
    if (i + 1) % 50 == 0:
        print(f"  {i+1}/300 done")
    time.sleep(0.12)

cnt = Counter(o["center_name"] for o in out)
print("\ncenter_name distribution (n=300):")
for k, v in cnt.most_common(30):
    print(f"  {k!r}: {v}")
print("fetch fails:", fails)

# classification
ebi_itself = cnt.get("EBI", 0) + cnt.get("European Bioinformatics Institute", 0) + cnt.get("EBIRNA", 0)
print("\ncenter_name == EBI itself:", ebi_itself, f"({100*ebi_itself/300:.1f}%)")
print("center_name = true national submitting centre:", 300 - ebi_itself - fails,
      f"({100*(300-ebi_itself-fails)/300:.1f}%)")

json.dump(dict(n=len(out), fails=fails, distribution=cnt.most_common(),
               ebi_itself=ebi_itself, detail=out),
          open(DST + r"\owner_semantics_check.json", "w"), indent=1)
print("saved owner_semantics_check.json")

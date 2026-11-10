# -*- coding: utf-8 -*-
"""Standalone AFT fit on no-window frame; prints full summary index."""
import csv, sys, json, io, os
from datetime import datetime
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.path.insert(0, "材料3_统计建模代码及运行结果/code")
from map_analyze import parse_date, sample_country
import numpy as np, pandas as pd

HUB = ("EBI", "European Bioinformatics Institute")
OWNER_CTRY = {"NIBSC": "GBR", "Tata Institute for Genetics and Society": "IND",
              "CSIR-National Chemical Laboratory, Pune": "IND",
              "Christian Medical College and Hospital": "IND",
              "University of Louisville": "USA", "Arizona State University": "USA",
              "University of Wisconsin": "USA", "Broad Institute of MIT and Harvard": "USA",
              "University of Illinois Urbana-Champaign": "USA",
              "Central Michigan University College of Medicine": "USA",
              "The University of Arizona": "USA", "Ginkgo Bioworks": "USA",
              "CNRS": "FRA", "Ethiopian Public Health Institute": "ETH",
              "Facultad de Ciencias, UdelaR": "URY", "University of Manitoba": "CAN",
              "Public Health Agency of Canada": "CAN",
              "Universitat de Barcelona": "ESP",
              "Vall d'Hebron Institut de Recerca": "ESP"}
rows = []
for r in csv.DictReader(open("材料1_归档审计数据集及数据字典/data/master_biosample.csv",
                             encoding="utf-8-sig")):
    d1, d2 = parse_date(r["submission_date"]), parse_date(r["collection_date"])
    c = sample_country(r["geo_location"]); c = "SWE" if c == "Sweden" else c
    if r["owner"] in HUB: pool = "hub_routed"
    elif OWNER_CTRY.get(r["owner"], "?") == c and c != "(missing)":
        pool = "indep_LMIC" if c in ("IND", "ETH") else "indep_HIC_dom"
    else: continue
    if not (d1 and d2): continue
    lag = (d1 - d2).days
    if lag <= -60: continue
    rows.append((max(lag, 0.5), 1, pool, d2.year))
f = pd.DataFrame(rows, columns=["dur", "event", "pool", "coll_year"])
f["pool_hub"] = (f.pool == "hub_routed").astype(float)
f["pool_lmic"] = (f.pool == "indep_LMIC").astype(float)
f["coll_year_c"] = (f.coll_year - 2021).astype(float)
print("n:", len(f), flush=True)

from lifelines import WeibullAFTFitter, LogNormalAFTFitter
out = {}
for name, Fitter in [("weibull", WeibullAFTFitter), ("lognormal", LogNormalAFTFitter)]:
    try:
        fit = Fitter()
        fit.fit(f[["dur", "event", "pool_hub", "pool_lmic", "coll_year_c"]],
                duration_col="dur", event_col="event")
        smry = fit.summary
        print(name, "INDEX:", list(smry.index), flush=True)
        keep = {}
        for k in smry.index:
            if "pool_" in str(k):
                keep[str(k)] = {"coef": float(smry.loc[k, "coef"]),
                                "exp": float(np.exp(smry.loc[k, "coef"])),
                                "lo": float(np.exp(smry.loc[k, "coef lower 95%"])),
                                "hi": float(np.exp(smry.loc[k, "coef upper 95%"])),
                                "p": float(smry.loc[k, "p"])}
        out[name] = keep
    except Exception as e:
        import traceback; traceback.print_exc()
        out[name] = {"error": str(e)[:300]}
json.dump(out, open("aft_results.json", "w"), indent=1)
print(json.dumps(out, indent=1))

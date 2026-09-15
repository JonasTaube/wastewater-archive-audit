# -*- coding: utf-8 -*-
"""Round-2 audit: flow reconciliation, WHO mortality ranking, all-sample burden,
rho re-check, back-fill sensitivity. Outputs JSON + reconciliation table."""
import csv, json, io, sys
from collections import defaultdict
from datetime import datetime

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.path.insert(0, r"D:\南开CDC\task\1、监测殖民主义")
import pandas as pd
import scipy.stats as st
from map_analyze import parse_date

DST = r"D:\南开CDC\task\1、监测殖民主义"

rows = list(csv.DictReader(open(DST + r"\master_biosample.csv", encoding="utf-8-sig")))
print("total records:", len(rows))

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
HUB_OWNERS = ("EBI", "European Bioinformatics Institute")

def sc(geo):
    from map_analyze import sample_country
    c = sample_country(geo)
    return "SWE" if c == "Sweden" else c

# ---------------- flow reconciliation ----------------
flow = defaultdict(int)
for r in rows:
    flow["total"] += 1
    c = sc(r["geo_location"])
    is_hub = r["owner"] in HUB_OWNERS
    d1, d2 = parse_date(r["submission_date"]), parse_date(r["collection_date"])
    valid = bool(d1 and d2) and (-60 < (d1 - d2).days <= 1460)
    if not valid:
        flow["invalid_dates"] += 1
        flow["inv_" + ("hub" if is_hub else "ind")] += 1
        continue
    if is_hub:
        flow["valid_hub"] += 1
    else:
        oc = OWNER_CTRY.get(r["owner"], "?")
        if oc == c:
            flow["valid_ind_dom"] += 1
            if c in ("IND", "ETH"):
                flow["valid_ind_lmic"] += 1
            flow["dom_by_" + c] += 1
        elif oc == "?":
            flow["valid_ind_unknown_owner"] += 1
        else:
            flow["valid_ind_foreign_owner"] += 1

hub_total = sum(1 for r in rows if r["owner"] in HUB_OWNERS)
print(json.dumps(flow, indent=1))
print("hub total (owner=EBI):", hub_total, "| valid hub:", flow["valid_hub"],
      "| hub invalid-date:", flow["inv_hub"])
print("union of pools (valid only):", flow["valid_hub"] + flow["valid_ind_dom"],
      "| unassigned:", flow["total"] - flow["valid_hub"] - flow["valid_ind_dom"])
print("independent-domestic by country:", {k: v for k, v in sorted(flow.items()) if k.startswith("dom_by_")})

# independent domestic by country from full data (no date filter)
dom_all = defaultdict(int)
for r in rows:
    c = sc(r["geo_location"])
    if r["owner"] not in HUB_OWNERS and OWNER_CTRY.get(r["owner"], "?") == c:
        dom_all[c] += 1
print("independent domestic ALL (no date filter):", dict(sorted(dom_all.items(), key=lambda kv: -kv[1])),
      "sum:", sum(dom_all.values()))

# ---------------- WHO mortality ranking ----------------
who = pd.read_csv(DST + r"\who_covid_global.csv", encoding="utf-8-sig")
latest = who.sort_values("Date_reported").groupby("Country_code").tail(1)
lat = latest[["Country", "Country_code", "Cumulative_deaths"]].copy()
lat = lat.sort_values("Cumulative_deaths", ascending=False).reset_index(drop=True)
lat["rank"] = lat.index + 1
print("\nWHO top-12 mortality:")
print(lat.head(12).to_string())
world_total = lat["Cumulative_deaths"].sum()
print("WHO sum of national cumulative deaths:", world_total)

ISO2 = {"GBR": "GB", "NLD": "NL", "DNK": "DK", "AUT": "AT", "CHE": "CH", "SVN": "SI",
        "USA": "US", "SWE": "SE", "NOR": "NO", "FRA": "FR", "LIE": "LI", "ITA": "IT",
        "DEU": "DE", "CHL": "CL", "CAN": "CA", "PER": "PE", "IND": "IN", "ETH": "ET",
        "URY": "UY", "ESP": "ES"}
REC = {"GBR": 49646, "NLD": 17405, "DNK": 8108, "AUT": 3497, "CHE": 1765, "SVN": 616,
       "USA": 563, "IND": 311, "SWE": 220, "NOR": 214, "FRA": 147, "ETH": 79,
       "LIE": 31, "ITA": 24, "DEU": 22, "URY": 10, "CHL": 9, "ESP": 3, "CAN": 2, "PER": 1}

# join: all WHO entities with record counts (0 if absent)
rec_by_cc2 = {ISO2[k]: v for k, v in REC.items()}
lat["records"] = lat["Country_code"].map(rec_by_cc2).fillna(0).astype(int)
n_entities = len(lat)
print("WHO entities:", n_entities)
matched = lat[lat["records"] > 0]
print("entities with records:", len(matched))

# all-sample Spearman (239 entities incl. zeros)
rho_all, p_all = st.spearmanr(lat["Cumulative_deaths"], lat["records"])
print(f"\nALL-SAMPLE Spearman (n={n_entities}): rho={rho_all:.4f} p={p_all:.4g}")

# among-contributors re-check
lat20 = lat[lat["records"] > 0].copy()
rho20, p20 = st.spearmanr(lat20["Cumulative_deaths"], lat20["records"])
print(f"Contributors-only Spearman (n={len(lat20)}): rho={rho20:.4f} p={p20:.4g}")

# HIC subset re-check
HIC = {"GB", "NL", "DK", "AT", "CH", "SI", "US", "SE", "NO", "FR", "LI", "IT", "DE",
       "UY", "CL", "ES", "CA"}
lat17 = lat20[lat20["Country_code"].isin(HIC)]
rho17, p17 = st.spearmanr(lat17["Cumulative_deaths"], lat17["records"])
print(f"HIC subset Spearman (n={len(lat17)}): rho={rho17:.4f} p={p17:.4g}")

# zero-record high-mortality countries
zero_hi = lat[(lat["records"] == 0) & (lat["Cumulative_deaths"] > 50000)]
print("\nZero-record countries with >50k deaths:")
print(zero_hi[["rank", "Country", "Cumulative_deaths"]].to_string())

# Peru rank
peru = lat[lat["Country_code"] == "PE"]
print("\nPeru rank (WHO cumulative reported deaths):", int(peru["rank"].iloc[0]),
      "deaths:", int(peru["Cumulative_deaths"].iloc[0]))

# presence vs deaths (logistic-ish: rank-biserial of deaths by presence)
g1 = lat[lat["records"] > 0]["Cumulative_deaths"]; g0 = lat[lat["records"] == 0]["Cumulative_deaths"]
u, pu = st.mannwhitneyu(g1, g0, alternative="two-sided")
print(f"Mann-Whitney presence vs deaths: n1={len(g1)} n0={len(g0)} U={u} p={pu:.3g}")

# ---------------- lag sensitivity ----------------
lags = []
for r in rows:
    d1, d2 = parse_date(r["submission_date"]), parse_date(r["collection_date"])
    if not d1 or not d2:
        continue
    lag = (d1 - d2).days
    if not (-60 < lag <= 1460):
        continue
    c = sc(r["geo_location"])
    is_hub = r["owner"] in HUB_OWNERS
    oc = OWNER_CTRY.get(r["owner"], "?")
    pool = "hub" if is_hub else ("dom" if oc == c else "other")
    coll_year = d2.year
    lags.append(dict(lag=lag, pool=pool, c=c, coll_year=coll_year, wave=r["submission_date"][:7]))

df = pd.DataFrame(lags)
def pool_stats(d, label):
    h = d[d["pool"] == "hub"]["lag"]; i = d[d["pool"] == "dom"]["lag"]
    H, p = st.kruskal(i, h)
    print(f"{label}: hub n={len(h)} med={h.median():.0f} | ind-dom n={len(i)} med={i.median():.0f} "
          f"| H={H:,.1f} p={p:.3g}")
    return dict(label=label, hub_n=len(h), hub_med=float(h.median()),
                dom_n=len(i), dom_med=float(i.median()), H=H, p=p)

sens = []
sens.append(pool_stats(df, "full window"))
sens.append(pool_stats(df[df["coll_year"] >= 2022], "collected >= 2022"))
sens.append(pool_stats(df[df["lag"] <= 365], "lag <= 365d only"))
sens.append(pool_stats(df[df["coll_year"] >= 2023], "collected >= 2023"))

# hub lag by collection year (back-fill signature)
hub = df[df["pool"] == "hub"]
by_year = hub.groupby("coll_year")["lag"].agg(["count", "median"])
print("\nhub lag by collection year:")
print(by_year.to_string())

# wave back-fill check: three bulk waves -> collection years
wq = {"2022Q3": ("2022-07", "2022-09"), "2024Q3": ("2024-07", "2024-09"), "2025Q4": ("2025-10", "2025-12")}
hubsub = hub[hub["lag"] > 730]
print("\nhub records with lag>730d:", len(hubsub), "of", len(hub),
      f"({100*len(hubsub)/len(hub):.1f}%)")

json_out = dict(flow=dict(flow), hub_total=hub_total,
                dom_all=dict(dom_all), dom_all_sum=sum(dom_all.values()),
                who_rank_top12=lat.head(12)[["rank", "Country", "Country_code", "Cumulative_deaths", "records"]].to_dict("records"),
                n_entities=int(n_entities),
                all_sample=dict(n=int(n_entities), rho=float(rho_all), p=float(p_all)),
                contributors=dict(n=len(lat20), rho=float(rho20), p=float(p20)),
                hic=dict(n=len(lat17), rho=float(rho17), p=float(p17)),
                peru_rank=int(peru["rank"].iloc[0]), peru_deaths=int(peru["Cumulative_deaths"].iloc[0]),
                zero_hi=zero_hi[["rank", "Country", "Cumulative_deaths"]].to_dict("records"),
                mannwhitney=dict(U=float(u), p=float(pu), n1=len(g1), n0=len(g0)),
                sensitivity=sens,
                hub_by_year={int(k): [int(v["count"]), float(v["median"])] for k, v in by_year.iterrows()})
json.dump(json_out, open(DST + r"\round2_analysis.json", "w"), indent=1)
print("\nsaved round2_analysis.json")

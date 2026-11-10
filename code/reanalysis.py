# -*- coding: utf-8 -*-
"""Reviewer-response reanalysis for Sci Data wastewater audit manuscript.

Covers:
  M1  administrative censoring framework -> recomputed Cox HR
  M2  lag-window (1,460 d) exclusions + extended-window (2,920 d) sensitivity
  M3  Firth penalised logistic + drop-Peru sensitivity + full income ORs
  M4  AFT (Weibull / log-normal) + RMST differences
  M5  exact record accounting (flow-diagram numbers)
  M8  metagenome extrapolation arithmetic check
  minors: negative-lag counts, strata 2/3 window counts, cluster bootstrap, seed
"""
import csv, json, os, sys, io
from collections import Counter
from datetime import datetime

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(BASE, "材料3_统计建模代码及运行结果", "code"))
DATA = os.path.join(BASE, "材料1_归档审计数据集及数据字典", "data")

import numpy as np
import pandas as pd
from map_analyze import parse_date, sample_country

ASC_SARS = datetime(2026, 9, 15)   # SARS-CoV-2 ascertainment date
ASC_MULTI = datetime(2026, 11, 2)  # strata 2/3 ascertainment date

HUB_OWNERS = ("EBI", "European Bioinformatics Institute")
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
HIC_COUNTRIES = {"USA", "GBR", "DNK", "SWE", "CHE", "NLD", "DEU", "FRA", "ITA",
                 "ESP", "CAN", "AUS", "JPN", "KOR", "SGP", "NZL", "IRL", "NOR",
                 "FIN", "AUT", "BEL", "ISR", "ARE", "SAU", "QAT", "HKG", "TWN",
                 "PRI", "URY", "CHL", "LIE"}

def sc(geo):
    c = sample_country(geo)
    return "SWE" if c == "Sweden" else c

def build_records(path, asc_date):
    """Replicates upgrade_stats.build_lag_df pool logic, but keeps every record
    with full audit fields for accounting."""
    recs = []
    for r in csv.DictReader(open(path, encoding="utf-8-sig")):
        d1, d2 = parse_date(r["submission_date"]), parse_date(r["collection_date"])
        c = sc(r["geo_location"])
        is_hub = r["owner"] in HUB_OWNERS
        oc = OWNER_CTRY.get(r["owner"], "?")
        if is_hub:
            pool = "hub_routed"
        elif oc == c and c != "(missing)":
            pool = "indep_LMIC" if c in ("IND", "ETH") else "indep_HIC_dom"
        else:
            pool = "foreign_indep"
        recs.append(dict(accession=r["accession"], owner=r["owner"],
                         sample_country=c, pool=pool,
                         sd=d1, cd=d2,
                         lag=(d1 - d2).days if (d1 and d2) else None,
                         coll_year=d2.year if d2 else None))
    return pd.DataFrame(recs)

# ---------------- Part 1: SARS stratum accounting ----------------
print("=" * 60); print("PART 1: SARS stratum accounting (M5/M2)")
df = build_records(os.path.join(DATA, "master_biosample.csv"), ASC_SARS)
N = len(df)
acc = {"total": N}
acc["pool_all"] = df["pool"].value_counts().to_dict()
acc["no_geo"] = int((df["sample_country"] == "(missing)").sum())
acc["no_geo_by_pool"] = df[df["sample_country"] == "(missing)"]["pool"].value_counts().to_dict()

# ownership split
acc["hub_owned"] = int((df["pool"] == "hub_routed").sum())
acc["independent"] = int(df["pool"].isin(["indep_HIC_dom", "indep_LMIC", "foreign_indep"]).sum())

# foreign attribution (owner country != sampling country), replicated from map_analyze
rows_raw = list(csv.DictReader(open(os.path.join(DATA, "master_biosample.csv"), encoding="utf-8-sig")))
foreign = 0; pair = Counter()
from map_analyze import resolve_country, UNKNOWN
for r in rows_raw:
    s = sample_country(r["geo_location"]); o = resolve_country(r["owner"])
    if s == "(missing)" or o == UNKNOWN:
        continue
    if s != o:
        foreign += 1; pair[(s, o)] += 1
acc["foreign_attributed"] = foreign
acc["foreign_to_GBR"] = int(sum(v for k, v in pair.items() if k[1] == "GBR"))

# date validity
acc["valid_pair"] = int(df["lag"].notna().sum())
acc["invalid_pair"] = int(df["lag"].isna().sum())
inv = df[df["lag"].isna()]
inv_by_pool = inv["pool"].value_counts().to_dict()
acc["invalid_pair_by_pool"] = inv_by_pool
# which date missing (parse_date returns None; empty string check)
miss_sd = miss_cd = both = 0
for r in rows_raw:
    s_raw = (r["submission_date"] or "").strip(); c_raw = (r["collection_date"] or "").strip()
    d1, d2 = parse_date(r["submission_date"]), parse_date(r["collection_date"])
    if not d1 and not d2: both += 1
    elif not d1: miss_sd += 1
    elif not d2: miss_cd += 1
acc["invalid_missing_both"] = both
acc["invalid_missing_submission"] = miss_sd
acc["invalid_missing_collection"] = miss_cd

# window exclusions in ORIGINAL (−60, 1460] analysis
lag_ok = df[df["lag"].notna()]
acc["lag_gt_1460"] = int((lag_ok["lag"] > 1460).sum())
acc["lag_le_minus60"] = int((lag_ok["lag"] <= -60).sum())
acc["lag_neg_retained"] = int(((lag_ok["lag"] > -60) & (lag_ok["lag"] < 0)).sum())
acc["max_lag"] = int(lag_ok["lag"].max())
acc["min_lag"] = int(lag_ok["lag"].min())
# original velocity set check
orig = lag_ok[(lag_ok["lag"] > -60) & (lag_ok["lag"] <= 1460) & (lag_ok["pool"] != "foreign_indep")]
acc["orig_velocity_n"] = int(len(orig))
acc["orig_velocity_pools"] = orig["pool"].value_counts().to_dict()
fi = df[df["pool"] == "foreign_indep"]
acc["foreign_indep_detail"] = fi[["accession", "sample_country", "lag", "sd", "cd"]].to_dict("records")
# reconcile 3,392 = invalid pair + window exclusions
acc["reconcile"] = {"valid_pair": acc["valid_pair"], "invalid_pair": acc["invalid_pair"],
                    "window_excluded": acc["lag_gt_1460"] + acc["lag_le_minus60"],
                    "check": acc["valid_pair"] - 1 - (acc["lag_gt_1460"] + acc["lag_le_minus60"])}
print("ACCOUNTING:", json.dumps(acc, indent=1, default=str))

# ---------------- Part 2: M1/M2 Cox refits ----------------
print("=" * 60); print("PART 2: M1/M2 Cox refits")
from lifelines import CoxPHFitter
from lifelines.utils import restricted_mean_survival_time
from lifelines import KaplanMeierFitter

# Frames:
#  (A) no-window frame: every record with valid collection date is an EVENT at
#      observed lag (max 1,677 d); the old (−60,1460] truncation is dropped.
#  (B) complete-follow-up frame: collections <= 2022-09-15, so every record has
#      >= 1,460 days of potential follow-up before ascertainment (closed cohort).
#  (C) original windowed frame (−60,1460] for reproduction.

def make_frame(min_coll=None, max_coll=None, lo=-60, hi=None):
    recs = []
    for r in rows_raw:
        d1, d2 = parse_date(r["submission_date"]), parse_date(r["collection_date"])
        c = sc(r["geo_location"]); is_hub = r["owner"] in HUB_OWNERS
        oc = OWNER_CTRY.get(r["owner"], "?")
        if is_hub: pool = "hub_routed"
        elif oc == c and c != "(missing)":
            pool = "indep_LMIC" if c in ("IND", "ETH") else "indep_HIC_dom"
        else:
            pool = "foreign_indep"
        if pool == "foreign_indep":
            continue
        if not (d1 and d2):
            continue
        lag = (d1 - d2).days
        if lag <= lo:
            continue
        if hi is not None and lag > hi:
            continue
        if min_coll and d2 < min_coll:
            continue
        if max_coll and d2 > max_coll:
            continue
        recs.append(dict(dur=max(lag, 0.5), event=1, pool=pool,
                         coll_year=d2.year, country=c,
                         pool_hub=int(pool == "hub_routed"),
                         pool_lmic=int(pool == "indep_LMIC")))
    f = pd.DataFrame(recs)
    f["coll_year_c"] = f["coll_year"] - 2021
    return f

def fit_cox(f, label):
    out = {"n": int(len(f)), "pools": f["pool"].value_counts().to_dict()}
    cph = CoxPHFitter()
    cph.fit(f[["dur", "event", "pool_hub", "pool_lmic", "coll_year_c"]],
            duration_col="dur", event_col="event", fit_options={"step_size": 0.5})
    s = cph.summary
    out["cox"] = {k: {"HR": float(s.loc[k, "exp(coef)"]),
                      "lo": float(s.loc[k, "exp(coef) lower 95%"]),
                      "hi": float(s.loc[k, "exp(coef) upper 95%"]),
                      "p": float(s.loc[k, "p"])} for k in s.index}
    print(label, json.dumps(out, indent=1))
    return out, cph

COX_CUT = datetime(2022, 9, 15)  # ascertainment − 1,460 days
fA = make_frame()                       # no window
fB = make_frame(min_coll=None, max_coll=COX_CUT)  # complete follow-up
fC = make_frame(hi=1460)                # original window
resA, cphA = fit_cox(fA, "A no-window:")
resB, cphB = fit_cox(fB, "B complete-follow-up (<=2022-09-15):")
resC, cphC = fit_cox(fC, "C original window:")
m1 = {"A_nowindow": resA, "B_complete_followup": resB, "C_original_window": resC,
      "coll_range_A": [str(fA["coll_year"].min()), str(fA["coll_year"].max())],
      "coll_range_B": [str(fB["coll_year"].min()), str(fB["coll_year"].max())]}

# M2 extended window (−60, 2920] on original non-censoring frame
fE = make_frame(hi=2920)
m2 = {"n_2920": int(len(fE))}
m2["equals_nowindow"] = bool(len(fE) == len(fA))

# ---------------- Part 3: M4 AFT + RMST ----------------
print("=" * 60); print("PART 4: M4 AFT + RMST")
from lifelines import WeibullAFTFitter, LogNormalAFTFitter
aft_out = {}
for frame_name, frame in [("nowindow", fA), ("complete_followup", fB)]:
    for name, Fitter in [("weibull", WeibullAFTFitter), ("lognormal", LogNormalAFTFitter)]:
        try:
            f = Fitter()
            f.fit(frame[["dur", "event", "pool_hub", "pool_lmic", "coll_year_c"]],
                  duration_col="dur", event_col="event")
            smry = f.summary
            keep = [k for k in smry.index if "pool_" in k]
            aft_out[f"{frame_name}_{name}"] = {k: {"coef": float(smry.loc[k, "coef"]),
                                         "exp": float(np.exp(smry.loc[k, "coef"])),
                                         "lo": float(np.exp(smry.loc[k, "coef lower 95%"])),
                                         "hi": float(np.exp(smry.loc[k, "coef upper 95%"])),
                                         "p": float(smry.loc[k, "p"])} for k in keep}
        except Exception as e:
            aft_out[f"{frame_name}_{name}"] = {"error": str(e)[:200]}

km_out = {}
for frame_name, frame in [("nowindow", fA), ("complete_followup", fB)]:
    for tau in (365, 730, 1095):
        rmsts = {}
        for pool in ["indep_HIC_dom", "indep_LMIC", "hub_routed"]:
            sub = frame[frame["pool"] == pool]
            kmf = KaplanMeierFitter(); kmf.fit(sub["dur"], sub["event"])
            rmsts[pool] = float(restricted_mean_survival_time(kmf, t=tau))
        km_out[f"{frame_name}_tau_{tau}"] = rmsts
print("AFT:", json.dumps(aft_out, indent=1))
print("RMST:", json.dumps(km_out, indent=1))

# bootstrap moved to bootstrap_fast.py

# ---------------- Part 5: strata 2/3 window + censoring counts ----------------
print("=" * 60); print("PART 5: strata 2/3")
strata23 = {}
for name, fn in [("multipathogen", "multipathogen_biosample.csv"),
                 ("metagenome", "metagenome_sample.csv")]:
    d = build_records(os.path.join(DATA, fn), ASC_MULTI)
    st = {"n": int(len(d)),
          "valid_pair": int(d["lag"].notna().sum()),
          "invalid_pair": int(d["lag"].isna().sum()),
          "lag_gt_2920": int((d["lag"] > 2920).sum()),
          "lag_le_minus60": int((d["lag"] <= -60).sum()),
          "lag_neg_retained": int(((d["lag"] > -60) & (d["lag"] < 0)).sum()),
          "max_lag": int(d["lag"].max()) if d["lag"].notna().any() else None,
          "pool_all": d["pool"].value_counts().to_dict()}
    inv = d[d["lag"].isna()]
    st["invalid_by_pool"] = inv["pool"].value_counts().to_dict()
    strata23[name] = st
print(json.dumps(strata23, indent=1))

# metagenome country counts
meta = pd.read_csv(os.path.join(DATA, "metagenome_sample.csv"))
meta["c"] = meta["geo_location"].map(lambda g: sc(g) if isinstance(g, str) else "(missing)")
cc = meta["c"].value_counts()
us_share = cc.get("USA", 0) / 3000; cn_share = cc.get("CHN", 0) / 3000
meta_countries = int(meta[meta["c"] != "(missing)"]["c"].nunique())
m8 = {"us_share": us_share, "us_extrapolated": us_share * 373151,
      "us_ci_half": 1.96 * np.sqrt(us_share * (1 - us_share) / 3000) * 373151,
      "cn_share": cn_share, "cn_extrapolated": cn_share * 373151,
      "cn_ci_half": 1.96 * np.sqrt(cn_share * (1 - cn_share) / 3000) * 373151,
      "n_countries_meta": meta_countries}
multi = pd.read_csv(os.path.join(DATA, "multipathogen_biosample.csv"))
multi["c"] = multi["geo_location"].map(lambda g: sc(g) if isinstance(g, str) else "(missing)")
m8["n_countries_multi"] = int(multi[multi["c"] != "(missing)"]["c"].nunique())
print("M8 metagenome:", json.dumps(m8, indent=1))

# ---------------- Part 6: M3 presence model: Firth + drop-Peru ----------------
print("=" * 60); print("PART 6: M3 Firth logistic + drop-Peru")
import statsmodels.formula.api as smf

def firth_logit(X, y, max_iter=200, tol=1e-8):
    """Firth penalised logistic regression (Jeffreys prior), returns beta."""
    X = np.asarray(X, float); y = np.asarray(y, float)
    n, p = X.shape
    X = np.column_stack([np.ones(n), X]) if not np.allclose(X[:, 0], 1) else X
    n, p = X.shape
    beta = np.zeros(p)
    for it in range(max_iter):
        eta = X @ beta
        pr = 1 / (1 + np.exp(-eta))
        W = pr * (1 - pr)
        XW = X * W[:, None]
        I = X.T @ XW
        I_inv = np.linalg.pinv(I)
        H = np.sum((X @ I_inv) * XW, axis=1)  # hat diagonal
        U = X.T @ (y - pr + H * (0.5 - pr))
        step = I_inv @ U
        beta_new = beta + step
        if np.max(np.abs(step)) < tol:
            beta = beta_new; break
        beta = beta_new
    return beta, I_inv, it

def firth_se(X, y):
    Xd = np.asarray(X, float)
    b, Iinv, _ = firth_logit(Xd, y)
    return b, np.sqrt(np.diag(Iinv))

who = pd.read_csv(os.path.join(DATA, "who_covid_global.csv"))
latest = who.sort_values("Date_reported").groupby("Country_code").tail(1)
wb = json.load(open(os.path.join(BASE, "worldbank_cache.json"), encoding="utf-8"))
lat = latest[["Country", "Country_code", "Cumulative_deaths"]].copy()
lat["iso3"] = lat["Country_code"].map(lambda x: wb.get(x, {}).get("iso3"))
lat["income"] = lat["Country_code"].map(lambda x: wb.get(x, {}).get("income"))
lat["pop"] = lat["Country_code"].map(lambda x: wb.get(x, {}).get("pop2021"))
REC = {"GB": 49646, "NL": 17405, "DK": 8108, "AT": 3497, "CH": 1765, "SI": 616,
       "US": 563, "IN": 311, "SE": 220, "NO": 214, "FR": 147, "ET": 79,
       "LI": 31, "IT": 24, "DE": 22, "UY": 10, "CL": 9, "ES": 3, "CA": 2, "PE": 1}
lat["records"] = lat["Country_code"].map(REC).fillna(0).astype(float)
lat = lat.dropna(subset=["iso3", "income", "pop"]).copy()
lat["log_deaths"] = np.log(lat["Cumulative_deaths"].astype(float) + 1)
lat["log_pop"] = np.log(lat["pop"].astype(float))
lat["income"] = lat["income"].replace({"High income": "HIC", "Upper middle income": "UMIC",
                                       "Lower middle income": "LMIC", "Low income": "LIC"})
lat["present"] = (lat["records"] > 0).astype(int)

def design(d):
    X = pd.DataFrame({
        "log_deaths": d["log_deaths"], "log_pop": d["log_pop"],
        "UMIC": (d["income"] == "UMIC").astype(float),
        "LMIC": (d["income"] == "LMIC").astype(float),
        "LIC": (d["income"] == "LIC").astype(float)})
    return X, d["present"].values

names = ["Intercept", "log_deaths", "log_pop", "UMIC", "LMIC", "LIC"]
def run_firth(d):
    X, y = design(d)
    b, se = firth_se(X, y)
    out = {}
    for i, nm in enumerate(names):
        out[nm] = {"OR": float(np.exp(b[i])),
                   "lo": float(np.exp(b[i] - 1.96 * se[i])),
                   "hi": float(np.exp(b[i] + 1.96 * se[i]))}
    return out

m3 = {"n": int(len(lat)), "n_events": int(lat["present"].sum())}
m3["firth_full"] = run_firth(lat)
# drop-Peru
m3["firth_drop_peru"] = run_firth(lat[lat["Country"] != "Peru"])
# ordinary logit drop-Peru
lg = smf.logit("present ~ log_deaths + log_pop + C(income, Treatment('HIC'))",
               data=lat[lat["Country"] != "Peru"]).fit(disp=0)
smry = lg.summary2().tables[1]
m3["logit_drop_peru"] = {k: {"OR": float(np.exp(smry.loc[k, "Coef."])),
                             "lo": float(np.exp(smry.loc[k, "Coef."] - 1.96 * smry.loc[k, "Std.Err."])),
                             "hi": float(np.exp(smry.loc[k, "Coef."] + 1.96 * smry.loc[k, "Std.Err."]))}
                         for k in smry.index}
# how many UMIC events / LMIC events
m3["events_by_income"] = lat[lat["present"] == 1].groupby("income")["present"].count().to_dict()
print(json.dumps(m3, indent=1, default=str))

# ---------------- save ----------------
result = {"M1_cox_refits": m1, "M2_window": m2, "M3_presence": m3,
          "M4_aft": aft_out, "M4_rmst": km_out, "M5_accounting": acc,
          "M8_metagenome": m8, "strata23": strata23,
          "cluster_bootstrap": json.load(open(os.path.join(BASE, "cluster_bootstrap.json"))) if os.path.exists(os.path.join(BASE, "cluster_bootstrap.json")) else None,
          "versions": {"lifelines": __import__("lifelines").__version__,
                       "statsmodels": __import__("statsmodels").__version__,
                       "scipy": __import__("scipy").__version__,
                       "numpy": np.__version__, "pandas": pd.__version__}}
with open(os.path.join(BASE, "reviewer_response_stats.json"), "w", encoding="utf-8") as f:
    json.dump(result, f, ensure_ascii=False, indent=2, default=str)
print("saved reviewer_response_stats.json")

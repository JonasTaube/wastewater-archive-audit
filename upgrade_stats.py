# -*- coding: utf-8 -*-
"""统计升级 (v2.0): Cox比例风险模型 + 零膨胀负二项回归
Part A: 提交滞后(lag) ~ 通道(pool) 的Cox模型, 控制采样年份
Part B: WHO 239实体 记录数 ~ 死亡负担 + 人口 + 收入组 的ZINB模型
输出: stats_upgrade.json
口径: 与 PH_submission_final.md / round2_analysis.json 完全一致
  - hub: owner in (EBI, European Bioinformatics Institute)
  - indep_HIC_dom: 独立账户注册国==采样国 且采样国为HIC
  - indep_LMIC: 独立账户注册国==采样国 且采样国为IND/ETH
  - lag窗口: (-60, 1460] 天
"""
import csv, json, os, sys, io, time
from collections import defaultdict
from datetime import datetime

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.path.insert(0, r"D:\南开CDC\task\1、监测殖民主义")
import pandas as pd
import numpy as np
import requests
from map_analyze import parse_date

DST = r"D:\南开CDC\task\1、监测殖民主义"

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
HIC_COUNTRIES = {"USA", "GBR", "DNK", "SWE", "CHE", "NLD", "DEU", "FRA", "ITA",
                 "ESP", "CAN", "AUS", "JPN", "KOR", "SGP", "NZL", "IRL", "NOR",
                 "FIN", "AUT", "BEL", "ISR", "ARE", "SAU", "QAT", "HKG", "TWN",
                 "PRI", "URY", "CHL", "LIE"}

def sc(geo):
    from map_analyze import sample_country
    c = sample_country(geo)
    return "SWE" if c == "Sweden" else c

# ============ Part A: Cox 模型 ============
def build_lag_df():
    rows = list(csv.DictReader(open(DST + r"\master_biosample.csv", encoding="utf-8-sig")))
    recs = []
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
        if is_hub:
            pool = "hub_routed"
        elif oc == c:
            pool = "indep_LMIC" if c in ("IND", "ETH") else "indep_HIC_dom"
        else:
            continue  # 外国所有者独立提交(仅1条), 排除以对齐论文互斥池
        recs.append(dict(lag=max(lag, 0.5), pool=pool, coll_year=d2.year, country=c))
    return pd.DataFrame(recs)

def run_cox():
    from lifelines import CoxPHFitter
    from lifelines.statistics import proportional_hazard_test
    df = build_lag_df()
    out = {"n": int(len(df)),
           "pool_counts": df["pool"].value_counts().to_dict()}
    # 编码: ref = indep_HIC_dom
    df["pool_hub"] = (df["pool"] == "hub_routed").astype(int)
    df["pool_lmic"] = (df["pool"] == "indep_LMIC").astype(int)
    df["coll_year_c"] = df["coll_year"] - 2021
    # 主模型: lag ~ pool + coll_year
    cph = CoxPHFitter()
    fit_df = df[["lag", "pool_hub", "pool_lmic", "coll_year_c"]]
    fit_df = fit_df.assign(event=np.ones(len(fit_df)))
    cph.fit(fit_df, duration_col="lag", event_col="event")
    s = cph.summary
    out["cox_main"] = {
        "coef": {k: float(s.loc[k, "coef"]) for k in s.index},
        "HR": {k: float(s.loc[k, "exp(coef)"]) for k in s.index},
        "HR_ci_lower": {k: float(s.loc[k, "exp(coef) lower 95%"]) for k in s.index},
        "HR_ci_upper": {k: float(s.loc[k, "exp(coef) upper 95%"]) for k in s.index},
        "p": {k: float(s.loc[k, "p"]) for k in s.index},
    }
    # 分层敏感性: 按2022前后分层 (子集需含全部池, 否则记缺失)
    sens = {}
    for label, sub in [("collected<=2021", df[df["coll_year"] <= 2021]),
                       ("collected>=2022", df[df["coll_year"] >= 2022])]:
        pools = sub["pool"].value_counts().to_dict()
        if "indep_HIC_dom" not in pools or "hub_routed" not in pools:
            sens[label] = {"n": int(len(sub)), "skipped": "pool missing", "pool_counts": pools}
            continue
        if "indep_LMIC" not in pools:
            f2 = sub[["lag", "pool_hub"]].assign(event=1)
        else:
            f2 = sub[["lag", "pool_hub", "pool_lmic"]].assign(event=1)
        try:
            c2 = CoxPHFitter()
            c2.fit(f2, duration_col="lag", event_col="event")
            s2 = c2.summary
            sens[label] = {"n": int(len(sub)), "pool_counts": pools,
                           "HR_hub": float(s2.loc["pool_hub", "exp(coef)"]),
                           "HR_hub_ci": [float(s2.loc["pool_hub", "exp(coef) lower 95%"]),
                                         float(s2.loc["pool_hub", "exp(coef) upper 95%"])],
                           "p_hub": float(s2.loc["pool_hub", "p"])}
            if "indep_LMIC" in pools:
                sens[label]["HR_lmic"] = float(s2.loc["pool_lmic", "exp(coef)"])
        except Exception as e:
            sens[label] = {"n": int(len(sub)), "pool_counts": pools, "error": str(e)[:150]}
    out["cox_stratified"] = sens
    return out

# ============ Part B: ZINB 模型 ============
def fetch_worldbank():
    """收入分组 + 2021人口, 缓存。"""
    cache = os.path.join(DST, "worldbank_cache.json")
    if os.path.exists(cache):
        return json.load(open(cache, encoding="utf-8"))
    s = requests.Session(); s.trust_env = False
    # 收入分组
    r = s.get("https://api.worldbank.org/v2/country?format=json&per_page=400", timeout=60)
    wb = r.json()[1]
    info = {c["iso2Code"]: {"iso3": c["id"], "name": c["name"],
                            "income": c["region"]["value"] and c["incomeLevel"]["value"]}
            for c in wb if c["region"]["value"] != "Aggregates"}
    time.sleep(1)
    # 人口 2021
    r = s.get("https://api.worldbank.org/v2/country/all/indicator/SP.POP.TOTL"
              "?format=json&per_page=20000&date=2021", timeout=120)
    for item in r.json()[1]:
        iso3 = item["countryiso3code"]
        val = item["value"]
        if val is not None:
            for d in info.values():
                if d["iso3"] == iso3:
                    d["pop2021"] = int(val)
    json.dump(info, open(cache, "w", encoding="utf-8"), ensure_ascii=False)
    return info

def run_zinb():
    import statsmodels.formula.api as smf
    import statsmodels.api as smapi
    who = pd.read_csv(DST + r"\who_covid_global.csv", encoding="utf-8-sig")
    latest = who.sort_values("Date_reported").groupby("Country_code").tail(1)
    lat = latest[["Country", "Country_code", "Cumulative_deaths"]].copy()
    wb = fetch_worldbank()
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
    lat["income"] = lat["income"].replace({"High income": "HIC",
                                           "Upper middle income": "UMIC",
                                           "Lower middle income": "LMIC",
                                           "Low income": "LIC"})
    n_before = 239
    out = {"n_entities_input": n_before, "n_entities_modeled": int(len(lat)),
           "n_with_records": int((lat["records"] > 0).sum())}
    # ---- ZINB ----
    try:
        zinb = smapi.ZeroInflatedNegativeBinomialP.from_formula(
            "records ~ log_deaths + log_pop + C(income)",
            exog_infl=lat[["log_deaths", "log_pop"]],
            data=lat, infl_formula=None).fit(maxiter=200, disp=False)
        # IRR (incidence rate ratios)
        params = zinb.params
        conf = zinb.conf_int()
        names = list(params.index)
        irr = {n: float(np.exp(params[n])) for n in names if not n.startswith("inflate")}
        irrci = {n: [float(np.exp(conf.loc[n, 0])), float(np.exp(conf.loc[n, 1]))]
                 for n in names if not n.startswith("inflate")}
        out["zinb"] = {"converged": bool(zinb.mle_retvals.get("converged", True)),
                       "IRR": irr, "IRR_CI": irrci,
                       "pvalues": {n: float(zinb.pvalues[n]) for n in names},
                       "llf": float(zinb.llf), "aic": float(zinb.aic)}
    except Exception as e:
        out["zinb_error"] = str(e)
        # 回退: 负二项 GLM
        try:
            nb = smf.glm("records ~ log_deaths + log_pop + C(income)",
                         data=lat, family=smapi.families.NegativeBinomial(alpha=1.0)).fit()
            out["nb_glm_fallback"] = {
                "IRR": {n: float(np.exp(nb.params[n])) for n in nb.params.index},
                "p": {n: float(nb.pvalues[n]) for n in nb.params.index}}
        except Exception as e2:
            out["nb_error"] = str(e2)
    # ---- 两部模型: presence logit (更稳, 直接可解释) ----
    lat["present"] = (lat["records"] > 0).astype(int)
    logit = smf.logit("present ~ log_deaths + log_pop + C(income)", data=lat).fit(disp=0)
    out["presence_logit"] = {
        "n": int(len(lat)), "n_events": int(lat["present"].sum()),
        "OR": {n: float(np.exp(logit.params[n])) for n in logit.params.index},
        "OR_CI": {n: [float(np.exp(logit.conf_int().loc[n, 0])),
                      float(np.exp(logit.conf_int().loc[n, 1]))] for n in logit.params.index},
        "p": {n: float(logit.pvalues[n]) for n in logit.params.index}}
    return out

if __name__ == "__main__":
    result = {"ascertained": "2026-11-02"}
    print("=== Part A: Cox ===", flush=True)
    result["cox"] = run_cox()
    print(json.dumps(result["cox"]["cox_main"], indent=1), flush=True)
    print("=== Part B: ZINB ===", flush=True)
    result["zinb"] = run_zinb()
    with open(os.path.join(DST, "stats_upgrade.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
    print("saved stats_upgrade.json", flush=True)

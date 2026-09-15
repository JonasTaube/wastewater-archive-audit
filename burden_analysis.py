# -*- coding: utf-8 -*-
"""
Burden-surveillance decoupling analysis for Science Policy Article (v2, WHO source).
Joins 82,757 BioSample wastewater records with WHO COVID-19 cumulative deaths.
Outputs: burden_decoupling.json + burden_country.csv + console summary.
"""
import json, sys
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import requests

sys.stdout.reconfigure(encoding="utf-8")
BASE = r"D:\南开CDC\task\1、监测殖民主义"

# ---- 1. Load master records, parse sampling country ----
df = pd.read_csv(BASE + r"\master_biosample.csv", encoding="utf-8-sig", low_memory=False)
print(f"records total: {len(df)}")

def parse_country(g):
    if pd.isna(g):
        return None
    g = str(g).strip()
    if g.lower().startswith("not provided") or g == "":
        return None
    return g.split(":")[0].strip()

df["sample_country"] = df["geo_location"].map(parse_country)

HIC = {"GBR","NLD","DNK","AUT","CHE","SVN","USA","SWE","NOR","FRA","LIE","ITA","DEU","URY","CHL","ESP","CAN"}
LMIC_WORLD = {"IND","ETH"}
NAME2ISO = {
    "United Kingdom":"GBR","Netherlands":"NLD","Denmark":"DNK","Austria":"AUT","Switzerland":"CHE",
    "Slovenia":"SVN","USA":"USA","Sweden":"SWE","Norway":"NOR","France":"FRA","Liechtenstein":"LIE",
    "Italy":"ITA","Germany":"DEU","Uruguay":"URY","Chile":"CHL","Spain":"ESP","Canada":"CAN",
    "India":"IND","Ethiopia":"ETH","Peru":"PER",
}
df["iso"] = df["sample_country"].map(lambda c: NAME2ISO.get(c))
missing_geo = int(df["iso"].isna().sum())
print(f"records without geolocation: {missing_geo}")

rec = df.dropna(subset=["iso"]).groupby("iso").size().rename("records").reset_index()
rec["income"] = rec["iso"].map(lambda c: "HIC" if c in HIC else ("LMIC" if c in LMIC_WORLD else "UMIC"))

# ---- 2. WHO COVID-19 cumulative deaths (official global CSV) ----
who = pd.read_csv(BASE + r"\who_covid_global.csv", encoding="utf-8-sig")
who.columns = [c.strip() for c in who.columns]
print("WHO columns:", list(who.columns))
who["Date_reported"] = pd.to_datetime(who["Date_reported"])
latest = who.sort_values("Date_reported").groupby("Country_code").tail(1)
burden = latest[["Country_code","Country","Cumulative_deaths","Date_reported"]].copy()
burden = burden.rename(columns={"Cumulative_deaths":"total_deaths"})
print(f"WHO countries: {len(burden)}; as-of date: {burden['Date_reported'].max().date()}")
print(f"WHO world cumulative deaths: {burden['total_deaths'].sum():,.0f}")

# ---- 3. ISO2 -> ISO3 via World Bank API (also validates income groups) ----
s = requests.Session(); s.trust_env = False
r = s.get("https://api.worldbank.org/v2/country", params={"format":"json","per_page":400}, timeout=60)
wb = pd.json_normalize(r.json()[1])
wb = wb[wb["region.value"] != "Aggregates"].copy()
wb["iso2"] = wb["iso2Code"]
wb["iso3"] = wb["id"]
wb["wb_income"] = wb["incomeLevel.value"]
map_df = wb[["iso2","iso3","wb_income"]]
burden = burden.merge(map_df, left_on="Country_code", right_on="iso2", how="left")
print(f"countries with WB iso3 mapping: {burden['iso3'].notna().sum()}/{len(burden)}")

# ---- 4. Join & decoupling stats ----
m = rec.merge(burden[["iso3","Country","total_deaths","wb_income"]], left_on="iso", right_on="iso3", how="left")
world_deaths = burden["total_deaths"].sum()
world_records = int(rec["records"].sum())
print(f"joined: {len(m)} countries; audited records: {world_records:,}")

m["deaths_share_pct"] = 100*m["total_deaths"]/world_deaths
m["records_share_pct"] = 100*m["records"]/world_records
m["records_per_100k_deaths"] = 1e5*m["records"]/m["total_deaths"]
m["mismatch_ratio"] = m["records_share_pct"]/m["deaths_share_pct"].replace(0, np.nan)

r_all, p_all = spearmanr(m["total_deaths"], m["records"])
mh = m[m["income"]=="HIC"]
r_hic, p_hic = spearmanr(mh["total_deaths"], mh["records"])
print("\n--- Spearman deaths vs records ---")
print(f"ALL joined countries (n={len(m)}): rho={r_all:.3f}, P={p_all:.3g}")
print(f"HIC only (n={len(mh)}):            rho={r_hic:.3f}, P={p_hic:.3g}")

print("\n--- Burden top-10 countries vs their wastewater records ---")
top10 = m.nlargest(10,"total_deaths")
for _,r_ in top10.iterrows():
    print(f"  {r_['Country']:<28} deaths={r_['total_deaths']:>12,.0f} ({r_['deaths_share_pct']:5.2f}% world)  records={int(r_['records']):>6,} ({r_['records_share_pct']:6.3f}% archive)  mismatch x{r_['mismatch_ratio']:.3f}")

print("\n--- Key rows ---")
for iso in ["USA","GBR","BRA","IND","ETH","PER","DEU","FRA"]:
    sub = m[m["iso"]==iso]
    if len(sub):
        r_ = sub.iloc[0]
        print(f"  {iso}: {r_['Country']:<25} deaths={r_['total_deaths']:>12,.0f} ({r_['deaths_share_pct']:5.2f}%)  records={int(r_['records']):>6,} ({r_['records_share_pct']:6.3f}%)  rec/100k deaths={r_['records_per_100k_deaths']:.3f}  mismatch x{r_['mismatch_ratio']:.3f}")

# ---- 5. Save ----
m_out = m[["iso","Country","income","wb_income","records","total_deaths",
           "deaths_share_pct","records_share_pct","records_per_100k_deaths","mismatch_ratio"]]
m_out = m_out.sort_values("total_deaths", ascending=False)
m_out.to_csv(BASE + r"\burden_country.csv", index=False)

stats = {
    "source": "WHO COVID-19 Dashboard, global data CSV (data.who.int/dashboards/covid19), downloaded 15 Sep 2026",
    "who_asof_date": str(burden["Date_reported"].max().date()),
    "who_world_cumulative_deaths": float(world_deaths),
    "n_records_audited": world_records,
    "n_records_missing_geo": missing_geo,
    "n_countries_joined": int(len(m)),
    "spearman_all": {"rho": float(r_all), "p": float(p_all), "n": int(len(m))},
    "spearman_hic": {"rho": float(r_hic), "p": float(p_hic), "n": int(len(mh))},
}
for iso in ["USA","GBR","BRA","IND","ETH","PER"]:
    sub = m[m["iso"]==iso]
    if len(sub):
        r_ = sub.iloc[0]
        stats[iso.lower()] = {
            "country": r_["Country"], "deaths": float(r_["total_deaths"]),
            "deaths_share_pct": float(r_["deaths_share_pct"]),
            "records": int(r_["records"]), "records_share_pct": float(r_["records_share_pct"]),
            "records_per_100k_deaths": float(r_["records_per_100k_deaths"]),
            "mismatch_ratio": None if pd.isna(r_["mismatch_ratio"]) else float(r_["mismatch_ratio"]),
        }
with open(BASE + r"\burden_decoupling.json","w",encoding="utf-8") as f:
    json.dump(stats,f,indent=2,ensure_ascii=False)
print("\nSaved burden_country.csv + burden_decoupling.json")

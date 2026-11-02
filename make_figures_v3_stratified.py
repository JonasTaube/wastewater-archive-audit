# -*- coding: utf-8 -*-
"""Figure rendering for stratified-audit manuscript v2 (2026-11-02)
Fig 2c  : Kaplan-Meier time-to-submission curves, three pools + Cox HR
Fig 3b  : two-part presence model - predicted presence probability by income group
          + contributors-only scatter (deaths vs records)
Fig 4   : stratified custody architecture (stacked bars / country counts & sizes / lag)
Spec: Times New Roman, pdf.fonttype 42, white background, no in-figure titles,
      journal double-column width 7.2 in, panel letters a/b/c.
Output: D:\\南开CDC\\2026课题\\青年\\figures_v2\\ (PDF + 300 dpi PNG)
"""
import csv, json, os, sys, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.path.insert(0, r"D:\南开CDC\task\1、监测殖民主义")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.stats as st
from matplotlib import font_manager

for f in ("times.ttf", "timesbd.ttf", "timesi.ttf", "timesbi.ttf"):
    p = os.path.join(r"C:\Windows\Fonts", f)
    if os.path.exists(p):
        font_manager.fontManager.addfont(p)

# ---------- global style ----------
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
})

DST = r"D:\南开CDC\task\1、监测殖民主义"
OUT = r"D:\南开CDC\2026课题\青年\figures_v2"
os.makedirs(OUT, exist_ok=True)

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

C_HUB = "#B2182B"      # gateway-routed (red)
C_DOM = "#2166AC"      # domestic (blue)
C_LMIC = "#E08214"     # LMIC independent (amber)
C_META = "#4393C3"
C_PATH = "#92C5DE"
C_GREY = "#878787"


# ================= Fig 2c : Kaplan-Meier =================
def build_lag_df():
    from map_analyze import parse_date, sample_country
    rows = list(csv.DictReader(open(DST + r"\master_biosample.csv", encoding="utf-8-sig")))
    recs = []
    for r in rows:
        d1, d2 = parse_date(r["submission_date"]), parse_date(r["collection_date"])
        if not d1 or not d2:
            continue
        lag = (d1 - d2).days
        if not (-60 < lag <= 1460):
            continue
        c = sample_country(r["geo_location"])
        c = "SWE" if c == "Sweden" else c
        is_hub = r["owner"] in HUB_OWNERS
        oc = OWNER_CTRY.get(r["owner"], "?")
        if is_hub:
            pool = "Hub-routed"
        elif oc == c:
            pool = "Indep. LMIC" if c in ("IND", "ETH") else "Indep. HIC domestic"
        else:
            continue
        recs.append(dict(lag=lag, pool=pool))
    return pd.DataFrame(recs)


def fig2c():
    from lifelines import KaplanMeierFitter
    df = build_lag_df()
    fig, ax = plt.subplots(figsize=(3.5, 2.7))  # single column
    order = [("Indep. HIC domestic", C_DOM, "n = 543"),
             ("Indep. LMIC", C_LMIC, "n = 389"),
             ("Hub-routed", C_HUB, "n = 78,432")]
    for label, color, nn in order:
        sub = df[df["pool"] == label]
        kmf = KaplanMeierFitter()
        kmf.fit(sub["lag"], event_observed=np.ones(len(sub)),
                label=f"{label} ({nn})")
        kmf.plot_survival_function(ax=ax, ci_show=False, color=color,
                                   linewidth=1.2, show_censors=False)
    ax.set_ylim(0, 1.0)
    ax.set_xlim(0, 1460)
    ax.set_xlabel("Days from collection to public submission")
    ax.set_ylabel("Proportion not yet submitted")
    ax.invert_yaxis()  # curves descend as records arrive
    meds = df.groupby("pool")["lag"].median()
    for label, color, _ in order:
        m = meds[label]
        ax.axvline(m, color=color, lw=0.6, ls=":", alpha=0.7)
    txt = ("Cox model (adjusted for collection year):\n"
           "Hub-routed  HR 0.093 (95% CI 0.085\u20130.101)\n"
           "Indep. LMIC  HR 0.245 (95% CI 0.215\u20130.279)")
    ax.text(0.97, 0.05, txt, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=6.2,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", lw=0.5))
    ax.legend(frameon=False, loc="upper right")
    fig.savefig(os.path.join(OUT, "fig2c.pdf"))
    fig.savefig(os.path.join(OUT, "fig2c.png"), dpi=300)
    plt.close(fig)
    print("saved: fig2c (n = %d)" % len(df), flush=True)


# ================= Fig 3b : presence model =================
def fig3b():
    import statsmodels.formula.api as smf
    who = pd.read_csv(DST + r"\who_covid_global.csv", encoding="utf-8-sig")
    latest = who.sort_values("Date_reported").groupby("Country_code").tail(1)
    lat = latest[["Country", "Country_code", "Cumulative_deaths"]].copy()
    wb = json.load(open(DST + r"\worldbank_cache.json", encoding="utf-8"))
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
    lat["present"] = (lat["records"] > 0).astype(int)
    logit = smf.logit("present ~ log_deaths + log_pop + C(income, "
                      "Treatment(reference='HIC'))", data=lat).fit(disp=0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 2.7))
    # --- panel b (left): predicted presence probability by income group
    grid = pd.DataFrame({
        "log_deaths": [lat["log_deaths"].median()] * 4,
        "log_pop": [lat["log_pop"].median()] * 4,
        "income": ["HIC", "UMIC", "LMIC", "LIC"]})
    pred = logit.predict(grid)
    ci = np.asarray(logit.get_prediction(grid).conf_int())
    cols = {"HIC": C_DOM, "UMIC": C_LMIC, "LMIC": "#F4A582", "LIC": C_GREY}
    for i, g in enumerate(grid["income"]):
        lo, hi = ci[i, 0], ci[i, 1]
        ax1.vlines(i, lo, hi, color=cols[g], lw=1.3)
        ax1.plot(i, pred.iloc[i], "o", color=cols[g], ms=5.5,
                 mec="black", mew=0.4)
    ax1.set_xticks(np.arange(4))
    ax1.set_xticklabels(grid["income"])
    ax1.set_ylabel("Predicted probability of\narchive presence")
    ax1.set_xlabel("World Bank income group")
    ax1.set_ylim(-0.03, 1.08)
    ax1.text(0.03, 0.97, "at median mortality and population\n(n = 213 WHO entities)",
             transform=ax1.transAxes, fontsize=6.2, style="italic", va="top")
    ax1.text(0.03, 0.72,
             "log deaths  OR 2.50 (1.13\u20135.52)\nUMIC vs HIC  OR 0.032 (0.004\u20130.29)",
             transform=ax1.transAxes, fontsize=6.2, va="top",
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", lw=0.5))
    # --- panel b (right): contributors scatter
    sub = lat[lat["records"] > 0].copy()
    rho, p = st.spearmanr(sub["Cumulative_deaths"], sub["records"])
    ax2.scatter(sub["Cumulative_deaths"], sub["records"], s=15, color=C_DOM,
                alpha=0.85, edgecolor="black", linewidth=0.3)
    annot = {"GB": "United Kingdom", "IN": "India", "US": "United States",
             "PE": "Peru", "ET": "Ethiopia", "DK": "Denmark", "NL": "Netherlands"}
    for cc, name in annot.items():
        row = sub[sub["Country_code"] == cc]
        if len(row):
            dx, dy = (4, 3)
            if cc == "IN":
                dx, dy = (-46, -7)
            elif cc == "NL":
                dx, dy = (6, -10)
            ax2.annotate(name, (row["Cumulative_deaths"].iloc[0], row["records"].iloc[0]),
                         textcoords="offset points", xytext=(dx, dy), fontsize=6)
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_xlabel("Cumulative reported COVID-19 deaths (log scale)")
    ax2.set_ylabel("Wastewater sequencing\nrecords (log scale)")
    # Spearman box -> 右面板左下空白区(低死亡+低记录区无数据点)
    ax2.text(0.03, 0.04,
             f"contributors only (n = {len(sub)})\nSpearman \u03c1 = \u22120.08 "
             f"(95% CI \u22120.50 to 0.38)",
             transform=ax2.transAxes, fontsize=6.2, va="bottom", ha="left",
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", lw=0.5))
    ax1.text(-0.30, 1.06, "b", transform=ax1.transAxes, fontsize=10,
             fontweight="bold", va="top")
    fig.tight_layout(w_pad=2.5)
    fig.savefig(os.path.join(OUT, "fig3b.pdf"))
    fig.savefig(os.path.join(OUT, "fig3b.png"), dpi=300)
    plt.close(fig)
    print("saved: fig3b", flush=True)


# ================= Fig 4 : stratified architecture =================
def fig4():
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(7.2, 2.9))

    strata = ["SARS-CoV-2", "Other\npathogens", "Metagenome"]
    # --- panel a: stacked ownership architecture
    gateway = np.array([98.9, 1.2, 19.8])
    domestic = np.array([1.1, 76.5, 60.7])
    foreign = 100.0 - gateway - domestic
    x = np.arange(3)
    b1 = ax1.bar(x, gateway, width=0.55, color=C_HUB, label="Gateway-routed (EBI)")
    b2 = ax1.bar(x, domestic, width=0.55, bottom=gateway, color=C_DOM,
                 label="Domestic account")
    b3 = ax1.bar(x, foreign, width=0.55, bottom=gateway + domestic, color=C_GREY,
                 label="Foreign / unresolved")
    for i in range(3):
        for val, bottom in [(gateway[i], 0.0), (domestic[i], gateway[i])]:
            if val > 9:
                ax1.text(i, bottom + val / 2, f"{val:.1f}%", ha="center",
                         va="center", fontsize=6.5, color="white")
    # 1.1% 太窄, 引线标注
    ax1.annotate("1.1%", xy=(0.28, 99.0), xytext=(0.62, 92.0), fontsize=6.0,
                 arrowprops=dict(arrowstyle="-", lw=0.5, color="0.3"))
    ax1.set_xticks(x); ax1.set_xticklabels(strata, fontsize=7)
    ax1.set_ylabel("Share of records (%)")
    ax1.set_ylim(0, 100)
    ax1.set_yticks([0, 25, 50, 75, 100])
    ax1.text(0.02, 1.08, "a", transform=ax1.transAxes, fontsize=10,
             fontweight="bold", va="top")

    # --- panel b: countries represented + stratum sizes
    n_ctry = [20, 9, 57]
    sizes = [82757, 1426, 373151]
    cols3 = [C_DOM, C_PATH, C_META]
    ax2.bar(x - 0.19, n_ctry, width=0.38, color=cols3, edgecolor="black",
            linewidth=0.4)
    ax2.set_xticks(x); ax2.set_xticklabels(strata, fontsize=7)
    ax2.set_ylabel("Sampling countries (count)")
    ax2.set_ylim(0, 70)
    ax2b = ax2.twinx()
    ax2b.bar(x + 0.19, sizes, width=0.38, color=cols3, alpha=0.35,
             hatch="//", edgecolor="black", linewidth=0.4)
    ax2b.set_yscale("log")
    ax2b.set_yticks([1e2, 1e3, 1e4, 1e5, 1e6])
    ax2b.set_yticklabels(["10\u00b2", "10\u00b3", "10\u2074", "10\u2075", "10\u2076"],
                         fontsize=7)
    ax2b.set_ylim(100, 3e6)
    ax2b.set_ylabel("Stratum size (records)")
    ax2b.spines["right"].set_visible(True)
    ax2b.spines["top"].set_visible(False)
    for i, v in enumerate(n_ctry):
        ax2.text(i - 0.19, v + 1.5, str(v), ha="center", fontsize=6.5)
    for i, v in enumerate(sizes):
        ax2b.text(i + 0.19, v * 0.55, f"{v:,}", ha="center", fontsize=6.0,
                  rotation=90, color="0.25")
    ax2.text(0.02, 1.08, "b", transform=ax2.transAxes, fontsize=10,
             fontweight="bold", va="top")

    # --- panel c: median lag by stratum / pathway
    labels = ["SARS2 HIC dom.", "SARS2 LMIC", "SARS2 hub", "Other path.",
              "Metagenome"]
    lags = [99, 349, 890, 197, 267]
    cols5 = [C_DOM, C_LMIC, C_HUB, C_PATH, C_META]
    xpos = np.arange(5)
    ax3.bar(xpos, lags, width=0.62, color=cols5, edgecolor="black",
            linewidth=0.4)
    for i, v in enumerate(lags):
        ax3.text(i, v + 20, str(v), ha="center", fontsize=6.5)
    ax3.set_xticks(xpos)
    ax3.set_xticklabels(labels, fontsize=6.2, rotation=25, ha="right")
    ax3.set_ylabel("Median collection-to-\nsubmission lag (days)")
    ax3.set_ylim(0, 1010)
    ax3.text(0.02, 1.08, "c", transform=ax3.transAxes, fontsize=10,
             fontweight="bold", va="top")

    # 图例移到全图下方, 避免遮挡
    fig.tight_layout(w_pad=2.2)
    handles = [b1[0], b2[0], b3[0]]
    fig.legend(handles, ["Gateway-routed (EBI)", "Domestic account",
                         "Foreign / unresolved"],
               loc="lower center", ncol=3, fontsize=6.5, frameon=False,
               bbox_to_anchor=(0.5, -0.035))
    fig.savefig(os.path.join(OUT, "fig4.pdf"))
    fig.savefig(os.path.join(OUT, "fig4.png"), dpi=300)
    plt.close(fig)
    print("saved: fig4", flush=True)


if __name__ == "__main__":
    fig2c()
    fig3b()
    fig4()
    print("ALL FIGURES DONE ->", OUT, flush=True)

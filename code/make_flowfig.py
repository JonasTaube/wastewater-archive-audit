# -*- coding: utf-8 -*-
"""Supplementary Figure 3: STROBE-style record flow / accounting diagram."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

RED = "#b2182b"; BLUE = "#2166ac"; GREY = "#878787"; DARK = "#1a1a1a"
plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42})

fig, ax = plt.subplots(figsize=(7.6, 6.8))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

# main column x 32..99 ; side column x 1..25
MX, MW = 30, 69
SX, SW = 1, 25

def box(x, y, w, h, text, fc="#f2f2f2", ec=DARK, fs=8.0, bold=False, align="center"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.55,rounding_size=1.1",
                                fc=fc, ec=ec, lw=1.0))
    tx = x + w/2 if align == "center" else x + 1.6
    ax.text(tx, y + h/2, text, ha="center" if align == "center" else "left", va="center",
            fontsize=fs, color=DARK, fontweight="bold" if bold else "normal")

def arrow(x1, y1, x2, y2, style="-|>"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=9, color=DARK, lw=0.9))

# Level 1
box(MX, 89, MW, 8.5, "SARS-CoV-2 wastewater census (NCBI BioSample, 15 Sep 2026)\nn = 82,757",
    fc="#e8eef7", bold=True, fs=8.4)
box(SX, 78, SW, 7, "No geolocation\n84 (0.10%; all EBI-owned)", fc="#fafafa", ec=GREY, fs=7.2)
arrow(MX + 8, 89, SX + SW/2, 85.8)
# Level 2
box(MX, 75.5, MW, 10.5, "Ownership attribution\n"
    "EBI-owned (gateway) 81,813 (98.9%)\n"
    "independent accounts 944 (1.1%) = domestic 934\n"
    "(HIC 544 · LMIC 390) + foreign-owned 10 (Chile 9 · Peru 1)", fs=7.4)
arrow(65, 89, 65, 86.8)
# Level 3
box(SX, 63.5, SW, 7.5, "Foreign-attributed\n32,099 (38.79%)\n99.97% to the GBR", fc="#fafafa", ec=GREY, fs=7.2)
arrow(MX + 10, 76.5, SX + SW/2, 71.8)
# Level 4
box(MX, 59.5, MW, 12.5, "Date-pair validity\n"
    "valid collection\u2013submission pair: 82,616\n"
    "invalid pair: 141  (gateway 130 · foreign-owned Chile 9 ·\n"
    "indep. HIC 1 · indep. LMIC 1)\n"
    "foreign-owned independent with valid dates: 1 (Peru)\n"
    "\u2192 outside the three pools", fs=7.4)
arrow(65, 75.5, 65, 72.8)
# Level 5
box(SX, 48.5, SW, 7.5, "Lag > 1,460 days\n3,251 excluded\n(all gateway-routed)", fc="#fdf0f0", ec=RED, fs=7.2)
arrow(MX + 10, 60, SX + SW/2, 56.8)
box(MX, 45.5, MW, 10.5, "Lag window (\u221260, 1,460] days: 79,365\n"
    "pre-collection submissions (lag \u2264 \u221260 days): 0\n"
    "negative lags (\u221260 < lag < 0): 0\n"
    "(minimum observed lag: 8 days)", fs=7.4)
arrow(65, 59.5, 65, 56.8)
# Level 6
box(MX, 30.5, MW, 9.5, "Velocity-analysis set\nn = 79,364\n"
    "gateway-routed 78,432 · independent high-income domestic 543\n"
    "independent LMIC 389",
    fc="#e8f0e8", fs=7.6)
arrow(65, 45.5, 65, 40.8)
# Note
box(MX, 10, MW, 14, "Censoring note\n"
    "The census is a snapshot: records not yet submitted\n"
    "by the ascertainment date are unobservable, and the\n"
    "velocity analysis therefore conditions on submission\n"
    "having occurred. All 82,616 records with a valid\n"
    "collection date also carry a submission date, so no\n"
    "within-census administrative censoring arises.", fc="#ffffff", ec=GREY, fs=7.2)
arrow(65, 30.5, 65, 24.8, style="-")

fig.savefig("../submission_package/figures/Supplementary_Figure_3.pdf")
fig.savefig("../submission_package/figures/Supplementary_Figure_3.png", dpi=300)
print("saved Supplementary_Figure_3 v2")

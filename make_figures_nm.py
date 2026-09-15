# -*- coding: utf-8 -*-
"""Render Figs 1-4 for the Nature Medicine submission.

Nature Medicine figure specification:
  - sans-serif, 5-7 pt final size
  - max width 180 mm (= 7.09 in); height max 240 mm
  - 300 dpi PNG (for DOCX embedding) + PDF vector (for production)
  - pdf.fonttype 42 (TrueType embedded, editable text)
  - white background, English labels, NO in-figure titles (captions live in the manuscript)
  - panel letters bold lowercase

Fig 1  geographic concentration   | a world map (Europe inset) | b top-10 | c per-capita gradient
Fig 2  custodial concentration    | a chord: sampled -> submitting country | b domestic vs hub | c lag by pathway
Fig 3  decoupling + hiatus        | a quarterly release | b burden scatter | c records per 100k deaths
Fig 4  governance models          | four custodial arrangements and their failure modes
"""
import os, sys, csv, json, io, time
from collections import defaultdict
import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.path.insert(0, r"D:\南开CDC\task\1、监测殖民主义")

import matplotlib
import matplotlib.patches as mpatches
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Wedge, Polygon as MPoly, FancyBboxPatch
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from matplotlib import font_manager, gridspec
import matplotlib.patheffects as pe
import scipy.stats as st

DST = r"D:\南开CDC\task\1、监测殖民主义"
FIGD = os.path.join(DST, "figures")
os.makedirs(FIGD, exist_ok=True)

# ---- Arial (NM sans-serif requirement) ----
for f in ("arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"):
    p = os.path.join(r"C:\Windows\Fonts", f)
    if os.path.exists(p):
        font_manager.fontManager.addfont(p)
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 6,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 1.8, "ytick.major.size": 1.8,
    "axes.unicode_minus": False,
    "savefig.facecolor": "white",
})
print("font:", font_manager.findfont("Arial"))

from map_analyze import sample_country, parse_date

NAME = {"GBR": "United Kingdom", "NLD": "Netherlands", "DNK": "Denmark", "AUT": "Austria",
        "CHE": "Switzerland", "SVN": "Slovenia", "SWE": "Sweden", "NOR": "Norway",
        "USA": "United States", "LIE": "Liechtenstein", "FRA": "France", "ITA": "Italy",
        "DEU": "Germany", "CHL": "Chile", "CAN": "Canada", "PER": "Peru", "IND": "India",
        "ETH": "Ethiopia", "URY": "Uruguay", "ESP": "Spain"}
ISO2 = {"GBR": "UK", "NLD": "NL", "DNK": "DK", "AUT": "AT", "CHE": "CH", "SVN": "SI",
        "SWE": "SE", "NOR": "NO", "USA": "US", "LIE": "LI", "FRA": "FR", "ITA": "IT",
        "DEU": "DE", "CHL": "CL", "CAN": "CA", "PER": "PE", "IND": "IN", "ETH": "ET",
        "URY": "UY", "ESP": "ES"}

RED = "#B2182B"; BLUE = "#2166AC"; GREEN = "#1B7837"; GREY = "#878787"; LGREY = "#E9E9E9"
HALO = [pe.withStroke(linewidth=1.4, foreground="white")]

def sc(geo):
    c = sample_country(geo)
    return "SWE" if c == "Sweden" else c

rows = list(csv.DictReader(open(os.path.join(DST, "master_biosample.csv"), encoding="utf-8-sig")))
cstats = list(csv.DictReader(open(os.path.join(DST, "country_stats.csv"), encoding="utf-8-sig")))
for r in cstats:
    r["n_samples"] = int(r["n_samples"]); r["misaligned_n"] = int(r["misaligned_n"])
    if r["sample_country"] == "Sweden":
        r["sample_country"] = "SWE"
cstats.sort(key=lambda r: -r["n_samples"])
TOT = 82757
print(f"loaded {len(rows)} records, {len(cstats)} countries")

burden = list(csv.DictReader(open(os.path.join(DST, "burden_country.csv"), encoding="utf-8-sig")))
for b in burden:
    b["records"] = int(b["records"]); b["total_deaths"] = float(b["total_deaths"])
    b["records_per_100k_deaths"] = float(b["records_per_100k_deaths"])
    b["mismatch_ratio"] = float(b["mismatch_ratio"])
BUR = {b["iso"]: b for b in burden}
print("burden rows:", len(burden))

POP = {"GBR": 67.65, "NLD": 17.59, "DNK": 5.87, "AUT": 8.93, "CHE": 8.70, "SVN": 2.11,
       "USA": 331.9, "SWE": 10.42, "NOR": 5.41, "FRA": 67.8, "LIE": 0.039, "ITA": 59.24,
       "DEU": 83.2, "ESP": 47.4, "CAN": 38.2, "CHL": 19.6, "URY": 3.43,
       "IND": 1407.6, "ETH": 120.3, "PER": 34.0}   # UN WPP 2022, 2021 estimates, millions

def letter(ax, s, dx=-0.052, dy=1.035, fs=7):
    ax.text(dx, dy, s, transform=ax.transAxes, fontweight="bold", fontsize=fs,
            ha="left", va="bottom")

def save(fig, name):
    fig.savefig(os.path.join(FIGD, name + ".png"), dpi=300)
    fig.savefig(os.path.join(FIGD, name + ".pdf"))
    plt.close(fig)
    print("saved", name)

# ================================================================= Fig 1
def fig1():
    t0 = time.time()
    fig = plt.figure(figsize=(7.09, 6.4))   # 180 mm wide
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.30, 1.0], hspace=0.42, wspace=0.32,
                           left=0.105, right=0.975, top=0.975, bottom=0.085)
    axm = fig.add_subplot(gs[0, :]); axb = fig.add_subplot(gs[1, 0]); axc = fig.add_subplot(gs[1, 1])
    vals = {r["sample_country"]: r["n_samples"] for r in cstats}
    vmax = max(vals.values())
    cmap = LinearSegmentedColormap.from_list("wbe", ["#C6DBEF", "#6BAED6", "#2171B5", "#08306B"])
    norm = LogNorm(vmin=1, vmax=vmax)

    gjd = json.load(open(os.path.join(FIGD, "ne_50m_countries.geojson"), encoding="utf-8"))

    def polys_of(feat):
        geom = feat.get("geometry")
        if not geom:
            return []
        if geom["type"] == "Polygon":
            return [geom["coordinates"]]
        if geom["type"] == "MultiPolygon":
            return geom["coordinates"]
        return []

    def draw_world(ax, lw):
        for feat in gjd["features"]:
            a3 = (feat.get("properties") or {}).get("ADM0_A3", "")
            n = vals.get(a3, 0)
            fc = LGREY if n == 0 else cmap(norm(n))
            for rings in polys_of(feat):
                ext = np.asarray(rings[0], dtype=float)
                if len(ext) < 3:
                    continue
                ax.add_patch(MPoly(ext, closed=True, facecolor=fc,
                                   edgecolor="white", linewidth=lw, zorder=1))
        for feat in gjd["features"]:
            a3 = (feat.get("properties") or {}).get("ADM0_A3", "")
            if vals.get(a3, 0) > 0:
                continue
            for rings in polys_of(feat):
                for hole in rings[1:]:
                    h = np.asarray(hole, dtype=float)
                    if len(h) >= 3:
                        ax.add_patch(MPoly(h, closed=True, facecolor=LGREY,
                                           edgecolor="white", linewidth=lw * 0.6, zorder=2))

    draw_world(axm, 0.3)
    axm.set_xlim(-168, 186); axm.set_ylim(-58, 84)
    axm.set_aspect("equal"); axm.set_xticks([]); axm.set_yticks([])
    for s in axm.spines.values():
        s.set_visible(False)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cax = axm.inset_axes([0.115, 0.115, 0.20, 0.028])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_ticks([1, 10, 100, 1000, 10000])
    cb.ax.tick_params(labelsize=5, length=1.5, width=0.4)
    cb.set_label("Records per country (log scale)", fontsize=5.6)
    axi = axm.inset_axes([0.415, 0.36, 0.34, 0.60])
    draw_world(axi, 0.25)
    axi.set_xlim(-11, 33); axi.set_ylim(35.5, 66)
    axi.set_aspect("equal"); axi.set_xticks([]); axi.set_yticks([])
    for s in axi.spines.values():
        s.set_linewidth(0.5); s.set_edgecolor("#666666")
    eu = {"GBR": (-3.6, 53.2), "NLD": (5.6, 52.6), "DNK": (9.4, 56.4),
          "AUT": (14.6, 47.5), "CHE": (8.2, 46.6), "SVN": (14.9, 45.9),
          "SWE": (16.5, 63.0), "NOR": (7.0, 61.6), "FRA": (2.2, 46.8),
          "ITA": (12.8, 42.3), "DEU": (10.0, 51.0), "LIE": (17.8, 47.7)}
    for k, (x, y) in eu.items():
        axi.text(x, y, ISO2[k], fontsize=5, ha="center", va="center",
                 color="#1a1a1a", path_effects=HALO, zorder=5)
    axm.text(0.775, 0.90, "Europe (inset)", fontsize=5.2, color="#555555",
             transform=axm.transAxes)
    axm.text(0.012, 0.06, "n = 82,757 records\n20 countries\n84 lack geolocation",
             transform=axm.transAxes, fontsize=5.2, va="bottom",
             bbox=dict(boxstyle="round,pad=0.26", fc="white", ec="#BBBBBB", lw=0.4, alpha=0.93))
    letter(axm, "a", dx=-0.045, dy=1.005)

    top10 = cstats[:10]
    y = np.arange(len(top10))[::-1].astype(float)
    axb.barh(y, [r["n_samples"] for r in top10], color=BLUE, height=0.62)
    axb.set_yticks(y)
    axb.set_yticklabels([NAME[r["sample_country"]] for r in top10], fontsize=6)
    axb.set_xscale("log"); axb.set_xlim(90, 900000)
    for yi, r in zip(y, top10):
        share = 100.0 * r["n_samples"] / TOT
        axb.text(r["n_samples"] * 1.35, yi, f"{r['n_samples']:,}  ({share:.1f}%)",
                 va="center", fontsize=5.4, color="#222222")
    axb.set_xlabel("Records (log scale)", fontsize=6)
    axb.tick_params(axis="x", labelsize=5.4)
    axb.spines[["top", "right"]].set_visible(False)
    letter(axb, "b")

    pc = [(c, cnt / POP[c]) for c, cnt in vals.items()]
    keep = ["LIE", "DNK", "NLD", "GBR", "SVN", "AUT", "CHE", "NOR", "ETH", "IND", "PER", "USA"]
    pcs = sorted([(c, dict(pc)[c]) for c in keep], key=lambda t: t[1])
    xs = np.arange(len(pcs))
    barcol = [RED if c in ("GBR", "NLD", "DNK") else "#762A83" if c in ("IND", "ETH", "PER") else BLUE
              for c, _ in pcs]
    axc.bar(xs, [v for _, v in pcs], color=barcol, width=0.62)
    axc.set_yscale("log"); axc.set_ylim(0.015, 6e4)
    axc.set_xticks(xs)
    axc.set_xticklabels([ISO2[c] for c, _ in pcs], fontsize=5.4)
    axc.set_ylabel("Records per million residents (log)", fontsize=6)
    axc.tick_params(axis="y", labelsize=5.4)
    axc.spines[["top", "right"]].set_visible(False)
    d = dict(pcs)
    axc.text(0.03, 0.955, f"UK/India {d['GBR']/d['IND']:,.0f}\u00d7\nDK/India {d['DNK']/d['IND']:,.0f}\u00d7\nclinical genomes \u2248 100\u00d7",
             transform=axc.transAxes, fontsize=5.2, color="#67000D", va="top",
             bbox=dict(facecolor="white", edgecolor="#CCCCCC", lw=0.4, pad=1.4))
    axc.text(0.03, 0.38, "LMIC-world total: 391 records (0.47%)",
             transform=axc.transAxes, fontsize=5.2, color="#4D004B", va="bottom")
    letter(axc, "c")
    save(fig, "NM_Fig1")
    print("fig1 done in %.1fs" % (time.time() - t0))

# ================================================================= Fig 2
def fig2():
    t0 = time.time()
    flows = []
    for r in csv.DictReader(open(os.path.join(DST, "misalignment_matrix.csv"), encoding="utf-8-sig")):
        s = r["sampled_country"].strip(); s = "SWE" if s == "Sweden" else s
        flows.append((s, r["submitted_by_country"].strip(), int(r["n_samples"])))
    tot = defaultdict(float)
    for s, t, v in flows:
        tot[s] += v; tot[t] += v
    node_order = ["GBR", "CAN", "PER", "DEU", "ITA", "CHL", "FRA", "LIE", "USA",
                  "NOR", "SWE", "SVN", "CHE", "AUT", "DNK", "NLD"]
    total_flow = sum(v for _, _, v in flows)
    gap_units = total_flow * 0.004
    circ = sum(tot[n] for n in node_order) + gap_units * (len(node_order) - 1)
    dpu = 360.0 / circ
    arcs = {}
    a = 90.0
    for n in node_order:
        wdeg = tot[n] * dpu
        arcs[n] = (a, a - wdeg)
        a -= wdeg + gap_units * dpu

    node_ribbons = defaultdict(list)
    for i, (s, t, v) in enumerate(flows):
        node_ribbons[s].append(((0.5 * (arcs[t][0] + arcs[t][1]) - arcs[s][0]) % 360.0, i))
        node_ribbons[t].append(((0.5 * (arcs[s][0] + arcs[s][1]) - arcs[t][0]) % 360.0, i))
    spans = defaultdict(list)
    for n, lst in node_ribbons.items():
        lst.sort()
        off = 0.0
        for _, i in lst:
            w = flows[i][2] * dpu
            spans[n].append((i, off, off + w))
            off += w

    fig = plt.figure(figsize=(7.09, 7.9))
    ax = fig.add_axes([0.015, 0.545, 0.545, 0.44])
    ax.set_xlim(-1.58, 1.80); ax.set_ylim(-2.20, 1.50)
    ax.set_aspect("equal"); ax.axis("off")

    def polar(r, deg):
        rad = np.deg2rad(deg)
        return (r * np.cos(rad), r * np.sin(rad))

    R, WNODE, RC = 1.0, 0.038, 0.42
    rib_col = {"NLD": "#1F77B4", "DNK": "#FF7F0E", "AUT": "#2CA02C", "CHE": "#D62728",
               "SVN": "#9467BD", "SWE": "#8C564B", "NOR": "#E377C2", "USA": "#7F7F7F",
               "LIE": "#BCBD22", "FRA": "#17BECF", "ITA": "#AEC7E8", "DEU": "#FFBB78",
               "CHL": "#98DF8A", "PER": "#C5B0D5"}
    for i, (s, t, v) in enumerate(flows):
        sa0, sa1 = arcs[s]; da0, da1 = arcs[t]
        s_off0, s_off1 = [sp for sp in spans[s] if sp[0] == i][0][1:]
        d_off0, d_off1 = [sp for sp in spans[t] if sp[0] == i][0][1:]
        a_src_hi, a_src_lo = sa0 - s_off0, sa0 - s_off1
        a_dst_hi, a_dst_lo = da0 - d_off0, da0 - d_off1
        verts = [polar(R, a_src_hi)]
        verts += [polar(RC, a_src_hi), polar(RC, a_dst_hi), polar(R, a_dst_hi)]
        k = max(2, int(abs(a_dst_hi - a_dst_lo) / 1.2) + 1)
        for aa in np.linspace(a_dst_hi, a_dst_lo, k)[1:]:
            verts.append(polar(R, aa))
        verts += [polar(RC, a_dst_lo), polar(RC, a_src_lo), polar(R, a_src_lo)]
        k2 = max(2, int(abs(a_src_hi - a_src_lo) / 1.2) + 1)
        for aa in np.linspace(a_src_lo, a_src_hi, k2)[1:]:
            verts.append(polar(R, aa))
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        codes += [Path.LINETO] * (k - 1)
        codes += [Path.CURVE4, Path.CURVE4, Path.CURVE4]
        codes += [Path.LINETO] * (k2 - 1) + [Path.CLOSEPOLY]
        verts.append(verts[0])
        ax.add_patch(PathPatch(Path(verts, codes), facecolor=rib_col[s],
                               edgecolor=rib_col[s], lw=0.15, alpha=0.40, zorder=1))
    for n in node_order:
        a0, a1 = arcs[n]
        if n == "GBR":
            fc = "#67000D"
        elif n in ("CAN", "USA"):
            fc = "#525252"
        else:
            fc = rib_col[n]
        ax.add_patch(Wedge((0, 0), R, a1, a0, width=WNODE, facecolor=fc,
                           edgecolor="white", lw=0.2, zorder=3))
    placed = []

    def collides(x0, y0, w, h):
        for (px0, px1, py0, py1) in placed:
            if not (x0 > px1 + 0.02 or x0 + w < px0 - 0.02 or y0 > py1 + 0.015 or y0 + h < py0 - 0.015):
                return True
        return False

    small = sorted([n for n in node_order if abs(arcs[n][0] - arcs[n][1]) < 10.0],
                   key=lambda n: -((0.5 * (arcs[n][0] + arcs[n][1]) - 90.0) % 360.0))
    for n in small:
        a0, a1 = arcs[n]
        mid = 0.5 * (a0 + a1)
        txt = f"{ISO2[n]}  {int(tot[n]):,}"
        w = 0.14 + 0.0115 * len(txt)
        h = 0.085
        side = 1 if np.cos(np.deg2rad(mid)) >= 0 else -1
        rl = 1.045
        while rl < 2.6:
            x, yy = polar(rl, mid)
            x0 = x if side > 0 else x - w
            if not collides(x0, yy - h / 2, w, h):
                break
            rl += 0.048
        placed.append((x0, x0 + w, yy - h / 2, yy + h / 2))
        lx, ly = polar(R + 0.006, mid)
        ax.plot([lx, x - side * 0.012], [ly, yy], color="#AAAAAA", lw=0.35, zorder=2)
        ax.text(x + side * 0.006, yy, txt, ha="left" if side > 0 else "right",
                va="center", fontsize=5, path_effects=HALO, zorder=4)
    for n in node_order:
        a0, a1 = arcs[n]
        if abs(a0 - a1) < 10.0:
            continue
        mid = 0.5 * (a0 + a1)
        x, yy = polar(1.035, mid)
        ang = mid % 360
        if 90 < ang < 270:
            rot, ha = ang - 180, "right"
        else:
            rot, ha = ang, "left"
        rot = ((rot + 90) % 180) - 90
        ax.text(x, yy, f"{ISO2[n]}  {int(tot[n]):,}", rotation=rot, rotation_mode="anchor",
                ha=ha, va="center", fontsize=5.6, path_effects=HALO, zorder=4)
    ax.text(0, 0.06, "32,099", ha="center", fontsize=8, fontweight="bold", color="#67000D")
    ax.text(0, -0.07, "foreign-attributed\nrecords", ha="center", va="top", fontsize=5,
            color="#444444")
    letter(ax, "a", dx=0.02, dy=0.99)

    axb = fig.add_axes([0.655, 0.545, 0.305, 0.44])
    cs = cstats
    yb = np.arange(len(cs))[::-1].astype(float)
    dom = [100.0 * (r["n_samples"] - r["misaligned_n"]) / r["n_samples"] for r in cs]
    forg = [100.0 * r["misaligned_n"] / r["n_samples"] for r in cs]
    axb.barh(yb, dom, color=GREEN, height=0.66, label="Domestic submission account")
    axb.barh(yb, forg, left=dom, color=RED, height=0.66, label="Hub account (EBI/ENA)")
    axb.set_yticks(yb)
    axb.set_yticklabels([f"{ISO2[r['sample_country']]} {r['n_samples']:,}" for r in cs], fontsize=5)
    axb.set_xlim(0, 100)
    axb.set_xticks([0, 20, 40, 60, 80, 100])
    axb.set_xlabel("Domestically stewarded (%)", fontsize=5.8)
    axb.tick_params(axis="x", labelsize=5)
    axb.spines[["top", "right"]].set_visible(False)
    axb.legend(fontsize=4.6, frameon=False, loc="upper center",
               bbox_to_anchor=(0.5, -0.22), ncol=2, handlelength=1.1,
               columnspacing=1.2, handletextpad=0.5)
    letter(axb, "b", dx=-0.36, dy=1.015)

    lags = defaultdict(list)
    for r in rows:
        d1, d2 = parse_date(r["submission_date"]), parse_date(r["collection_date"])
        if not d1 or not d2:
            continue
        lag = (d1 - d2).days
        if -60 < lag <= 1460:
            lags[sc(r["geo_location"])].append(lag)
    groups = [
        ("Independent domestic", ["USA", "FRA"], BLUE),
        ("Independent LMIC-world", ["IND", "ETH"], "#EE7733"),
        ("Hub-submitted (EBI/ENA)", ["SWE", "AUT", "CHE", "GBR", "SVN", "NOR", "DNK", "NLD"], RED),
    ]
    data, labels, cols, medians, gset = [], [], [], [], []
    for gname, cts, col in groups:
        cc = sorted(cts, key=lambda c: np.median(lags[c]) if lags[c] else 1e9)
        gset.append([])
        for c in cc:
            data.append(lags[c]); labels.append(c); cols.append(col)
            gset[-1] += lags[c]
            medians.append(float(np.median(lags[c])))
    axc = fig.add_axes([0.078, 0.095, 0.892, 0.315])
    bp = axc.boxplot(data, patch_artist=True, widths=0.55,
                     medianprops=dict(color="black", lw=0.9),
                     boxprops=dict(lw=0.5), whiskerprops=dict(lw=0.5),
                     capprops=dict(lw=0.5),
                     flierprops=dict(marker=".", markersize=1.2, alpha=0.18,
                                     markeredgecolor="none", markerfacecolor="#555555"))
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.75); patch.set_edgecolor("#333333")
    for i, m in enumerate(medians):
        axc.text(i + 1, m + 40, f"{m:.0f}", ha="center", fontsize=5.2)
    axc.set_xticklabels([f"{ISO2[c]}\nn={len(d):,}" for c, d in zip(labels, data)], fontsize=5.6)
    axc.set_ylabel("Collection-to-submission\nlag (days)", fontsize=6)
    axc.tick_params(axis="y", labelsize=5.4)
    axc.set_ylim(-40, 1660)
    axc.spines[["top", "right"]].set_visible(False)
    bounds = [0]
    for _, cts, _ in groups:
        bounds.append(bounds[-1] + len(cts))
    for j, (gname, cts, col) in enumerate(groups):
        x0, x1 = bounds[j] + 0.5, bounds[j + 1] + 0.5
        axc.plot([x0, x0, x1, x1], [1420, 1455, 1455, 1420], color="#666666", lw=0.6, clip_on=False)
        axc.text((x0 + x1) / 2, 1470, gname, ha="center", va="bottom", fontsize=5.4, color="#333333")
    H3, p3 = st.kruskal(*gset)
    # submitter-defined pools (matches manuscript locked numbers):
    # independent domestic n=932 med=203 / independent LMIC-world n=389 med=349 / hub n=78,432 med=890
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
    g_ind, g_lmic, g_hub = [], [], []
    for r in rows:
        d1, d2 = parse_date(r["submission_date"]), parse_date(r["collection_date"])
        if not d1 or not d2:
            continue
        lag = (d1 - d2).days
        if not (-60 < lag <= 1460):
            continue
        c = sc(r["geo_location"])
        if r["owner"] in ("EBI", "European Bioinformatics Institute"):
            g_hub.append(lag)
        elif OWNER_CTRY.get(r["owner"], "?") == c:
            g_ind.append(lag)
            if c in ("IND", "ETH"):
                g_lmic.append(lag)
    Hs, ps = st.kruskal(g_ind, g_lmic, g_hub)
    print("country-group KW: H=%.1f p=%.3g | submitter-pool KW: n=(%d,%d,%d) H=%.1f p=%.3g medians=%s"
          % (H3, p3, len(g_ind), len(g_lmic), len(g_hub), Hs, ps,
             [float(np.median(g)) for g in (g_ind, g_lmic, g_hub)]))
    axc.text(0.997, 0.02,
             "Kruskal\u2013Wallis H(2) = %s, $P<10^{-300}$   (medians %d / %d / %d d; submitter-defined pools)"
             % (f"{Hs:,.1f}", np.median(g_ind), np.median(g_lmic), np.median(g_hub)),
             transform=axc.transAxes, ha="right", va="bottom", fontsize=5.4)
    letter(axc, "c", dx=-0.04, dy=1.10)
    save(fig, "NM_Fig2")
    print("fig2 done in %.1fs" % (time.time() - t0))

# ================================================================= Fig 3
def fig3():
    t0 = time.time()
    fig = plt.figure(figsize=(7.09, 5.05))
    axq = fig.add_axes([0.075, 0.395, 0.535, 0.575])
    axs = fig.add_axes([0.700, 0.395, 0.285, 0.575])
    axc = fig.add_axes([0.075, 0.075, 0.910, 0.190])

    qrec = defaultdict(lambda: {"hub": 0, "other": 0})
    for r in rows:
        d = parse_date(r["submission_date"])
        if not d:
            continue
        q = f"{d.year}Q{(d.month - 1) // 3 + 1}"
        is_hub = r["owner"] in ("EBI", "European Bioinformatics Institute")
        qrec[q]["hub" if is_hub else "other"] += 1
    allq = [f"{y}Q{q}" for y in range(2020, 2027) for q in range(1, 5)]
    allq = [q for q in allq if "2020Q1" <= q <= "2026Q2"]
    hub = np.array([qrec[q]["hub"] for q in allq], dtype=float)
    oth = np.array([qrec[q]["other"] for q in allq], dtype=float)
    x = np.arange(len(allq))
    axq.bar(x, hub, color=RED, width=0.78, label="Hub-submitted (EBI)")
    axq.bar(x, oth, bottom=hub, color=BLUE, width=0.78, label="Independently submitted")
    axq.set_xticks(x)
    axq.set_xticklabels(allq, rotation=90, fontsize=4.6)
    axq.set_ylabel("Records submitted per quarter", fontsize=6)
    axq.tick_params(axis="y", labelsize=5.4)
    axq.set_xlim(-0.8, len(allq) - 0.2)
    axq.spines[["top", "right"]].set_visible(False)
    top3 = np.sort(hub)[-3:]
    hubtot = hub.sum()
    zero_hub = int((hub == 0).sum())
    for xi, (h, o) in enumerate(zip(hub, oth)):
        if h >= 17000:
            axq.text(xi, h + 700, f"{int(h):,}", ha="center", fontsize=4.8, color="#67000D")
    axq.text(0.025, 0.965,
             f"three bulk waves = {100*top3.sum()/hubtot:.1f}% of hub submissions\n"
             f"hub released nothing in {zero_hub} of {len(allq)} quarters",
             transform=axq.transAxes, va="top", ha="left", fontsize=5.2, color="#333333",
             bbox=dict(facecolor="white", edgecolor="#CCCCCC", lw=0.4, pad=1.6))
    axq.legend(loc="center right", fontsize=5.2, frameon=True, framealpha=0.93,
               edgecolor="#CCCCCC", borderpad=0.4, bbox_to_anchor=(1.0, 0.42))
    letter(axq, "a", dx=-0.058, dy=1.02)
    zero_q = [q for q in allq if qrec[q]["hub"] == 0]
    print("top3 waves:", dict(zip(["2024Q3", "2022Q3", "2025Q4"], sorted(hub)[-3:])),
          "share=%.1f%%" % (100 * top3.sum() / hubtot), "zero-q:", len(zero_q))

    # b: burden scatter
    INCOME = {"GBR": "H", "NLD": "H", "DNK": "H", "AUT": "H", "CHE": "H", "SVN": "H",
              "USA": "H", "SWE": "H", "NOR": "H", "FRA": "H", "LIE": "H", "ITA": "H",
              "DEU": "H", "URY": "H", "CHL": "H", "ESP": "H", "CAN": "H",
              "IND": "L", "ETH": "L", "PER": "L"}
    for b in burden:
        c = RED if INCOME[b["iso"]] == "H" else "#762A83"
        axs.scatter(b["total_deaths"], b["records"], s=14, c=c, alpha=0.88,
                    edgecolor="white", linewidth=0.3, zorder=3)
    for iso, dx, dy in [("GBR", 0.70, 1.5), ("USA", 0.80, 1.35), ("IND", 0.78, 1.45),
                        ("PER", 0.72, 0.62), ("ETH", 1.22, 0.72)]:
        if iso in BUR:
            axs.text(BUR[iso]["total_deaths"] * dx, BUR[iso]["records"] * dy, NAME[iso],
                     fontsize=5.2, color="#111111", style="italic")
    axs.set_xscale("log"); axs.set_yscale("log")
    axs.set_xlim(5e3, 3e6); axs.set_ylim(0.55, 2.4e5)
    axs.set_xlabel("Cumulative COVID-19 deaths\n(WHO, log)", fontsize=5.8)
    axs.set_ylabel("Wastewater records (log)", fontsize=5.8)
    axs.tick_params(labelsize=5.2)
    axs.spines[["top", "right"]].set_visible(False)
    axs.text(0.04, 0.955, "Spearman \u03c1 = \u22120.08\nP = 0.75 (n = 20)", transform=axs.transAxes,
             fontsize=5.2, va="top",
             bbox=dict(boxstyle="round,pad=0.24", fc="white", ec="#BBBBBB", lw=0.4))
    from matplotlib.lines import Line2D
    axs.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor=RED,
                               markersize=4, label="High-income"),
                        Line2D([0], [0], marker="o", color="w", markerfacecolor="#762A83",
                               markersize=4, label="LMIC/UMIC")],
               fontsize=5, loc="lower left", frameon=False, handletextpad=0.2)
    letter(axs, "b", dx=-0.26, dy=1.02)

    # c: records per 100k deaths (bottom strip, full width)
    sel = ["DNK", "NLD", "LIE", "GBR", "ETH", "NOR", "SWE", "FRA", "IND", "USA", "ESP", "PER"]
    sel = [s for s in sel if s in BUR]
    ys = np.arange(len(sel))[::-1].astype(float)
    vv = [BUR[s]["records_per_100k_deaths"] for s in sel]
    axc.barh(ys, vv, color=[RED if s in ("GBR", "NLD", "DNK") else
                            "#762A83" if s in ("IND", "ETH", "PER") else BLUE for s in sel],
             height=0.62)
    axc.set_yticks(ys); axc.set_yticklabels([NAME[s] for s in sel], fontsize=5.2)
    axc.set_xscale("log"); axc.set_xlim(0.2, 2e6)
    for yi, v in zip(ys, vv):
        axc.text(v * 1.4, yi, f"{v:,.0f}" if v >= 10 else f"{v:.2f}", va="center", fontsize=4.9)
    axc.set_xlabel("Records per 100,000 reported COVID-19 deaths (log)", fontsize=5.8)
    axc.tick_params(axis="x", labelsize=5)
    axc.set_xticks([1, 10, 100, 1000, 10000, 100000])
    axc.spines[["top", "right"]].set_visible(False)
    letter(axc, "c", dx=-0.050, dy=1.16)
    save(fig, "NM_Fig3")
    print("fig3 done in %.1fs" % (time.time() - t0))

# ================================================================= Fig 4 (data-flow architecture)
def fig4():
    """Icon-based data-flow architecture with audit metrics and O/V/E badges."""
    t0 = time.time()
    fig = plt.figure(figsize=(7.09, 5.60))
    ax = fig.add_axes([0.005, 0.005, 0.99, 0.99])
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

    BROWN = "#8C510A"; DGREEN = "#1B7837"; CRIM = "#A50F15"; BLU = "#2171B5"
    COL_CLOSED = "#F5EFE6"; COL_OPEN = "#EFF6FB"
    BADGE = {"g": "#238B45", "a": "#E69F00", "r": "#CB181D", "n": "#B3B3B3"}

    rows = [
        dict(num="1", name="Poliovirus\nsurveillance", col=BROWN,
             icons=["drop", "flask", "gov", "folderlock"],
             labels=["Sewershed", "National\nlaboratory", "WHO / national\nstructures", "Primary records\n(not deposited)"],
             fills=[COL_CLOSED] * 4,
             metrics="0 records in the public\narchive; primary records\nnot openly deposited",
             outcome="Gated", badges="rnn"),
        dict(num="2", name="China CWSS\nplatform", col=BROWN,
             icons=["drop", "flask", "server", "redx"],
             labels=["Sewershed", "City / national\nlaboratory", "National platform\n(internal analysis)", "No public\nport"],
             fills=[COL_CLOSED] * 4,
             metrics="0 records in the public\narchive; >200,000 data\npoints held nationally",
             outcome="Closed", badges="rnn"),
        dict(num="3", name="Open INSDC\n(SARS-CoV-2)", col=CRIM,
             icons=["drop", "flask", "globe", "dbopen"],
             labels=["Sewershed", "National\nlaboratory", "EBI / ENA gateway\n(UK-registered;\n98.9%)", "Public archive"],
             fills=[COL_OPEN] * 4,
             metrics="82,757 records; 98.9%\nvia one account; median\nlag 890 d; LMIC 0.47%",
             outcome="Open but routed", badges="grr", highlight=True),
        dict(num="4", name="GISAID\nrepository", col=BLU,
             icons=["drop", "flask", "doclock", "dbkey"],
             labels=["Sewershed", "Laboratory", "Registration +\ndata-access agreement", "GISAID\n(approval)"],
             fills=[COL_CLOSED] * 4,
             metrics="Few wastewater records;\napproval-gated access",
             outcome="Gated", badges="ann"),
        dict(num="5", name="Proposed\nregional hub", col=DGREEN,
             icons=["drop", "flask", "house", "dbcheck"],
             labels=["Sewershed", "Regional\nsequencing hub", "Local submission\naccount + attribution\nmetadata", "Public archive"],
             fills=["#E8F4E8"] * 4,
             metrics="Local account holding;\nmodel: Ethiopia, 79\nrecords, median 672 d",
             outcome="Open and sovereign", badges="ggg"),
    ]

    row_y = {0: 87.0, 1: 70.6, 2: 54.2, 3: 37.8, 4: 19.4}
    row_h = 12.8
    widths = [12.5, 12.5, 15.5, 13.5]
    xs_left = [13.0, 28.5, 45.5, 66.5]
    CB_X, CB_W = 82.0, 17.6          # callout box

    # ---- icon primitives (line art, centred at cx, cy, size s ~ half-width)
    def ic_drop(cx, cy, s, c):
        th = np.deg2rad(np.linspace(130, 410, 60))
        r = s * 0.55
        cyc = cy - s * 0.18
        xs = cx + r * np.cos(th)
        ys = cyc + r * np.sin(th)
        ax.add_patch(plt.Polygon(np.column_stack([np.r_[xs, cx], np.r_[ys, cy + s * 0.85]]),
                                 closed=True, fill=False, edgecolor=c, lw=0.9, zorder=4,
                                 joinstyle="round"))

    def ic_flask(cx, cy, s, c):
        ax.plot([cx - s * 0.28, cx - s * 0.28, cx - s * 0.75, cx + s * 0.75,
                 cx + s * 0.28, cx + s * 0.28, cx - s * 0.28],
                [cy + s * 0.8, cy + s * 0.15, cy - s * 0.7, cy - s * 0.7,
                 cy + s * 0.15, cy + s * 0.8, cy + s * 0.8],
                color=c, lw=0.9, zorder=4)
        ax.plot([cx - s * 0.5, cx + s * 0.5], [cy - s * 0.28, cy - s * 0.28],
                color=c, lw=0.9, zorder=4)

    def ic_gov(cx, cy, s, c):
        ax.add_patch(plt.Polygon([(cx - s * 0.85, cy + s * 0.25), (cx, cy + s * 0.85),
                                  (cx + s * 0.85, cy + s * 0.25)], closed=True,
                                 fill=False, edgecolor=c, lw=0.9, zorder=4))
        for k in (-0.55, -0.18, 0.18, 0.55):
            ax.plot([cx + k * s, cx + k * s], [cy + s * 0.2, cy - s * 0.55],
                    color=c, lw=0.9, zorder=4)
        ax.plot([cx - s * 0.85, cx + s * 0.85], [cy - s * 0.7, cy - s * 0.7],
                color=c, lw=0.9, zorder=4)

    def ic_server(cx, cy, s, c):
        ax.add_patch(plt.Rectangle((cx - s * 0.65, cy - s * 0.8), s * 1.3, s * 1.6,
                                   fill=False, edgecolor=c, lw=0.9, zorder=4))
        for yy in (cy + s * 0.42, cy):
            ax.plot([cx - s * 0.65, cx + s * 0.65], [yy, yy], color=c, lw=0.9, zorder=4)
        ax.plot([cx - s * 0.35], [cy - s * 0.38], marker="o", ms=1.3, color=c, zorder=4)
        ax.plot([cx + 0.05 * s, cx + 0.4 * s], [cy - s * 0.38] * 2, color=c, lw=0.9, zorder=4)

    def ic_globe(cx, cy, s, c):
        ax.add_patch(plt.Circle((cx, cy), s * 0.78, fill=False, edgecolor=c, lw=0.9, zorder=4))
        ax.add_patch(mpatches.Ellipse((cx, cy), s * 0.75, s * 1.56, fill=False,
                                 edgecolor=c, lw=0.7, zorder=4))
        ax.plot([cx - s * 0.78, cx + s * 0.78], [cy, cy], color=c, lw=0.7, zorder=4)

    def ic_doclock(cx, cy, s, c):
        ax.plot([cx - s * 0.45, cx - s * 0.45, cx + s * 0.05, cx + s * 0.5,
                 cx + s * 0.5, cx - s * 0.45],
                [cy - s * 0.75, cy + s * 0.55, cy + s * 0.55, cy + s * 0.1,
                 cy - s * 0.75, cy - s * 0.75], color=c, lw=0.9, zorder=4)
        ax.plot([cx + 0.05 * s, cx + 0.5 * s], [cy + s * 0.55, cy + s * 0.1],
                color=c, lw=0.7, zorder=4)
        ax.add_patch(plt.Rectangle((cx - 0.02 * s, cy - s * 0.45), s * 0.34, s * 0.3,
                                   fill=False, edgecolor=c, lw=0.8, zorder=4))
        ax.add_patch(mpatches.Arc((cx + 0.15 * s, cy - s * 0.15), s * 0.24, s * 0.24,
                             theta1=0, theta2=180, edgecolor=c, lw=0.8, zorder=4))

    def ic_house(cx, cy, s, c):
        ax.add_patch(plt.Polygon([(cx - s * 0.75, cy - s * 0.75), (cx - s * 0.75, cy + s * 0.05),
                                  (cx, cy + s * 0.75), (cx + s * 0.75, cy + s * 0.05),
                                  (cx + s * 0.75, cy - s * 0.75)], closed=True,
                                 fill=False, edgecolor=c, lw=0.9, zorder=4))
        ax.plot([cx], [cy - s * 0.05], marker="o", ms=2.0, color=c, zorder=4)
        for ang in (45, 135, 225, 315):
            dx = 0.42 * s * np.cos(np.deg2rad(ang)); dy = 0.42 * s * np.sin(np.deg2rad(ang))
            ax.plot([cx + dx * 0.45, cx + dx], [cy - s * 0.05 + dy * 0.45, cy - s * 0.05 + dy],
                    color=c, lw=0.7, zorder=4)

    def ic_folderlock(cx, cy, s, c):
        ax.plot([cx - s * 0.75, cx - s * 0.75, cx - s * 0.3, cx - s * 0.15, cx + s * 0.75,
                 cx + s * 0.75, cx - s * 0.75],
                [cy - s * 0.55, cy + s * 0.55, cy + s * 0.55, cy + s * 0.35, cy + s * 0.35,
                 cy - s * 0.55, cy - s * 0.55], color=c, lw=0.9, zorder=4)
        ax.add_patch(plt.Rectangle((cx - s * 0.18, cy - s * 0.35), s * 0.4, s * 0.34,
                                   fill=False, edgecolor=c, lw=0.8, zorder=4))
        ax.add_patch(mpatches.Arc((cx, cy - s * 0.01), s * 0.26, s * 0.26,
                             theta1=0, theta2=180, edgecolor=c, lw=0.8, zorder=4))

    def ic_redx(cx, cy, s, c):
        ax.plot([cx - s * 0.6, cx + s * 0.6], [cy - s * 0.6, cy + s * 0.6],
                color=CRIM, lw=1.6, zorder=4, solid_capstyle="round")
        ax.plot([cx - s * 0.6, cx + s * 0.6], [cy + s * 0.6, cy - s * 0.6],
                color=CRIM, lw=1.6, zorder=4, solid_capstyle="round")

    def _cyl(cx, cy, s, c):
        ax.add_patch(plt.Rectangle((cx - s * 0.6, cy - s * 0.5), s * 1.2, s * 1.0,
                                   fill=False, edgecolor=c, lw=0.9, zorder=4))
        ax.add_patch(mpatches.Ellipse((cx, cy + s * 0.5), s * 1.2, s * 0.42,
                                 fill=False, edgecolor=c, lw=0.9, zorder=4))
        ax.add_patch(mpatches.Arc((cx, cy - s * 0.5), s * 1.2, s * 0.42,
                             theta1=180, theta2=360, edgecolor=c, lw=0.9, zorder=4))

    def ic_dbopen(cx, cy, s, c):
        _cyl(cx, cy, s, c)
        for yy in (cy + s * 0.1, cy - s * 0.2):
            ax.add_patch(mpatches.Arc((cx, yy), s * 1.2, s * 0.42, theta1=20,
                                 theta2=160, edgecolor=c, lw=0.6, zorder=4))

    def ic_dbkey(cx, cy, s, c):
        _cyl(cx, cy, s, c)
        ax.add_patch(plt.Rectangle((cx + s * 0.15, cy - s * 0.62), s * 0.42, s * 0.34,
                                   fill=False, edgecolor=c, lw=0.8, zorder=5))
        ax.add_patch(mpatches.Arc((cx + s * 0.36, cy - s * 0.28), s * 0.3, s * 0.3,
                             theta1=0, theta2=180, edgecolor=c, lw=0.8, zorder=5))

    def ic_dbcheck(cx, cy, s, c):
        _cyl(cx, cy, s, c)
        ax.plot([cx - s * 0.35, cx - s * 0.08, cx + s * 0.42],
                [cy - s * 0.05, cy - s * 0.38, cy + s * 0.25],
                color=c, lw=1.2, zorder=5, solid_capstyle="round")

    ICON_FN = {"drop": ic_drop, "flask": ic_flask, "gov": ic_gov, "server": ic_server,
               "globe": ic_globe, "doclock": ic_doclock, "house": ic_house,
               "folderlock": ic_folderlock, "redx": ic_redx, "dbopen": ic_dbopen,
               "dbkey": ic_dbkey, "dbcheck": ic_dbcheck}

    for ri, m in enumerate(rows):
        yc = row_y[ri]
        y0 = yc - row_h / 2
        if ri == 4:
            sep = y0 + row_h + 2.9
            ax.plot([0.5, 99.5], [sep, sep], color=DGREEN, lw=0.8, linestyle=(0, (4, 2)))
            ax.text(99.0, sep + 0.8, "proposed architecture", fontsize=5.4,
                    color=DGREEN, ha="right", va="bottom", style="italic")
        # number chip + name
        ax.add_patch(FancyBboxPatch((1.0, y0 + 2.0), 3.6, row_h - 4.0,
                                    boxstyle="round,pad=0.3,rounding_size=0.9",
                                    fc=m["col"], ec=m["col"], lw=0.8))
        ax.text(2.8, yc, m["num"], fontsize=6.5, fontweight="bold", color="white",
                ha="center", va="center", zorder=5)
        ax.text(5.6, yc + 1.0, m["name"], fontsize=5.1, va="center", ha="left",
                color="#222222", linespacing=1.35)
        # stage nodes: icon above label
        for ni, (lab, fc) in enumerate(zip(m["labels"], m["fills"])):
            x0 = xs_left[ni]; w = widths[ni]
            lw = 1.2 if (ri == 2 and ni == 2) else 0.7
            ax.add_patch(FancyBboxPatch((x0, y0), w, row_h,
                                        boxstyle="round,pad=0.3,rounding_size=1.1",
                                        fc=fc, ec=m["col"], lw=lw, zorder=2))
            cx = x0 + w / 2
            ICON_FN[m["icons"][ni]](cx, yc + row_h * 0.22, w * 0.36, m["col"])
            ax.text(cx, y0 + 2.5, lab, fontsize=4.1, ha="center", va="center",
                    color="#222222", zorder=4, linespacing=1.25)
            if ri == 2 and ni == 2:
                ax.text(cx, y0 - 1.9, "single national gateway", fontsize=4.4,
                        ha="center", va="top", color=CRIM, style="italic", zorder=4)
            if ni < 3:
                x1 = x0 + w
                x2 = xs_left[ni + 1]
                ax.annotate("", xy=(x2 - 0.4, yc), xytext=(x1 + 0.4, yc),
                            arrowprops=dict(arrowstyle="-|>", color="#888888", lw=0.8),
                            zorder=1)
        # right callout box
        ax.add_patch(FancyBboxPatch((CB_X, y0), CB_W, row_h,
                                    boxstyle="round,pad=0.3,rounding_size=1.1",
                                    fc="white", ec=m["col"], lw=0.8, zorder=2))
        ax.text(CB_X + CB_W / 2, y0 + row_h - 2.6, m["metrics"], fontsize=3.9,
                ha="center", va="center", color="#333333", linespacing=1.3, zorder=4)
        ax.text(CB_X + CB_W / 2, y0 + 4.1, m["outcome"], fontsize=5.0,
                ha="center", va="center", color=m["col"], fontweight="bold", zorder=4)
        # badges
        bx0 = CB_X + CB_W / 2 - 3.4
        for bi, (letter_, code) in enumerate(zip(["O", "V", "E"], m["badges"])):
            bxc = bx0 + bi * 3.4
            ax.add_patch(plt.Circle((bxc, y0 + 1.6), 1.15, fc=BADGE[code],
                                    ec="white", lw=0.5, zorder=5))
            ax.text(bxc, y0 + 1.6, letter_, fontsize=3.6, ha="center", va="center",
                    color="white", fontweight="bold", zorder=6)
        ax.annotate("", xy=(CB_X - 0.4, yc), xytext=(66.5 + widths[3] + 0.4, yc),
                    arrowprops=dict(arrowstyle="-|>", color="#888888", lw=0.8), zorder=1)

    ax.text(13.0, 97.6, "DATA-FLOW ARCHITECTURE", fontsize=5.4, color="#555555",
            va="center", style="italic")
    ax.text(99.0, 97.6, "AUDIT METRICS AND ASSESSMENT", fontsize=5.4, color="#555555",
            va="center", ha="right", style="italic")
    # bottom legend
    ly = 4.2
    ax.text(1.0, ly + 2.2, "Badges:  O = openness   V = velocity   E = equity",
            fontsize=4.6, color="#222222", va="center")
    lx = 1.0
    for code, lab in [("g", "favourable"), ("a", "intermediate"),
                      ("r", "unfavourable"), ("n", "not applicable")]:
        ax.add_patch(plt.Circle((lx, ly), 1.05, fc=BADGE[code], ec="white", lw=0.5))
        ax.text(lx + 1.7, ly, lab, fontsize=4.4, color="#333333", va="center")
        lx += 1.7 + len(lab) * 0.95 + 3.2
    save(fig, "NM_Fig4")
    print("fig4 done in %.1fs" % (time.time() - t0))

# ================================================================= SI Fig S1 (hub lag distribution + sensitivity)
def figs1():
    t0 = time.time()
    fig = plt.figure(figsize=(7.09, 3.0))
    axa = fig.add_axes([0.065, 0.17, 0.30, 0.76])
    axb = fig.add_axes([0.415, 0.17, 0.24, 0.76])
    axc = fig.add_axes([0.715, 0.17, 0.25, 0.76])

    lags = {"hub": [], "dom": []}
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
    import scipy.stats as st
    for r in rows:
        d1, d2 = parse_date(r["submission_date"]), parse_date(r["collection_date"])
        if not d1 or not d2:
            continue
        lag = (d1 - d2).days
        if not (-60 < lag <= 1460):
            continue
        c = sc(r["geo_location"])
        if r["owner"] in ("EBI", "European Bioinformatics Institute"):
            lags["hub"].append(lag)
        elif OWNER_CTRY.get(r["owner"], "?") == c:
            lags["dom"].append(lag)

    bins = np.arange(0, 1470, 40)
    axa.hist(lags["hub"], bins=bins, color=RED, alpha=0.75, label="Hub-routed (n=78,432)")
    axa.hist(lags["dom"], bins=bins, color=BLUE, alpha=0.75, label="Independent domestic (n=932)")
    axa.set_yscale("log")
    axa.set_ylim(0.5, 30000)
    axa.set_xlabel("Collection-to-submission lag (days)", fontsize=5.6)
    axa.set_ylabel("Records (log)", fontsize=5.6)
    axa.tick_params(labelsize=5)
    axa.legend(fontsize=4.6, frameon=False, loc="upper right")
    axa.spines[["top", "right"]].set_visible(False)
    letter(axa, "a", dx=-0.115, dy=1.06)

    by_year = defaultdict(list)
    for r in rows:
        if r["owner"] not in ("EBI", "European Bioinformatics Institute"):
            continue
        d1, d2 = parse_date(r["submission_date"]), parse_date(r["collection_date"])
        if not d1 or not d2:
            continue
        lag = (d1 - d2).days
        if -60 < lag <= 1460:
            by_year[d2.year].append(lag)
    yrs = sorted(by_year)
    med = [np.median(by_year[y]) for y in yrs]
    xs = np.arange(len(yrs))
    axb.bar(xs, med, color=RED, width=0.62)
    for xi, (y, m) in enumerate(zip(yrs, med)):
        axb.text(xi, m + 22, f"{m:.0f}", ha="center", fontsize=5)
        axb.text(xi, -115, f"n={len(by_year[y]):,}", ha="center", fontsize=4.4, color="#555555")
    axb.set_xticks(xs); axb.set_xticklabels([str(y) for y in yrs], fontsize=5.2)
    axb.set_ylim(-160, 1050)
    axb.set_xlabel("Collection year (hub-routed)", fontsize=5.6)
    axb.set_ylabel("Median lag (days)", fontsize=5.6)
    axb.tick_params(labelsize=5)
    axb.spines[["top", "right"]].set_visible(False)
    axb.text(0.04, 0.92, "back-filling signature", transform=axb.transAxes,
             fontsize=5, color="#67000D", style="italic")
    letter(axb, "b", dx=-0.16, dy=1.06)

    # c) sensitivity of medians
    specs = [
        ("Full\nwindow", 890, 203, r"$<10^{-300}$"),
        ("Collected\n≥ 2022", 896, 351, r"$3.9\times10^{-58}$"),
        ("Collected\n≥ 2023", 825, 678, r"$1.1\times10^{-20}$"),
        ("Lag ≤ 365 d\nonly", 263, 153, r"$5.2\times10^{-85}$"),
    ]
    xs = np.arange(len(specs))
    w = 0.36
    axc.bar(xs - w / 2, [s[1] for s in specs], width=w, color=RED, label="Hub-routed")
    axc.bar(xs + w / 2, [s[2] for s in specs], width=w, color=BLUE, label="Independent domestic")
    for xi, s in enumerate(specs):
        axc.text(xi - w / 2, s[1] + 18, f"{s[1]:,}", ha="center", fontsize=4.8)
        axc.text(xi + w / 2, s[2] + 18, f"{s[2]:,}", ha="center", fontsize=4.8)
        axc.text(xi, 990, f"P = {s[3]}", ha="center", fontsize=4.1, color="#333333")
    axc.set_xticks(xs)
    axc.set_xticklabels([s[0] for s in specs], fontsize=5.0)
    axc.set_ylim(0, 1120)
    axc.set_ylabel("Median lag (days)", fontsize=5.6)
    axc.tick_params(labelsize=5)
    axc.legend(fontsize=4.6, frameon=False, loc="lower right",
               bbox_to_anchor=(1.0, 1.005), ncol=2)
    axc.spines[["top", "right"]].set_visible(False)
    letter(axc, "c", dx=-0.16, dy=1.06)
    save(fig, "NM_FigS1")
    print("figS1 done in %.1fs" % (time.time() - t0))

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "1"): fig1()
    if which in ("all", "2"): fig2()
    if which in ("all", "3"): fig3()
    if which in ("all", "4"): fig4()
    if which in ("all", "s1"): figs1()
    print("ALL DONE")

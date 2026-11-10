# -*- coding: utf-8 -*-
"""Fast cluster bootstrap for the no-window Cox model (M minor 8).

Breslow-ties Cox partial likelihood implemented directly in numpy (Newton-Raphson),
verified against lifelines on the full frame. Country-level pairs bootstrap, B=2000.
"""
import csv, sys, json, io
from datetime import datetime
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd

sys.path.insert(0, "材料3_统计建模代码及运行结果/code")
from map_analyze import parse_date, sample_country

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
    if r["owner"] in HUB:
        pool = "hub_routed"
    elif OWNER_CTRY.get(r["owner"], "?") == c and c != "(missing)":
        pool = "indep_LMIC" if c in ("IND", "ETH") else "indep_HIC_dom"
    else:
        continue
    if not (d1 and d2):
        continue
    lag = (d1 - d2).days
    if lag <= -60:
        continue
    rows.append((max(lag, 0.5), 1, pool, d2.year, c))
f = pd.DataFrame(rows, columns=["dur", "event", "pool", "coll_year", "country"])
X = np.column_stack([
    (f.pool == "hub_routed").astype(float),
    (f.pool == "indep_LMIC").astype(float),
    (f.coll_year - 2021).astype(float)])
t = f.dur.values; ev = f.event.values; cl = f.country.values

def make_static(t, ev, X):
    """Precompute time structure (invariant across Newton iterations)."""
    order = np.argsort(-t)
    Xs, ts, evs = X[order], t[order], ev[order]
    ut = np.unique(ts)
    pos = np.searchsorted(-ts, -ut, side="right") - 1  # last index of each unique time;
    # prefix sum up to pos = sum over {ts >= u} = risk set (descending-sorted array)
    evmask = (evs == 1)
    ut_ev = ut[evmask[np.searchsorted(-ts, -ut)]]  # not used; compute below properly
    # events per unique time
    dN = np.array([evmask[np.searchsorted(-ts, -u):].sum() if False else 0 for u in ut])
    # efficient: counts via bincount over event records mapped to unique-time index
    tidx = np.searchsorted(ut, ts)
    dN = np.bincount(tidx[evmask], minlength=len(ut))
    dX = np.zeros((len(ut), X.shape[1]))
    np.add.at(dX, tidx[evmask], Xs[evmask])
    has_ev = dN > 0
    return Xs, ts, ut, pos[has_ev], dN[has_ev], dX[has_ev]

def loglik(Xs, pos, dN, dX, beta):
    w = np.exp(Xs @ beta)
    csw = np.cumsum(w)
    ll1 = float(dX.sum(axis=0) @ beta)
    with np.errstate(divide="ignore", invalid="ignore"):
        ll2 = float((dN * np.log(csw[pos])).sum())
    return ll1 - ll2

def cox_breslow(X, t, ev, beta0=None, max_iter=60, tol=1e-8):
    n, p = X.shape
    Xs, ts, ut, pos, dN, dX = make_static(t, ev, X)
    beta = np.zeros(p) if beta0 is None else beta0.copy()
    ll = loglik(Xs, pos, dN, dX, beta)
    for it in range(max_iter):
        w = np.exp(Xs @ beta)
        csw = np.cumsum(w)
        cswX = np.cumsum(w[:, None] * Xs, axis=0)
        S0 = csw[pos]; S1 = cswX[pos]
        with np.errstate(divide="ignore", invalid="ignore"):
            U = (dX - dN[:, None] * (S1 / S0[:, None])).sum(axis=0)
            S2 = np.cumsum(w[:, None, None] * Xs[:, :, None] * Xs[:, None, :], axis=0)[pos]
            V = S2 / S0[:, None, None] - (S1[:, :, None] * S1[:, None, :] / (S0[:, None, None] * S0[:, None, None]))
        I = (V * dN[:, None, None]).sum(axis=0)
        try:
            step = np.linalg.solve(I, U)
        except np.linalg.LinAlgError:
            break
        if (not np.all(np.isfinite(step))) or np.max(np.abs(step)) < tol:
            break
        accepted = False
        for _ in range(20):
            cand = beta + step
            cand_ll = loglik(Xs, pos, dN, dX, cand)
            if np.isfinite(cand_ll) and cand_ll >= ll:
                accepted = True
                break
            step = step / 2
        if not accepted:
            break
        beta = cand
        ll = cand_ll
    return beta

beta_full = cox_breslow(X, t, ev)
print("fast Breslow full-frame beta:", beta_full, "-> HR:", np.exp(beta_full))

rng = np.random.default_rng(42)
countries = np.unique(cl)
B = 2000
boot = []
for b in range(B):
    pick = rng.choice(countries, size=len(countries), replace=True)
    idx = np.concatenate([np.where(cl == c)[0] for c in pick])
    Xb, tb, evb = X[idx], t[idx], ev[idx]
    try:
        bb = cox_breslow(Xb, tb, evb, beta0=beta_full)
        boot.append(bb)
    except Exception:
        pass
boot = np.array(boot)
out = {}
for i, nm in enumerate(["pool_hub", "pool_lmic", "coll_year_c"]):
    v = np.exp(boot[:, i])
    out[nm] = {"p2.5": float(np.percentile(v, 2.5)),
               "p97.5": float(np.percentile(v, 97.5)),
               "median": float(np.median(v)), "n_valid": int(len(v))}
print(json.dumps(out, indent=1))
json.dump(out, open("cluster_bootstrap.json", "w"), indent=1)
print("saved cluster_bootstrap.json")

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""48_benchmark_prep.py —— 对标模板的数据准备（PPI/ROC/DCA/PCA/热图）

输出（results/benchmark/）:
  string_edges.csv / string_nodes.csv   STRING PPI（核心 100 基因）
  per_gene_roc.csv                      3 签名基因 ROC 曲线坐标
  per_gene_auc.csv                      3 签名基因 AUC
  dca.csv                               列线图决策曲线（净获益）
  subtype_pca.csv                       核心基因 PCA（IgAN 样本）
  subtype_heatmap.csv                   核心基因 z-score（IgAN 样本 x 基因，按亚型排序）
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import requests
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, roc_curve

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"E:/sheng xin/ObstructiveNephropathy_MRG"
PROC = os.path.join(ROOT, "data", "processed")
OUT = os.path.join(ROOT, "results", "benchmark")
os.makedirs(OUT, exist_ok=True)

meta = pd.read_csv(os.path.join(PROC, "meta_human_main.csv"))
sig = pd.read_csv(os.path.join(OUT, "sig_expr.csv"))
core = pd.read_csv(os.path.join(OUT, "core_expr.csv"))
cal = pd.read_csv(os.path.join(PROC, "panel_export", "calibration.csv"))
sub = pd.read_csv(os.path.join(PROC, "subtype_assignment.csv"))

ctrl = set(meta.loc[meta["group"] == "other", "sample"])
iga = set(meta.loc[meta["group"].str.startswith("IgAN"), "sample"])


# ---- 1) 逐基因 ROC ----
def per_gene_roc():
    rows = []
    aucs = []
    for gene in ["NDNF", "PCDHB7", "RRAGB"]:
        g = sig[sig.gene == gene].iloc[0]
        vals = {c: v for c, v in g.items() if c != "gene"}
        # 只保留 IgAN 与 30 对照（排 FSGS）
        y, x = [], []
        for s, v in vals.items():
            if s in iga:
                y.append(1)
                x.append(v)
            elif s in ctrl:
                y.append(0)
                x.append(v)
        y = np.array(y)
        x = np.array(x)
        direction = "up" if roc_auc_score(y, x) >= 0.5 else "down"
        if direction == "down":
            x = -x
        auc = roc_auc_score(y, x)
        aucs.append({"gene": gene, "AUC": auc,
                     "direction_in_disease": direction,
                     "n_iga": int(y.sum()),
                     "n_ctrl": int(len(y) - y.sum())})
        fpr, tpr, _ = roc_curve(y, x)
        rows.append(pd.DataFrame({"gene": gene, "fpr": fpr, "tpr": tpr}))
    pd.concat(rows).to_csv(os.path.join(OUT, "per_gene_roc.csv"), index=False)
    pd.DataFrame(aucs).to_csv(os.path.join(OUT, "per_gene_auc.csv"), index=False)
    print("per-gene AUC:", {a["gene"]: round(a["AUC"], 3) for a in aucs})


# ---- 2) DCA ----
def dca():
    p = cal.pred.to_numpy()
    y = cal.obs.to_numpy()
    pts = np.linspace(0.01, 0.99, 99)
    n = len(y)
    rows = []
    for pt in pts:
        pred_pos = p >= pt
        tp = int(((pred_pos) & (y == 1)).sum())
        fp = int(((pred_pos) & (y == 0)).sum())
        nb = (tp - fp * pt / (1 - pt)) / n
        # 全部治疗/全部不治疗
        nb_all = (y.sum() - (n - y.sum()) * pt / (1 - pt)) / n
        rows.append({"threshold": pt, "net_benefit_model": nb,
                     "net_benefit_all": nb_all, "net_benefit_none": 0.0})
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "dca.csv"), index=False)
    print("DCA rows:", len(rows))


# ---- 3) 亚型 PCA 与热图 ----
def subtype_panels():
    submap = dict(zip(sub["sample"], sub["subtype"]))
    iga_meta = meta[meta["group"].str.startswith("IgAN")]
    samples = [s for s in iga_meta["sample"] if s in submap]
    core_sub = core[["gene"] + samples]
    X = core_sub[samples].T.to_numpy()
    Xs = (X - X.mean(axis=0)) / X.std(axis=0)
    pca = PCA(n_components=2, random_state=42)
    P = pca.fit_transform(Xs)
    pca_df = pd.DataFrame({"sample": samples, "PC1": P[:, 0], "PC2": P[:, 1],
                           "subtype": [submap[s] for s in samples]})
    pca_df.to_csv(os.path.join(OUT, "subtype_pca.csv"), index=False)
    print("PCA var ratio:", [round(v, 3) for v in pca.explained_variance_ratio_])

    # 热图：按亚型排序样本，行为基因（按各亚型间 t 值排序取 top 60）
    sub_order = pca_df.sort_values("subtype")
    cols = list(sub_order["sample"])
    z = core_sub.set_index("gene")[cols].T
    z = ((z - z.mean()) / z.std()).T
    tvals = []
    for gene in z.index:
        a = z.loc[gene, sub_order["sample"][sub_order["subtype"] == "C1"]].to_numpy()
        b = z.loc[gene, sub_order["sample"][sub_order["subtype"] == "C2"]].to_numpy()
        se = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
        tvals.append((gene, abs(a.mean() - b.mean()) / se if se > 0 else 0))
    tvals.sort(key=lambda r: -r[1])
    keep = [g for g, _ in tvals[:60]]
    hm = z.loc[keep, cols]
    hm.insert(0, "gene", hm.index)
    hm.to_csv(os.path.join(OUT, "subtype_heatmap.csv"), index=False)
    anno = pd.DataFrame({"sample": cols,
                         "subtype": [submap[s] for s in cols]})
    anno.to_csv(os.path.join(OUT, "subtype_heatmap_anno.csv"), index=False)
    print("heatmap genes:", len(keep))


# ---- 4) STRING PPI ----
def string_ppi():
    genes = core.gene.tolist()
    url = "https://string-db.org/api/tsv/network"
    try:
        r = requests.post(url, data={"identifiers": "\r".join(genes),
                                     "species": 9606,
                                     "required_score": 400},
                          timeout=90)
        r.raise_for_status()
        df = pd.DataFrame([ln.split("\t") for ln in
                           r.text.strip().splitlines()[1:]],
                          columns=r.text.splitlines()[0].split("\t"))
        df.to_csv(os.path.join(OUT, "string_edges.csv"), index=False)
        deg = pd.concat([df.preferredName_A, df.preferredName_B]
                        ).value_counts().rename("degree")
        nodes = pd.DataFrame({"gene": deg.index, "degree": deg.values})
        nodes.to_csv(os.path.join(OUT, "string_nodes.csv"), index=False)
        print("STRING edges:", len(df), "nodes:", len(nodes),
              "top hubs:", nodes.head(8).gene.tolist())
    except Exception as e:
        print("STRING failed:", repr(e))
        pd.DataFrame().to_csv(os.path.join(OUT, "string_edges.csv"))


per_gene_roc()
dca()
subtype_panels()
string_ppi()
print("DONE")

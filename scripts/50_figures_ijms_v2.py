#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""50_figures_ijms_v2.py —— 对标模板范文的 IJMS 7 图版（Fig1–Fig6）

Fig7（qPCR，仅描述性）由 42_qpcr_ijms.py 生成。
规格同 43：178 mm 双栏、Arial 7 pt、PDF 矢量 + 300 dpi PNG/TIFF。
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from PIL import Image
from scipy import stats as sps

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"E:/sheng xin/ObstructiveNephropathy_MRG"
PROC = os.path.join(ROOT, "data", "processed")
PANEL = os.path.join(PROC, "panel_export")
BENCH = os.path.join(ROOT, "results", "benchmark")
OUTDIR = os.path.join(ROOT, "submission", "ijms", "figures")
PREVDIR = os.path.join(OUTDIR, "preview")
FIG_W_IN = 7.0
DPI = 300

CB = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"]
BLUE, ORANGE = CB[0], CB[1]
GREY = "#666666"
TOL_MUTED = ["#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
             "#DDCC77", "#CC6677", "#882255", "#AA4499", "#6699CC"]
GRIDSPEC_ORIG = Figure.add_gridspec


def log(msg):
    print(f"[figv2] {msg}", flush=True)


def setup_fonts():
    for fam in ("Arial",):
        if any(f.name == fam for f in font_manager.fontManager.ttflist):
            plt.rcParams["font.family"] = fam
            break
    plt.rcParams.update({
        "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.linewidth": 0.6, "figure.dpi": DPI, "savefig.dpi": DPI,
        "axes.edgecolor": "#333333", "xtick.color": "#333333",
        "ytick.color": "#333333", "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def patch_left(delta):
    def wrapped(self, nrows=1, ncols=1, **kw):
        if "left" in kw and "right" in kw:
            kw["left"] = kw["left"] + delta
        return GRIDSPEC_ORIG(self, nrows, ncols, **kw)
    Figure.add_gridspec = wrapped


def panel_label(ax, letter, title=None):
    ax.text(0.0, 1.008, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom", ha="left", color="black")
    if title:
        ax.text(0.055, 1.008, title, transform=ax.transAxes, fontsize=8,
                va="bottom", ha="left", color="#222222")


def style_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=0.6, length=2.5)


def main_df():
    df = pd.read_csv(os.path.join(PANEL, "main_panel_df.csv"))
    ctrl = df[(df.disease == "control")
              & (df.status != "Focal Segmental Glomerulosclerosis")]
    iga = df[df.disease == "IgAN"]
    return df, ctrl, iga


def box_jitter(ax, vals, xpos, color, jitter=0.12, width=0.42):
    ax.boxplot(vals, positions=[xpos], widths=width, patch_artist=True,
               showfliers=False, zorder=2,
               boxprops=dict(facecolor=color, alpha=0.75, lw=0.6),
               medianprops=dict(color="black", lw=0.8),
               whiskerprops=dict(color=color, lw=0.6),
               capprops=dict(color=color, lw=0.6))
    rng = np.random.default_rng(42)
    jx = xpos + rng.uniform(-jitter, jitter, len(vals))
    ax.scatter(jx, vals, s=4, alpha=0.45, color=color, linewidths=0, zorder=3)


# ---------------- Fig1 ----------------
def workflow_axis(ax):
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    boxes = [
        (0.25, 1.2, 2.7, 8.6, "#EAF2F8", BLUE, "Public transcriptomes", [
            ["Human biopsies", "GSE115857 (n = 86)"],
            ["Mouse scRNA-seq", "GSE175412 (3 kidneys)"],
            ["Human spatial CosMx", "GSE282059 (523,855 cells)"],
            ["Mouse UUO time course", "GSE118339 (n = 15)"]]),
        (3.65, 1.2, 2.7, 8.6, "#FFF6E6", ORANGE, "Mechanical framework", [
            ["2,536 mechanical genes", "(3 GO categories)"],
            ["ssGSEA mechanical score"],
            ["WGCNA β = 4: 38 modules", "100-gene core (ME31)"],
            ["LASSO ∩ random forest", "NDNF, PCDHB7, RRAGB"]]),
        (7.05, 1.2, 2.7, 8.6, "#EDF7F1", CB[2], "Outcomes", [
            ["IgAN vs controls", "P = 7.01 × 10$^{-10}$"],
            ["ROC AUC 0.853", "C-index 0.848"],
            ["Subtypes C1/C2", "P = 0.0105"],
            ["CD hub: spatial ρ = 0.586", "time course ρ = 0.94"]]),
    ]
    pairs = []
    for x, y, w, h, face, edge, title, items in boxes:
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                                    boxstyle="round,pad=0.05,rounding_size=0.15",
                                    facecolor=face, edgecolor=edge, lw=1.0))
        t = ax.text(x + w / 2, y + h - 0.42, title, ha="center", va="top",
                    fontsize=7.5, fontweight="bold", color=edge)
        pairs.append((t, (x, y, w, h), ax))
        ty = y + h - 1.20
        for item in items:
            for line in item:
                t = ax.text(x + w / 2, ty, line, ha="center", va="top",
                            fontsize=7.0, color="#222222", linespacing=1.1)
                pairs.append((t, (x, y, w, h), ax))
                ty -= 0.72
            ty -= 0.26
    for xa in (2.95, 6.35):
        ax.add_patch(FancyArrowPatch((xa, 5.0), (xa + 0.68, 5.0),
                                     arrowstyle="-|>", mutation_scale=13,
                                     color="#333333", lw=1.0))
    return pairs


def volcano_axis(ax):
    deg = pd.read_csv(os.path.join(PROC, "deg_igag_vs_control.csv"))
    mrg = set(open(os.path.join(ROOT, "submission",
                                "ijms",
                                "supporting_information",
                                "S1_Table_MRG_gene_set.txt"),
                   encoding="utf-8").read().splitlines())
    deg["is_mrg"] = deg.gene.isin(mrg)
    up = (deg.logFC > 0.585) & (deg["adj.P.Val"] < 0.05)
    dn = (deg.logFC < -0.585) & (deg["adj.P.Val"] < 0.05)
    ns = ~(up | dn)
    ax.scatter(deg.loc[ns, "logFC"], -np.log10(deg.loc[ns, "adj.P.Val"]),
               s=4, color="#BBBBBB", lw=0, rasterized=True, zorder=2)
    ax.scatter(deg.loc[up, "logFC"], -np.log10(deg.loc[up, "adj.P.Val"]),
               s=5, color=ORANGE, lw=0, rasterized=True, zorder=3)
    ax.scatter(deg.loc[dn, "logFC"], -np.log10(deg.loc[dn, "adj.P.Val"]),
               s=5, color=BLUE, lw=0, rasterized=True, zorder=3)
    mr = deg[deg.is_mrg & (up | dn)]
    ax.scatter(mr.logFC, -np.log10(mr["adj.P.Val"]), s=8, facecolors="none",
               edgecolors="#111111", lw=0.6, zorder=4)
    label_off = {"NDNF": (0.08, 0.15), "PCDHB7": (0.08, 0.15),
                 "RRAGB": (0.08, 0.15)}
    for _, r in mr[mr.gene.isin(["NDNF", "PCDHB7", "RRAGB"])].iterrows():
        dx, dy = label_off[r.gene]
        ax.text(r.logFC + dx, -np.log10(r["adj.P.Val"]) + dy, r.gene,
                fontsize=6.5, ha="left", va="bottom", color="#111111",
                fontweight="bold")
    ax.axvline(0.585, color=GREY, lw=0.6, ls="--")
    ax.axvline(-0.585, color=GREY, lw=0.6, ls="--")
    ax.axhline(-np.log10(0.05), color=GREY, lw=0.6, ls="--")
    ax.set_xlabel("log2 fold change (IgAN vs. control)")
    ax.set_ylabel("−log10(adjusted P)")
    ax.set_xlim(deg.logFC.min() - 0.12, deg.logFC.max() + 0.12)
    ax.text(0.985, 0.97, f"58 up / 173 down\n36 MRG genes\n"
            "(open circles)",
            transform=ax.transAxes, fontsize=7, va="top", ha="right")
    style_axis(ax)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, prune="both"))


def make_fig1():
    df, ctrl, iga = main_df()
    fig = plt.figure(figsize=(FIG_W_IN, 6.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 0.85],
                          width_ratios=[1.0, 1.15],
                          left=0.08, right=0.96, top=0.92, bottom=0.07,
                          hspace=0.5, wspace=0.42)
    axa = fig.add_subplot(gs[0, :])
    pairs = workflow_axis(axa)
    panel_label(axa, "a")
    axb = fig.add_subplot(gs[1, 0])
    box_jitter(axb, ctrl.MRG_up.to_numpy(), 1, BLUE)
    box_jitter(axb, iga.MRG_up.to_numpy(), 2, ORANGE)
    axb.set_xticks([1, 2])
    axb.set_xticklabels(["Control\n(n = 30)", "IgAN\n(n = 55)"], fontsize=7)
    axb.set_xlim(0.45, 2.55)
    axb.set_ylabel("Mechanical score (ssGSEA, MRG-up)")
    ymax = max(ctrl.MRG_up.max(), iga.MRG_up.max())
    axb.set_ylim(min(ctrl.MRG_up.min(), iga.MRG_up.min()) - 0.15, ymax + 0.38)
    axb.plot([1, 1, 2, 2], [ymax + 0.09, ymax + 0.20, ymax + 0.20,
                            ymax + 0.09], color="black", lw=0.7)
    axb.text(1.5, ymax + 0.25, "P = 7.01 × 10$^{-10}$", ha="center",
             va="bottom", fontsize=7)
    style_axis(axb)
    panel_label(axb, "b")
    axc = fig.add_subplot(gs[1, 1])
    volcano_axis(axc)
    panel_label(axc, "c")
    return fig, pairs


# ---------------- Fig2 ----------------
def roc_panel(ax):
    roc = pd.read_csv(os.path.join(PANEL, "roc_curve.csv"))
    per = pd.read_csv(os.path.join(BENCH, "per_gene_roc.csv"))
    aucs = pd.read_csv(os.path.join(BENCH, "per_gene_auc.csv"))
    ax.plot(roc.spec, roc.sens, color=BLUE, lw=1.5, zorder=4,
            label="3-gene signature (AUC = 0.853)")
    cols = {"NDNF": CB[3], "PCDHB7": CB[1], "RRAGB": CB[2]}
    for gene, c in cols.items():
        d = per[per.gene == gene]
        a = aucs.loc[aucs.gene == gene, "AUC"].iloc[0]
        ax.plot(d.fpr, d.tpr, color=c, lw=0.9, zorder=3,
                label=f"{gene} (AUC = {a:.2f})")
    ax.plot([0, 1], [0, 1], ls="--", color=GREY, lw=0.8, zorder=2)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.04)
    ax.set_xticks(np.arange(0, 1.01, 0.2))
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_xlabel("1 − specificity")
    ax.set_ylabel("Sensitivity")
    ax.legend(frameon=False, fontsize=6.5, loc="lower right")
    style_axis(ax)


def nomogram_axis(ax):
    coefs = pd.read_csv(os.path.join(PANEL, "nomogram_coefs.csv"))
    ranges = pd.read_csv(os.path.join(PANEL, "nomogram_ranges.csv"))
    c = dict(zip(coefs.term, coefs.coef))
    lo = dict(zip(ranges.term, ranges["min"]))
    hi = dict(zip(ranges.term, ranges["max"]))
    ax.set_xlim(-6, 106)
    ax.set_ylim(-0.85, 6.65)
    ax.axis("off")
    lp_min = c["Intercept"] + c["NDNF"] * hi["NDNF"] + c["PCDHB7"] * lo["PCDHB7"]
    lp_max = c["Intercept"] + c["NDNF"] * lo["NDNF"] + c["PCDHB7"] * hi["PCDHB7"]
    span = lp_max - lp_min
    ndnf_pts = lambda v: 100.0 * c["NDNF"] * (v - lo["NDNF"]) / span
    pcdhb7_pts = lambda v: 100.0 * c["PCDHB7"] * (v - lo["PCDHB7"]) / span
    ndnf_off = -ndnf_pts(hi["NDNF"])

    def axis_line(y, label):
        ax.plot([0, 100], [y, y], color="black", lw=0.9)
        ax.text(0, y + 0.24, label, ha="left", va="bottom", fontsize=7.5)

    def ticks(values, y, transform, fmt=lambda v: f"{v:g}", every=1):
        for i, v in enumerate(values):
            x = transform(v)
            ax.plot([x, x], [y, y + 0.10], color="black", lw=0.7)
            if i % every == 0:
                ax.text(x, y - 0.34, fmt(v), ha="center", va="top",
                        fontsize=7)

    axis_line(5.9, "Points")
    ticks(np.arange(0, 101, 20), 5.9, lambda v: v)
    axis_line(4.2, "NDNF")
    ticks(np.round(np.arange(-0.5, 3.01, 0.5), 1), 4.2,
          lambda v: ndnf_off + ndnf_pts(v), every=2)
    axis_line(2.6, "PCDHB7")
    ticks(np.round(np.arange(0.25, 3.51, 0.5), 2), 2.6,
          lambda v: pcdhb7_pts(v), every=2)
    axis_line(1.0, "Total points")
    ticks(np.arange(0, 101, 20), 1.0, lambda v: v)
    axis_line(-0.2, "IgAN probability")
    for p in (0.1, 0.3, 0.5, 0.7, 0.9):
        lp = np.log(p / (1 - p))
        x = (lp - lp_min) / (span / 100.0)
        ax.plot([x, x], [-0.2, -0.10], color="black", lw=0.7)
        ax.text(x, -0.56, f"{p:g}", ha="center", va="top", fontsize=6)
    ax.text(104, 6.45, "C-index = 0.848", ha="right", va="top", fontsize=7)


def calibration_panel(ax):
    cal = pd.read_csv(os.path.join(PANEL, "calibration.csv"))
    bins = np.quantile(cal.pred, np.linspace(0, 1, 11))
    centers, obs = [], []
    for i in range(10):
        m = (cal.pred >= bins[i]) & (cal.pred <= bins[i + 1])
        if m.sum() == 0:
            continue
        centers.append(cal.loc[m, "pred"].mean())
        obs.append(cal.loc[m, "obs"].mean())
    ax.plot([0, 1], [0, 1], ls="--", color=GREY, lw=0.8)
    ax.plot(centers, obs, "-o", color=BLUE, ms=3, lw=1.1)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Predicted probability")
    ax.set_ylabel("Observed proportion")
    ax.text(0.02, 0.96, "calibration slope = 1.00",
            transform=ax.transAxes, fontsize=7, va="top")
    style_axis(ax)


def dca_panel(ax):
    d = pd.read_csv(os.path.join(BENCH, "dca.csv"))
    ax.plot(d.threshold, d.net_benefit_model, color=BLUE, lw=1.1,
            label="Nomogram")
    ax.plot(d.threshold, d.net_benefit_all, color=GREY, lw=0.8,
            label="Treat all")
    ax.plot(d.threshold, d.net_benefit_none, color="#999999", lw=0.8,
            ls=":", label="Treat none")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Threshold probability")
    ax.set_ylabel("Net benefit")
    ax.legend(frameon=False, fontsize=6.5, loc="upper right")
    style_axis(ax)


def make_fig2():
    fig = plt.figure(figsize=(FIG_W_IN, 5.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 0.8],
                          left=0.075, right=0.985, top=0.90, bottom=0.09,
                          hspace=0.72, wspace=0.48)
    axa = fig.add_subplot(gs[0, 0])
    roc_panel(axa)
    panel_label(axa, "a")
    axb = fig.add_subplot(gs[:, 1])
    nomogram_axis(axb)
    panel_label(axb, "b")
    axc = fig.add_subplot(gs[1, 0])
    calibration_panel(axc)
    panel_label(axc, "c")
    axd = fig.add_subplot(gs[1, 2])
    dca_panel(axd)
    panel_label(axd, "d")
    return fig, []


# ---------------- Fig3 ----------------
def subtype_box(ax, df, col, plabel, ylab):
    c1 = df.loc[df.subtype == "C1", col]
    c2 = df.loc[df.subtype == "C2", col]
    box_jitter(ax, c1.to_numpy(), 1, BLUE)
    box_jitter(ax, c2.to_numpy(), 2, ORANGE)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["C1 (n = 12)", "C2 (n = 43)"])
    ax.set_ylabel(ylab)
    ymax = max(c1.max(), c2.max())
    ymin = min(c1.min(), c2.min())
    ax.set_ylim(ymin - 0.12, ymax + 0.32)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, prune="both"))
    ax.plot([1, 1, 2, 2], [ymax + 0.04, ymax + 0.13, ymax + 0.13,
                           ymax + 0.04], color="black", lw=0.7)
    ax.text(1.5, ymax + 0.17, plabel, ha="center", va="bottom", fontsize=7)
    style_axis(ax)


def make_fig3():
    df, _, _ = main_df()
    sub = df[df.subtype.notna()].copy()
    hm = pd.read_csv(os.path.join(BENCH, "subtype_heatmap.csv")).iloc[:24]
    anno = pd.read_csv(os.path.join(BENCH, "subtype_heatmap_anno.csv"))
    pca = pd.read_csv(os.path.join(BENCH, "subtype_pca.csv"))
    fig = plt.figure(figsize=(FIG_W_IN, 5.9))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.05, 1.0, 1.05],
                          height_ratios=[1.0, 1.25],
                          left=0.095, right=0.975, top=0.90, bottom=0.09,
                          hspace=0.55, wspace=0.42)
    axa = fig.add_subplot(gs[0, 0])
    subtype_box(axa, sub, "MRG_up", "P = 0.0105", "Mechanical score")
    panel_label(axa, "a")
    axb = fig.add_subplot(gs[0, 1])
    subtype_box(axb, sub, "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                "P = 0.438", "EMT pathway score")
    panel_label(axb, "b")
    axc = fig.add_subplot(gs[0, 2])
    for ct, c in [("C1", BLUE), ("C2", ORANGE)]:
        d = pca[pca.subtype == ct]
        axc.scatter(d.PC1, d.PC2, s=7, color=c, lw=0, label=f"{ct} (n = {len(d)})")
    axc.set_xlabel("PC1 (61.5%)")
    axc.set_ylabel("PC2 (4.9%)")
    axc.set_ylim(-3.0, 6.0)
    axc.legend(frameon=False, fontsize=6.5)
    style_axis(axc)
    axc.xaxis.set_major_locator(MaxNLocator(nbins=4, prune="both"))
    axc.yaxis.set_major_locator(MaxNLocator(nbins=4, prune="both"))
    # 标注离群样本（PC2 远高于其余点）
    axc.annotate("off-scale\n(PC2 = 14.5)", xy=(10.44, 14.51),
                 xytext=(11.5, 5.2), fontsize=5.5, ha="left", va="center",
                 color="#333333")
    axc.set_ylim(-3.0, 6.0)
    panel_label(axc, "c")
    axd = fig.add_subplot(gs[1, :])
    Z = hm.set_index("gene").to_numpy()
    cmap = LinearSegmentedColormap.from_list("cbo", [BLUE, "#FFFFFF", ORANGE])
    vmax = float(np.nanmax(np.abs(Z)))
    im = axd.imshow(Z, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
    axd.set_yticks(np.arange(len(hm)))
    axd.set_yticklabels(hm.gene, fontsize=5.5)
    axd.set_xticks([])
    axd.set_xlabel("IgAN samples ordered by subtype (C1 | C2)")
    axd.tick_params(length=0)
    for sp in axd.spines.values():
        sp.set_visible(False)
    # 亚型颜色条
    colmap = {"C1": BLUE, "C2": ORANGE}
    bar = np.array([colmap[s] for s in anno.subtype])
    from matplotlib.colors import ListedColormap
    axbar = fig.add_axes([0.11, 0.035, 0.865, 0.012])
    axbar.imshow([np.arange(len(bar))], cmap=ListedColormap([BLUE, ORANGE]),
                 aspect="auto")
    axbar.set_xticks([])
    axbar.set_yticks([])
    axbar.text(0.0, 2.2, "C1", transform=axbar.transAxes, fontsize=6,
               ha="center")
    axbar.text(0.985, 2.2, "C2", transform=axbar.transAxes, fontsize=6,
               ha="right")
    cb = fig.colorbar(im, ax=axd, orientation="horizontal", fraction=0.04,
                      pad=0.06)
    cb.ax.tick_params(labelsize=6.5, width=0.6, length=2)
    cb.set_label("Row z-score", fontsize=7)
    panel_label(axd, "d")
    return fig, []


# ---------------- Fig4 ----------------
def make_fig4():
    imm = pd.read_csv(os.path.join(PANEL, "immune28_30ctrl_results.csv"))
    cor = pd.read_csv(os.path.join(PROC, "sig_immune_cell_cor.csv"))
    sig = pd.read_csv(os.path.join(BENCH, "sig_expr.csv"))
    bind = pd.read_csv(os.path.join(PROC, "immune_bindea_scores.csv"))
    pmat = np.zeros((3, 18))
    for gi, gene in enumerate(["NDNF", "PCDHB7", "RRAGB"]):
        gexpr = sig[sig.gene == gene].iloc[0]
        for ci, cell in enumerate(cor.columns):
            x, y = [], []
            for _, row in bind.iterrows():
                if row["sample"] in gexpr.index:
                    x.append(float(gexpr[row["sample"]]))
                    y.append(float(row[cell]))
            _, p = sps.spearmanr(x, y)
            pmat[gi, ci] = p
    fig = plt.figure(figsize=(FIG_W_IN, 4.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05],
                          left=0.105, right=0.965, top=0.90, bottom=0.16,
                          wspace=0.42)
    short = {
        "Activated B cell": "Act. B cell", "Activated CD4 T cell": "Act. CD4 T",
        "Activated CD8 T cell": "Act. CD8 T",
        "Activated dendritic cell": "Act. DC",
        "CD56bright natural killer cell": "CD56bright NK",
        "CD56dim natural killer cell": "CD56dim NK",
        "Central memory CD4 T cell": "Tcm CD4",
        "Central memory CD8 T cell": "Tcm CD8",
        "Effector memeory CD4 T cell": "Tem CD4",
        "Effector memeory CD8 T cell": "Tem CD8",
        "Gamma delta T cell": "γδ T cell", "Immature  B cell": "Imm. B cell",
        "Immature dendritic cell": "Imm. DC", "Memory B cell": "Mem. B cell",
        "Natural killer T cell": "NKT cell", "Natural killer cell": "NK cell",
        "Plasmacytoid dendritic cell": "pDC", "Regulatory T cell": "Treg",
        "T follicular helper cell": "Tfh", "Type 1 T helper cell": "Th1",
        "Type 17 T helper cell": "Th17", "Type 2 T helper cell": "Th2",
        "Macrophage": "Macro", "Mast cell": "Mast",
    }
    imm = imm.sort_values("padj", ascending=False).reset_index(drop=True)
    labels = [short.get(c, c) for c in imm.celltype]
    hm = imm[["IgAN", "control"]].to_numpy()
    vmax = float(np.nanmax(np.abs(hm)))
    cmap = LinearSegmentedColormap.from_list("cbo", [BLUE, "#FFFFFF", ORANGE])
    axa = fig.add_subplot(gs[0, 0])
    im = axa.imshow(hm, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
    axa.set_xticks([0, 1])
    axa.set_xticklabels(["IgAN (n = 55)", "Control (n = 30)"])
    axa.set_xlim(-0.5, 1.8)
    axa.set_yticks(np.arange(len(imm)))
    axa.set_yticklabels(labels, fontsize=6.5)
    for r, row in imm.iterrows():
        if row.padj < 0.05:
            axa.text(1.65, r, "*", ha="center", va="center", fontsize=8,
                     color="black")
    axa.tick_params(length=0)
    for sp in axa.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=axa, orientation="horizontal", fraction=0.045,
                      pad=0.08)
    cb.ax.tick_params(labelsize=6.5, width=0.6, length=2)
    cb.set_ticks([-vmax, 0.0, vmax])
    cb.set_ticklabels([f"{x:.2f}" for x in (-vmax, 0.0, vmax)])
    cb.set_label("ssGSEA score", fontsize=7)
    panel_label(axa, "a")

    axb = fig.add_subplot(gs[0, 1])
    C = cor.to_numpy()
    vmax2 = float(np.nanmax(np.abs(C)))
    im2 = axb.imshow(C, cmap=cmap, vmin=-vmax2, vmax=vmax2, aspect="auto")
    axb.set_xticks(np.arange(len(cor.columns)))
    short24 = {"B_cells": "B", "T_cells": "T", "T_helper_cells": "Th",
               "Tcm": "Tcm", "Tem": "Tem", "Th1_cells": "Th1",
               "Th2_cells": "Th2", "TFH": "Tfh", "CD8_T_cells": "CD8",
               "Cytotoxic_cells": "Cyto", "NK_cells": "NK",
               "NK_CD56dim_cells": "NK56d", "DC": "DC", "iDC": "iDC",
               "Eosinophils": "Eos", "Macrophages": "Macro",
               "Mast_cells": "Mast", "Neutrophils": "Neu"}
    axb.set_xticks(np.arange(len(cor.columns)))
    axb.set_xticklabels([short24[c] for c in cor.columns],
                        fontsize=5.5, rotation=90)
    axb.set_yticks([0, 1, 2])
    axb.set_yticklabels(["NDNF", "PCDHB7", "RRAGB"], fontstyle="italic")
    for i in range(C.shape[0]):
        for j in range(C.shape[1]):
            axb.text(j, i, f"{C[i, j]:.2f}", ha="center", va="center",
                     fontsize=4.5, color="#222222")
            if pmat[i, j] < 0.05:
                axb.text(j, i - 0.36, "*", ha="center", va="center",
                         fontsize=6, color="black")
    axb.tick_params(length=0)
    for sp in axb.spines.values():
        sp.set_visible(False)
    cb2 = fig.colorbar(im2, ax=axb, orientation="horizontal", fraction=0.045,
                       pad=0.08)
    cb2.ax.tick_params(labelsize=6.5, width=0.6, length=2)
    cb2.set_ticks([-vmax2, 0.0, vmax2])
    cb2.set_ticklabels([f"{x:.2f}" for x in (-vmax2, 0.0, vmax2)])
    cb2.set_label("Spearman ρ", fontsize=7)
    panel_label(axb, "b")
    fig.text(0.02, 0.012,
             "* BH-adjusted P < 0.05 (a: two-sided Wilcoxon rank-sum test, "
             "IgAN vs. control; b: two-sided Spearman correlation).",
             fontsize=6.5, va="bottom", ha="left")
    return fig, []


# ---------------- Fig6 ----------------
def make_fig6():
    go = pd.read_csv(os.path.join(BENCH, "go_bp_mrg_deg.csv"))
    kegg = pd.read_csv(os.path.join(BENCH, "kegg_mrg_deg.csv"))
    edges = pd.read_csv(os.path.join(BENCH, "string_edges.csv"))
    nodes = pd.read_csv(os.path.join(BENCH, "string_nodes.csv"))
    fig = plt.figure(figsize=(FIG_W_IN, 3.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 0.85],
                          left=0.24, right=0.96, top=0.90, bottom=0.13,
                          wspace=0.38)
    axa = fig.add_subplot(gs[0, 0])
    top = go.sort_values("pvalue").head(12).iloc[::-1]
    sc = axa.scatter(top.Count, top.Description,
                     c=-np.log10(top["p.adjust"]), cmap="viridis",
                     s=np.clip(top.Count * 12, 8, 120), edgecolors="#333333",
                     lw=0.4, zorder=3)
    axa.set_xlabel("Gene count")
    axa.set_xlim(0, max(top.Count) + 3)
    axa.set_yticks(np.arange(len(top)))
    desc = [d if len(d) <= 42 else d[:38] + "..." for d in top.Description]
    axa.set_yticklabels(desc, fontsize=5.5)
    cb = fig.colorbar(sc, ax=axa, orientation="horizontal", fraction=0.04,
                      pad=0.12)
    cb.set_label("−log10(adjusted P)", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    style_axis(axa)
    panel_label(axa, "a", "GO-BP (36 MRG ∩ DEG genes)")

    axb = fig.add_subplot(gs[0, 1])
    axb.axis("off")
    axb.set_xlim(-0.15, 1.15)
    axb.set_ylim(-0.15, 1.15)
    axb.set_xticks([])
    axb.set_yticks([])
    if len(edges) and "preferredName_A" in edges.columns:
        import networkx as nx
        G = nx.Graph()
        for _, r in edges.iterrows():
            G.add_edge(r.preferredName_A, r.preferredName_B)
        pos = nx.spring_layout(G, seed=42)
        deg = dict(G.degree())
        for u, v in G.edges():
            x1, y1 = pos[u]
            x2, y2 = pos[v]
            axb.plot([x1, x2], [y1, y2], color="#BBBBBB", lw=0.6, zorder=2)
        for n, (x, y) in pos.items():
            d = deg[n]
            axb.scatter(x, y, s=60 + 30 * d, color=BLUE if d >= 2 else "#88CCEE",
                        edgecolors="#333333", lw=0.4, zorder=3)
            if d >= 2:
                axb.text(x, y + 0.06, n, ha="center", va="bottom", fontsize=6)
        axb.text(0.5, 1.10, "STRING PPI of the 100 core genes\n"
                 f"({len(G.nodes())} nodes, {len(G.edges())} edges, "
                 "score ≥ 0.4)", ha="center", va="top", fontsize=6.5)
        axb.text(0.5, -0.06, "Node size scales with degree; "
                 "labeled hubs (degree ≥ 2)",
                 ha="center", va="bottom", fontsize=6, color="#555555")
    else:
        axb.text(5, 5, "PPI data not available", ha="center", fontsize=7)
    panel_label(axb, "b", "Protein–protein interaction")

    # KEGG 注释
    fig.text(0.02, 0.012,
             "KEGG enrichment of the 36 MRG ∩ DEG genes: Hippo signaling "
             "pathway (5 genes; adjusted P = 0.0058).",
             fontsize=6.5, va="bottom", ha="left")
    return fig, []


# ---------------- Fig5 复用（= 旧 Fig4） ----------------
def load_fig31():
    spec = importlib.util.spec_from_file_location(
        "fig31v2", os.path.join(ROOT, "scripts", "31_figures_plos.py"))
    m = importlib.util.module_from_spec(spec)
    sys.modules["fig31v2"] = m
    spec.loader.exec_module(m)
    return m


def fix_fig4_ticks(fig):
    axes = []
    for ax in fig.axes:
        axes.append(ax)
        axes.extend(getattr(ax, "child_axes", []) or [])
    for ax in axes:
        if ax.get_xlabel() == "Mechanical score (FOV)":
            x0, x1 = ax.get_xlim()
            y0, y1 = ax.get_ylim()
            ax.set_xticks(np.round(np.linspace(x0, x1, 3), 2))
            ax.set_yticks(np.round(np.linspace(y0, y1, 3), 2))


def annotate_fig5_stats(fig):
    """在 Fig5 的 sham vs UUO 细胞类型面板上补 Wilcoxon BH 显著性。"""
    uuo = pd.read_csv(os.path.join(PROC, "uuo_meta.csv"))
    order = ["PT", "TAL", "DCT", "CD", "Podo", "Endo", "Fibro", "Macro",
             "Tcell", "Bcell"]
    pvals = []
    for ct in order:
        a = uuo.loc[(uuo.celltype == ct) & (uuo.group == "sham"), "MRG1"]
        b = uuo.loc[(uuo.celltype == ct) & (uuo.group == "UUO"), "MRG1"]
        _, p = sps.mannwhitneyu(a, b, alternative="two-sided")
        pvals.append(float(p))
    m = len(pvals)
    order_idx = np.argsort(pvals)
    q = np.empty(m)
    prev = np.inf
    for i in range(m - 1, -1, -1):
        r = order_idx[i]
        cur = min(pvals[r] * m / (i + 1), prev)
        q[r] = cur
        prev = cur
    for ax in fig.axes:
        if ax.get_ylabel() != "Mechanical core module score":
            continue
        for i, ct in enumerate(order):
            if q[i] < 0.05:
                ax.text(i + 0.17, 0.115, "*", ha="center", va="center",
                        fontsize=8, color="black")
        break


def layout_report(fig, name, pairs):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.get_size_inches() * DPI
    issues = []
    axes = []
    for ax in fig.axes:
        axes.append(ax)
        axes.extend(getattr(ax, "child_axes", []) or [])
    for ax in axes:
        texts = list(ax.texts) + list(ax.get_xticklabels()) \
            + list(ax.get_yticklabels())
        for t in texts:
            if not t.get_visible():
                continue
            bb = t.get_window_extent(renderer=r)
            if bb.x0 < -1 or bb.y0 < -1 or bb.x1 > W + 1 or bb.y1 > H + 1:
                issues.append(f"outside: '{t.get_text()[:24]}'")
        for labels, axis in ((ax.get_yticklabels(), "y"),
                             (ax.get_xticklabels(), "x")):
            vis = [t for t in labels if t.get_visible() and t.get_text().strip()]
            vis.sort(key=lambda t: t.get_window_extent(renderer=r).y0)
            for a, b in zip(vis, vis[1:]):
                ba = a.get_window_extent(renderer=r)
                bb = b.get_window_extent(renderer=r)
                if (ba.x1 > bb.x0 + 1 and bb.x1 > ba.x0 + 1
                        and ba.y1 > bb.y0 + 1 and bb.y1 > ba.y0 + 1):
                    issues.append(f"tick-overlap({axis}): "
                                  f"'{a.get_text()[:12]}' vs '{b.get_text()[:12]}'")
        anns = [t for t in ax.texts if t.get_visible() and t.get_text().strip()]
        for i in range(len(anns)):
            for j in range(i + 1, len(anns)):
                ba = anns[i].get_window_extent(renderer=r)
                bb = anns[j].get_window_extent(renderer=r)
                if (ba.x1 > bb.x0 + 1 and bb.x1 > ba.x0 + 1
                        and ba.y1 > bb.y0 + 1 and bb.y1 > ba.y0 + 1):
                    issues.append(f"overlap: '{anns[i].get_text()[:16]}' vs "
                                  f"'{anns[j].get_text()[:16]}'")
    if issues:
        log(f"{name}: {len(issues)} issue(s)")
        for s in issues[:12]:
            log("   " + s)
    else:
        log(f"{name}: layout OK")
    return issues


def save_fig(fig, name, force):
    png_path = os.path.join(PREVDIR, f"{name}.png")
    pdf_path = os.path.join(OUTDIR, f"{name}.pdf")
    tif_path = os.path.join(OUTDIR, f"{name}.tif")
    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(PREVDIR, exist_ok=True)
    for p in (png_path, pdf_path, tif_path):
        if os.path.exists(p) and not force:
            raise SystemExit(f"refusing to overwrite: {p}")
    fig.savefig(pdf_path, format="pdf", dpi=DPI, facecolor="white")
    fig.savefig(png_path, format="png", dpi=DPI, facecolor="white")
    with Image.open(png_path) as im:
        rgb = np.asarray(im.convert("RGB"))
    tifffile.imwrite(tif_path, rgb, photometric="rgb", planarconfig="contig",
                     resolution=(DPI, DPI), resolutionunit="inch",
                     compression="lzw")
    plt.close(fig)
    with Image.open(png_path) as im:
        log(f"{name}: {im.size[0]}x{im.size[1]} px")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    setup_fonts()

    builders = [("Fig1", make_fig1), ("Fig2", make_fig2),
                ("Fig3", make_fig3), ("Fig4", make_fig4),
                ("Fig6", make_fig6)]
    for name, fn in builders:
        patch_left(0.02 if name == "Fig2" else 0.0)
        fig, pairs = fn()
        issues = layout_report(fig, name, pairs)
        save_fig(fig, name, args.force)

    # Fig5 = 旧 Fig4（scRNA/空间/时间）
    m = load_fig31()
    m.FIG_W_IN = 7.0
    m.panel_label = lambda ax, t: ax.text(0.0, 1.02, t,
                                          transform=ax.transAxes,
                                          fontsize=9, fontweight="bold",
                                          va="bottom", ha="left", color="black")
    patch_left(0.02)
    fig5, _ = m.make_fig4()
    fix_fig4_ticks(fig5)
    annotate_fig5_stats(fig5)
    layout_report(fig5, "Fig5", [])
    save_fig(fig5, "Fig5", args.force)
    log("DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

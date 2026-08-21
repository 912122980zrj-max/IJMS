#!/usr/bin/env python
# 58_rebuild_fig2.py —— 重建改进版 Fig2（a 签名ROC+逐基因、b 列线图、c 校准、d DCA）
# 逐基因 AUC/ROC 用 85 集口径（per_gene_auc.csv / per_gene_roc.csv 已更新）。
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from matplotlib import font_manager
from PIL import Image
from sklearn.metrics import roc_auc_score, roc_curve

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"E:/sheng xin/ObstructiveNephropathy_MRG"
PANEL = os.path.join(ROOT, "data", "processed", "panel_export")
BENCH = os.path.join(ROOT, "results", "benchmark")
OUTDIR = os.path.join(ROOT, "submission", "ijms", "figures")
PREVDIR = os.path.join(OUTDIR, "preview")
DPI = 300

CB = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"]
BLUE = CB[0]
ORANGE = CB[1]
GREY = "#666666"
GENE_COL = {"NDNF": CB[3], "PCDHB7": CB[1], "RRAGB": CB[2]}


def setup_fonts():
    for fam in ("Arial", "Helvetica", "DejaVu Sans"):
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


def style_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=0.6, length=2.5)


def panel_label(ax, text):
    ax.text(0.0, 1.02, text, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="bottom", ha="left", color="black")


def signature_roc_data(seed=7, n_boot=2000):
    s = pd.read_csv(os.path.join(PANEL, "signature_scores.csv"))
    y = s["y"].to_numpy()
    sc = s["score"].to_numpy()
    n = len(y)
    n_pos = int(y.sum())
    n_neg = int(n - n_pos)
    grid = np.linspace(0, 1, 201)
    fpr, tpr, _ = roc_curve(y, sc)
    tpr = np.interp(grid, fpr, tpr)
    tpr[0] = 0.0
    auc = roc_auc_score(y, sc)
    rng = np.random.default_rng(seed)
    boot_t, aucs = [], []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y[idx], sc[idx]))
        fb, tb, _ = roc_curve(y[idx], sc[idx])
        boot_t.append(np.interp(grid, fb, tb))
    boot_t = np.array(boot_t)
    return dict(grid=grid, tpr=tpr,
                lo=np.percentile(boot_t, 2.5, axis=0),
                hi=np.percentile(boot_t, 97.5, axis=0),
                auc=auc,
                auc_lo=np.percentile(aucs, 2.5), auc_hi=np.percentile(aucs, 97.5),
                n=n, n_pos=n_pos, n_neg=n_neg)


def roc_panel(ax):
    d = signature_roc_data()
    per_auc = pd.read_csv(os.path.join(BENCH, "per_gene_auc.csv"))
    per_roc = pd.read_csv(os.path.join(BENCH, "per_gene_roc.csv"))
    ax.fill_between(d["grid"], d["lo"], d["hi"], color=BLUE, alpha=0.15, lw=0, zorder=1)
    ax.plot(d["grid"], d["tpr"], color=BLUE, lw=1.6, zorder=4)
    for gene, c in GENE_COL.items():
        sub = per_roc[per_roc.gene == gene]
        aa = per_auc.loc[per_auc.gene == gene, "AUC"].iloc[0]
        ax.plot(sub.fpr, sub.tpr, color=c, lw=0.9, zorder=3)
    ax.plot([0, 1], [0, 1], ls="--", color=GREY, lw=0.8, zorder=2)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.03, 1.04)
    ax.set_xticks(np.arange(0, 1.01, 0.2))
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_xlabel("1 − specificity")
    ax.set_ylabel("Sensitivity")
    ax.text(0.60, 0.35,
            f"AUC = 0.853  (95% CI {d['auc_lo']:.3f}–{d['auc_hi']:.3f})\n"
            f"5-fold CV = 0.813 ± 0.079\n"
            f"n = {d['n']}  ({d['n_pos']} IgAN / {d['n_neg']} control)",
            fontsize=6.5, ha="left", va="center", linespacing=1.45)
    # legend: signature + 3 genes
    handles = [plt.Line2D([0], [0], color=BLUE, lw=1.6,
                          label="3-gene signature (AUC = 0.853)")]
    for gene in ["NDNF", "PCDHB7", "RRAGB"]:
        aa = per_auc.loc[per_auc.gene == gene, "AUC"].iloc[0]
        handles.append(plt.Line2D([0], [0], color=GENE_COL[gene], lw=0.9,
                                  label=f"{gene} (AUC = {aa:.2f})"))
    ax.legend(handles=handles, frameon=False, fontsize=6, loc="upper left",
              bbox_to_anchor=(0.02, 0.99))
    style_axis(ax)
    panel_label(ax, "a")


def load_nomogram():
    coefs = pd.read_csv(os.path.join(PANEL, "nomogram_coefs.csv"))
    ranges = pd.read_csv(os.path.join(PANEL, "nomogram_ranges.csv"))
    c = dict(zip(coefs.term, coefs.coef))
    lo = dict(zip(ranges.term, ranges["min"]))
    hi = dict(zip(ranges.term, ranges["max"]))
    lp_min = c["Intercept"] + c["NDNF"] * hi["NDNF"] + c["PCDHB7"] * lo["PCDHB7"]
    lp_max = c["Intercept"] + c["NDNF"] * lo["NDNF"] + c["PCDHB7"] * hi["PCDHB7"]
    return c, lo, hi, lp_min, lp_max


def nomogram_axis(ax, dense=False):
    c, lo, hi, lp_min, lp_max = load_nomogram()
    span = lp_max - lp_min
    ndnf_pts = lambda v: 100.0 * c["NDNF"] * (v - lo["NDNF"]) / span
    pcdhb7_pts = lambda v: 100.0 * c["PCDHB7"] * (v - lo["PCDHB7"]) / span
    ndnf_off = -ndnf_pts(hi["NDNF"])
    fs = 6.5
    ax.set_xlim(-6, 106)
    ax.set_ylim(-1.05, 7.05)
    ax.axis("off")
    ax.set_xticks([])
    ax.set_yticks([])

    def axis_line(y, label, style="normal"):
        ax.plot([0, 100], [y, y], color="black", lw=0.9)
        ax.text(0, y + 0.22, label, ha="left", va="bottom", fontsize=7.5, style=style)

    def ticks(values, y, transform, fmt=lambda v: f"{v:g}", every=1):
        for i, v in enumerate(values):
            x = transform(v)
            ax.plot([x, x], [y, y + 0.09], color="black", lw=0.7)
            if i % every == 0:
                ax.text(x, y - 0.30, fmt(v), ha="center", va="top", fontsize=fs)

    axis_line(6.5, "Points")
    ticks(np.arange(0, 101, 20), 6.5, lambda v: v)
    axis_line(5.0, "NDNF", "italic")
    ticks(np.round(np.arange(lo["NDNF"], hi["NDNF"] + 0.001, 1.0), 2), 5.0,
          lambda v: ndnf_off + ndnf_pts(v), every=1)
    axis_line(3.5, "PCDHB7", "italic")
    ticks(np.round(np.arange(0.25, hi["PCDHB7"] + 0.001, 0.5), 2), 3.5,
          lambda v: pcdhb7_pts(v), every=2)
    axis_line(2.0, "Total points")
    ticks(np.arange(0, 101, 20), 2.0, lambda v: v)
    axis_line(0.5, "Linear predictor")
    for lp in (t for t in [-3, -1, 1, 3, 5, 7] if lp_min - 0.2 <= t <= lp_max + 0.2):
        x = (lp - lp_min) / (span / 100.0)
        ax.plot([x, x], [0.5, 0.59], color="black", lw=0.7)
        ax.text(x, 0.5 - 0.30, f"{lp:g}", ha="center", va="top", fontsize=fs)
    axis_line(-0.85, "IgAN probability")
    for p in (0.1, 0.3, 0.5, 0.7, 0.9):
        lp = math.log(p / (1 - p))
        x = (lp - lp_min) / (span / 100.0)
        ax.plot([x, x], [-0.85, -0.76], color="black", lw=0.7)
        ax.text(x, -1.22, f"{p:g}", ha="center", va="top", fontsize=fs)


def calibration_panel(ax):
    cal = pd.read_csv(os.path.join(PANEL, "calibration.csv"))
    y = cal["obs"].to_numpy()
    p = cal["pred"].to_numpy()
    lp = np.log(p / (1 - p))
    X = np.column_stack([np.ones(len(lp)), lp])
    beta = np.zeros(2)
    for _ in range(60):
        eta = X @ beta
        mu = 1 / (1 + np.exp(-eta))
        W = mu * (1 - mu)
        z = eta + (y - mu) / np.maximum(W, 1e-9)
        beta = np.linalg.solve(X.T @ (W[:, None] * X), X.T @ (W * z))
    intercept, slope = beta
    brier = float(np.mean((p - y) ** 2))
    edges = np.quantile(p, np.linspace(0, 1, 11))
    centers, obs = [], []
    for i in range(10):
        m = (p >= edges[i]) & (p <= edges[i + 1])
        if m.sum() == 0:
            continue
        centers.append(p[m].mean())
        obs.append(y[m].mean())
    ax.plot([0, 1], [0, 1], ls="--", color=GREY, lw=0.8)
    ax.plot(centers, obs, "-o", color=BLUE, ms=3, lw=1.0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xticks(np.arange(0, 1.01, 0.2))
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_xlabel("Predicted probability")
    ax.set_ylabel("Observed proportion")
    ax.text(0.02, 0.96, f"slope = {slope:.2f}   intercept = {intercept:.2f}\nBrier = {brier:.3f}",
            transform=ax.transAxes, fontsize=6.5, va="top", linespacing=1.35)
    style_axis(ax)
    panel_label(ax, "c")


def dca_panel(ax):
    d = pd.read_csv(os.path.join(BENCH, "dca.csv"))
    ax.plot(d.threshold, d.net_benefit_model, color=BLUE, lw=1.1, label="Nomogram")
    ax.plot(d.threshold, d.net_benefit_all, color=GREY, lw=0.8, label="Treat all")
    ax.plot(d.threshold, d.net_benefit_none, color="#999999", lw=0.8, ls=":", label="Treat none")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.05, 0.75)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticks([0, 0.2, 0.4, 0.6])
    ax.set_xlabel("Threshold probability")
    ax.set_ylabel("Net benefit")
    ax.legend(frameon=False, fontsize=6, loc="upper right")
    style_axis(ax)
    panel_label(ax, "d")


def make_fig2():
    fig = plt.figure(figsize=(7.0, 5.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.05, 1.45, 1.0],
                          left=0.075, right=0.985, top=0.90, bottom=0.09,
                          hspace=0.50, wspace=0.42)
    roc_panel(fig.add_subplot(gs[0, 0]))
    nomogram_axis(fig.add_subplot(gs[:, 1]))
    panel_label(fig.axes[1], "b")
    calibration_panel(fig.add_subplot(gs[1, 0]))
    dca_panel(fig.add_subplot(gs[1, 2]))
    return fig


def save_fig(fig, name):
    png_path = os.path.join(PREVDIR, f"{name}.png")
    pdf_path = os.path.join(OUTDIR, f"{name}.pdf")
    tif_path = os.path.join(OUTDIR, f"{name}.tif")
    fig.savefig(pdf_path, format="pdf", dpi=DPI, facecolor="white")
    fig.savefig(png_path, format="png", dpi=DPI, facecolor="white")
    with Image.open(png_path) as im:
        rgb = np.asarray(im.convert("RGB"))
    tifffile.imwrite(tif_path, rgb, photometric="rgb", planarconfig="contig",
                     resolution=(DPI, DPI), resolutionunit="inch", compression="lzw")
    plt.close(fig)
    with Image.open(png_path) as im:
        print(f"{name}: {im.size[0]}x{im.size[1]} px -> {tif_path}")


def main():
    setup_fonts()
    save_fig(make_fig2(), "Fig2")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

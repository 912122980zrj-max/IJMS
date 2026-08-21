#!/usr/bin/env python
"""31_figures_plos.py

按 PLOS Computational Biology 图件规范重绘并合成 4 幅多面板图（2026-08-20 重排版）：
  - 单文件多面板 TIFF（LZW、RGB、无 alpha）；宽度 2250 px（19.05 cm @300 dpi），高度 <= 2625 px；
  - 字体 Arial 8-12 pt；色盲友好配色（Okabe-Ito / Tol muted）；
  - 布局宽松：Fig1 上下两行（workflow 全宽 + 评分图）；Fig3 加高容纳 28 行热图标签；
    Fig4 三行（UMAP/空间 - 细胞类型全宽 - 时间序列全宽），图例独立占列；
  - 输出前用 renderer 实测每个文本边界：文本不得越出画布，Fig1 条目文本不得越出方框。

输出：submission/plos/figures/FigN.tif + preview/FigN.png
用法：
    python scripts/31_figures_plos.py [--force]
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.ticker import MaxNLocator
from PIL import Image

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


ROOT = r"E:/sheng xin/ObstructiveNephropathy_MRG"
PANEL = os.path.join(ROOT, "data", "processed", "panel_export")
PROC = os.path.join(ROOT, "data", "processed")
OUTDIR = os.path.join(ROOT, "submission", "plos", "figures")
PREVDIR = os.path.join(OUTDIR, "preview")

FIG_W_IN = 7.5
DPI = 300

CB = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"]
BLUE = CB[0]
ORANGE = CB[1]
GREY = "#666666"

TOL_MUTED = ["#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
             "#DDCC77", "#CC6677", "#882255", "#AA4499", "#6699CC"]


def log(msg: str) -> None:
    print(f"[figures] {time.strftime('%H:%M:%S')} {msg}", flush=True)


def setup_fonts() -> None:
    for fam in ("Arial",):
        if any(f.name == fam for f in font_manager.fontManager.ttflist):
            plt.rcParams["font.family"] = fam
            break
    plt.rcParams.update({
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 9,
        "axes.linewidth": 0.6,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "axes.edgecolor": "#333333",
        "xtick.color": "#333333",
        "ytick.color": "#333333",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def panel_label(ax, text: str) -> None:
    """面板标签置于面板上方（va=bottom），避免压住数据区。"""
    ax.text(0.0, 1.02, text, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="bottom", ha="left", color="black")


def style_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=0.6, length=2.5)


def read_main_df():
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


# ---------------------------------------------------------------------------
# Fig 1（上下两行：a workflow 全宽，b 机械评分）
# ---------------------------------------------------------------------------

def workflow_axis(ax):
    """绘制全宽三栏工作流；返回 (文本, 方框数据坐标矩形, ax) 三元组列表供自检。"""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    boxes = [
        (0.25, 1.2, 2.7, 8.6, "#EAF2F8", BLUE, "Public transcriptomes", [
            ["Human biopsies", "GSE115857 (n = 86)"],
            ["Mouse scRNA-seq", "GSE175412 (3 kidneys)"],
            ["Human spatial CosMx", "GSE282059 (523,855 cells)"],
            ["Mouse UUO time course", "GSE118339 (n = 15)"],
        ]),
        (3.65, 1.2, 2.7, 8.6, "#FFF6E6", ORANGE, "Mechanical framework", [
            ["2,536 mechanical genes", "(3 GO categories)"],
            ["ssGSEA mechanical score"],
            ["WGCNA β = 4: 38 modules", "100-gene core (ME31)"],
            ["LASSO ∩ random forest", "NDNF, PCDHB7, RRAGB"],
        ]),
        (7.05, 1.2, 2.7, 8.6, "#EDF7F1", CB[2], "Outcomes", [
            ["IgAN vs controls", "P = 7.01 × 10$^{-10}$"],
            ["ROC AUC 0.853", "C-index 0.848"],
            ["Subtypes C1/C2", "P = 0.0105"],
            ["CD hub: spatial ρ = 0.586", "time course ρ = 0.94"],
        ]),
    ]
    pairs = []
    for x, y, w, h, face, edge, title, items in boxes:
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                                    boxstyle="round,pad=0.05,rounding_size=0.15",
                                    facecolor=face, edgecolor=edge, lw=1.0))
        t = ax.text(x + w / 2, y + h - 0.42, title, ha="center", va="top",
                    fontsize=8.5, fontweight="bold", color=edge)
        pairs.append((t, (x, y, w, h), ax))
        ty = y + h - 1.02
        for item in items:
            for line in item:
                t = ax.text(x + w / 2, ty, line, ha="center", va="top",
                            fontsize=8.0, color="#222222", linespacing=1.12)
                pairs.append((t, (x, y, w, h), ax))
                ty -= 0.74
            ty -= 0.26
    for xa in (2.95, 6.35):
        ax.add_patch(FancyArrowPatch((xa, 5.0), (xa + 0.68, 5.0),
                                     arrowstyle="-|>", mutation_scale=13,
                                     color="#333333", lw=1.0))
    return pairs


def make_fig1():
    df, ctrl, iga = read_main_df()
    fig = plt.figure(figsize=(FIG_W_IN, 5.05))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.22, 0.9],
                          left=0.09, right=0.985, top=0.90, bottom=0.08,
                          hspace=0.42)
    axa = fig.add_subplot(gs[0, 0])
    pairs = workflow_axis(axa)
    panel_label(axa, "a")

    axb = fig.add_subplot(gs[1, 0])
    box_jitter(axb, ctrl.MRG_up.to_numpy(), 1, BLUE)
    box_jitter(axb, iga.MRG_up.to_numpy(), 2, ORANGE)
    axb.set_xticks([1, 2])
    axb.set_xticklabels(["Control (n = 30)", "IgAN (n = 55)"])
    axb.set_xlim(0.45, 2.55)
    axb.set_ylabel("Mechanical score (ssGSEA, MRG-up)")
    ymax = max(ctrl.MRG_up.max(), iga.MRG_up.max())
    axb.set_ylim(min(ctrl.MRG_up.min(), iga.MRG_up.min()) - 0.15, ymax + 0.38)
    axb.yaxis.set_major_locator(MaxNLocator(nbins=5, prune="both"))
    axb.plot([1, 1, 2, 2], [ymax + 0.09, ymax + 0.20, ymax + 0.20, ymax + 0.09],
             color="black", lw=0.7)
    axb.text(1.5, ymax + 0.25, "P = 7.01 × 10$^{-10}$", ha="center",
             va="bottom", fontsize=8)
    style_axis(axb)
    panel_label(axb, "b")
    return fig, pairs


# ---------------------------------------------------------------------------
# Fig 2（a ROC，b nomogram）
# ---------------------------------------------------------------------------

def make_fig2():
    roc = pd.read_csv(os.path.join(PANEL, "roc_curve.csv"))
    coefs = pd.read_csv(os.path.join(PANEL, "nomogram_coefs.csv"))
    ranges = pd.read_csv(os.path.join(PANEL, "nomogram_ranges.csv"))
    c = dict(zip(coefs.term, coefs.coef))
    lo = dict(zip(ranges.term, ranges["min"]))
    hi = dict(zip(ranges.term, ranges["max"]))

    fig = plt.figure(figsize=(FIG_W_IN, 3.7))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.9, 1.35],
                          left=0.055, right=0.985, top=0.90, bottom=0.12,
                          wspace=0.38)

    axa = fig.add_subplot(gs[0, 0])
    axa.plot(roc.spec, roc.sens, color=BLUE, lw=1.2, zorder=3)
    axa.plot([0, 1], [0, 1], ls="--", color=GREY, lw=0.8, zorder=2)
    axa.set_xlim(-0.02, 1.02)
    axa.set_ylim(-0.02, 1.04)
    axa.set_xticks(np.arange(0, 1.01, 0.2))
    axa.set_yticks(np.arange(0, 1.01, 0.2))
    axa.set_xlabel("1 − specificity")
    axa.set_ylabel("Sensitivity")
    axa.text(0.60, 0.10, "AUC = 0.853\n5-fold CV = 0.813 ± 0.079",
             fontsize=8, ha="left", va="bottom")
    style_axis(axa)
    panel_label(axa, "a")

    lp_min = c["Intercept"] + c["NDNF"] * hi["NDNF"] + c["PCDHB7"] * lo["PCDHB7"]
    lp_max = c["Intercept"] + c["NDNF"] * lo["NDNF"] + c["PCDHB7"] * hi["PCDHB7"]
    span = lp_max - lp_min

    def ndnf_pts(v):
        return 100.0 * c["NDNF"] * (v - lo["NDNF"]) / span

    def pcdhb7_pts(v):
        return 100.0 * c["PCDHB7"] * (v - lo["PCDHB7"]) / span

    ndnf_off = -ndnf_pts(hi["NDNF"])
    axb = fig.add_subplot(gs[0, 1])
    axb.set_xlim(-6, 106)
    axb.set_ylim(-0.85, 6.65)
    axb.axis("off")
    axb.set_xticks([])
    axb.set_yticks([])

    def draw_axis(y, label):
        axb.plot([0, 100], [y, y], color="black", lw=0.9)
        axb.text(0, y + 0.24, label, ha="left", va="bottom", fontsize=8.5)

    def ticks(values, y, transform, fmt=lambda v: f"{v:g}", label_every=1):
        for i, v in enumerate(values):
            x = transform(v)
            axb.plot([x, x], [y, y + 0.10], color="black", lw=0.7)
            if i % label_every == 0:
                axb.text(x, y - 0.34, fmt(v), ha="center", va="top",
                         fontsize=8)

    draw_axis(5.9, "Points")
    ticks(np.arange(0, 101, 20), 5.9, lambda v: v)

    draw_axis(4.2, "NDNF")
    # 刻度每 0.5 一处，标签每 1.0 一处，避免标签互相重叠
    ticks(np.round(np.arange(-0.5, 3.01, 0.5), 1), 4.2,
          lambda v: ndnf_off + ndnf_pts(v), label_every=2)

    draw_axis(2.6, "PCDHB7")
    ticks(np.round(np.arange(0.25, 3.51, 0.5), 2), 2.6,
          lambda v: pcdhb7_pts(v))

    draw_axis(1.0, "Total points")
    ticks(np.arange(0, 101, 20), 1.0, lambda v: v)

    draw_axis(-0.2, "IgAN probability")
    for p in (0.1, 0.3, 0.5, 0.7, 0.9):
        lp = math.log(p / (1 - p))
        x = (lp - lp_min) / (span / 100.0)
        axb.plot([x, x], [-0.2, -0.10], color="black", lw=0.7)
        axb.text(x, -0.56, f"{p:g}", ha="center", va="top", fontsize=8)
    axb.text(104, 6.45, "C-index = 0.848\ncalibration slope = 1.00",
             ha="right", va="top", fontsize=8)
    panel_label(axb, "b")
    return fig, []


# ---------------------------------------------------------------------------
# Fig 3（a/b 亚型箱线图堆叠在左，c 28 细胞热图占右列）
# ---------------------------------------------------------------------------

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
    ax.plot([1, 1, 2, 2], [ymax + 0.04, ymax + 0.13, ymax + 0.13, ymax + 0.04],
            color="black", lw=0.7)
    ax.text(1.5, ymax + 0.17, plabel, ha="center", va="bottom", fontsize=8)
    style_axis(ax)


def make_fig3():
    df, _, _ = read_main_df()
    sub = df[df.subtype.notna()].copy()
    imm = pd.read_csv(os.path.join(PANEL, "immune28_30ctrl_results.csv"))

    fig = plt.figure(figsize=(FIG_W_IN, 5.35))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.55, 1.0],
                          height_ratios=[1.0, 1.0],
                          left=0.105, right=0.955, top=0.90, bottom=0.09,
                          hspace=0.62, wspace=0.40)

    axa = fig.add_subplot(gs[0, 0])
    subtype_box(axa, sub, "MRG_up", "P = 0.0105", "Mechanical score (ssGSEA)")
    panel_label(axa, "a")

    axb = fig.add_subplot(gs[1, 0])
    subtype_box(axb, sub, "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                "P = 0.438", "EMT pathway score")
    panel_label(axb, "b")

    short_labels = {
        "Activated B cell": "Act. B cell",
        "Activated CD4 T cell": "Act. CD4 T",
        "Activated CD8 T cell": "Act. CD8 T",
        "Activated dendritic cell": "Act. DC",
        "CD56bright natural killer cell": "CD56bright NK",
        "CD56dim natural killer cell": "CD56dim NK",
        "Central memory CD4 T cell": "Tcm CD4",
        "Central memory CD8 T cell": "Tcm CD8",
        "Effector memeory CD4 T cell": "Tem CD4",
        "Effector memeory CD8 T cell": "Tem CD8",
        "Gamma delta T cell": "γδ T cell",
        "Immature  B cell": "Imm. B cell",
        "Immature dendritic cell": "Imm. DC",
        "Memory B cell": "Mem. B cell",
        "Natural killer T cell": "NKT cell",
        "Natural killer cell": "NK cell",
        "Plasmacytoid dendritic cell": "pDC",
        "Regulatory T cell": "Treg",
        "T follicular helper cell": "Tfh",
        "Type 1 T helper cell": "Th1",
        "Type 17 T helper cell": "Th17",
        "Type 2 T helper cell": "Th2",
    }
    imm = imm.sort_values("padj", ascending=False).reset_index(drop=True)
    row_labels = [short_labels.get(c, c) for c in imm.celltype]
    hm = imm[["IgAN", "control"]].to_numpy()
    vmax = float(np.nanmax(np.abs(hm)))
    cmap = LinearSegmentedColormap.from_list("cbo", [BLUE, "#FFFFFF", ORANGE])
    axc = fig.add_subplot(gs[:, 1])
    im = axc.imshow(hm, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto",
                    interpolation="nearest")
    axc.set_xticks([0, 1])
    axc.set_xticklabels(["IgAN (n = 55)", "Control (n = 30)"], fontsize=8)
    axc.set_xlim(-0.5, 1.75)
    axc.set_yticks(np.arange(len(imm)))
    axc.set_yticklabels(row_labels, fontsize=8)
    for r, row in imm.iterrows():
        if row.padj < 0.05:
            axc.text(1.62, r, "*", ha="center", va="center", fontsize=9,
                     color="black")
    axc.tick_params(length=0)
    for spine in axc.spines.values():
        spine.set_visible(False)
    cb = fig.colorbar(im, ax=axc, orientation="horizontal", fraction=0.045,
                      pad=0.10)
    cb.ax.tick_params(labelsize=8, width=0.6, length=2)
    cb.set_ticks([-vmax, 0.0, vmax])
    cb.set_ticklabels([f"{x:.2f}" for x in (-vmax, 0.0, vmax)])
    cb.set_label("ssGSEA score", fontsize=8)
    panel_label(axc, "c")
    return fig, []


# ---------------------------------------------------------------------------
# Fig 4（三行：a UMAP / c 空间；b 细胞类型全宽；d 时间序列全宽）
# ---------------------------------------------------------------------------

def make_fig4():
    umap = pd.read_csv(os.path.join(PANEL, "umap.csv"))
    uuo = pd.read_csv(os.path.join(PROC, "uuo_meta.csv"))
    fov = pd.read_csv(os.path.join(PANEL, "spatial_fov.csv"))
    cells = pd.read_csv(os.path.join(PANEL, "spatial_cells_sample.csv"))
    tc = pd.read_csv(os.path.join(PANEL, "timecourse_samples.csv"))

    fig = plt.figure(figsize=(FIG_W_IN, 6.55))
    gs = fig.add_gridspec(3, 3, width_ratios=[1.0, 1.0, 0.40],
                          height_ratios=[1.05, 0.95, 0.95],
                          left=0.075, right=0.985, top=0.91, bottom=0.065,
                          hspace=0.58, wspace=0.36)

    axa = fig.add_subplot(gs[0, 0])
    ctypes = ["PT", "TAL", "DCT", "CD", "Podo", "Endo", "Fibro", "Macro",
              "Tcell", "Bcell"]
    for i, ct in enumerate(ctypes):
        d = umap[umap.celltype == ct]
        axa.scatter(d.umap_1, d.umap_2, s=0.3, c=TOL_MUTED[i], lw=0,
                    label=ct, rasterized=True)
    axa.set_xlabel("UMAP 1")
    axa.set_ylabel("UMAP 2")
    axa.set_xticks([])
    axa.set_yticks([])
    handles = [Line2D([0], [0], marker="o", color="none",
                      markerfacecolor=TOL_MUTED[i], markersize=5,
                      label=ctypes[i]) for i in range(len(ctypes))]
    ax_leg = fig.add_subplot(gs[0, 2])
    ax_leg.axis("off")
    ax_leg.legend(handles=handles, loc="center left", frameon=False,
                  handletextpad=0.4, borderaxespad=0, fontsize=8,
                  title="Cell type", title_fontsize=8)
    style_axis(axa)
    panel_label(axa, "a")

    outer = fig.add_subplot(gs[0, 1])
    outer.axis("off")
    axc_cell = outer.inset_axes([0.0, 0.03, 0.40, 0.92])
    axc_cell.scatter(cells.MRG, cells.ECM, s=0.3, c=BLUE, alpha=0.08, lw=0,
                     rasterized=True)
    axc_cell.set_xlabel("Mechanical score (cell)", fontsize=7.5)
    axc_cell.set_ylabel("Fibrotic ECM score (cell)", fontsize=7.5)
    axc_cell.tick_params(labelsize=7.5, pad=2)
    axc_cell.xaxis.set_major_locator(MaxNLocator(nbins=4, prune="both"))
    axc_cell.yaxis.set_major_locator(MaxNLocator(nbins=4, prune="both"))
    axc_cell.text(0.03, 0.95, "cell-level\nρ = 0.199",
                  transform=axc_cell.transAxes, fontsize=7.5, va="top",
                  ha="left", linespacing=1.1)
    style_axis(axc_cell)
    axc_fov = outer.inset_axes([0.60, 0.03, 0.40, 0.92])
    axc_fov.scatter(fov.MRG, fov.ECM, s=4, c=ORANGE, alpha=0.6, lw=0,
                    rasterized=True)
    axc_fov.set_xlabel("Mechanical score (FOV)", fontsize=7.5)
    axc_fov.set_ylabel("Fibrotic ECM score (FOV)", fontsize=7.5)
    axc_fov.tick_params(labelsize=7.5, pad=2)
    axc_fov.xaxis.set_major_locator(MaxNLocator(nbins=4, prune="both"))
    axc_fov.yaxis.set_major_locator(MaxNLocator(nbins=4, prune="both"))
    axc_fov.text(0.03, 0.95, "FOV-level\nρ = 0.586",
                 transform=axc_fov.transAxes, fontsize=7.5, va="top",
                 ha="left", linespacing=1.1)
    style_axis(axc_fov)
    panel_label(outer, "c")

    axb = fig.add_subplot(gs[1, :])
    order = ["PT", "TAL", "DCT", "CD", "Podo", "Endo", "Fibro", "Macro",
             "Tcell", "Bcell"]
    positions = np.arange(len(order))
    for i, ct in enumerate(order):
        for gi, (grp, color) in enumerate([("sham", BLUE), ("UUO", ORANGE)]):
            v = uuo.loc[(uuo.celltype == ct) & (uuo.group == grp), "MRG1"]
            x = positions[i] + (gi - 0.5) * 0.34
            axb.boxplot(v, positions=[x], widths=0.30, patch_artist=True,
                        showfliers=False,
                        boxprops=dict(facecolor=color, alpha=0.75, lw=0.5),
                        medianprops=dict(color="black", lw=0.7),
                        whiskerprops=dict(color=color, lw=0.5),
                        capprops=dict(color=color, lw=0.5))
    axb.set_xticks(positions)
    axb.set_xticklabels(order, fontsize=8, rotation=45, ha="right")
    axb.axhline(0, color=GREY, lw=0.6, ls="--")
    axb.set_ylabel("Mechanical core module score")
    axb.set_ylim(-0.07, 0.12)
    axb.set_yticks([-0.05, 0.0, 0.05, 0.1])
    handles_b = [Line2D([0], [0], marker="s", color="none",
                        markerfacecolor=BLUE, markersize=7, label="sham"),
                 Line2D([0], [0], marker="s", color="none",
                        markerfacecolor=ORANGE, markersize=7, label="UUO")]
    axb.legend(handles=handles_b, loc="upper right", frameon=False, fontsize=8)
    axb.text(3, 0.105, "*", ha="center", fontsize=9)
    style_axis(axb)
    panel_label(axb, "b")

    axd = fig.add_subplot(gs[2, :])
    days = [0, 3, 7, 14]
    colors = [CB[0], CB[1], CB[2], CB[4]]
    for i, day in enumerate(days):
        v = tc.loc[tc.day == day, "score"].to_numpy()
        axd.boxplot(v, positions=[i], widths=0.5, patch_artist=True,
                    showfliers=False,
                    boxprops=dict(facecolor=colors[i], alpha=0.75, lw=0.5),
                    medianprops=dict(color="black", lw=0.7),
                    whiskerprops=dict(color=colors[i], lw=0.5),
                    capprops=dict(color=colors[i], lw=0.5))
        rng = np.random.default_rng(42)
        axd.scatter(i + rng.uniform(-0.08, 0.08, len(v)), v, s=5, alpha=0.7,
                    color=colors[i], lw=0)
    axd.set_xticks(range(4))
    axd.set_xticklabels(["Normal", "Day 3", "Day 7", "Day 14"])
    axd.set_xlim(-0.6, 3.6)
    axd.set_xlabel("Days after UUO")
    axd.set_ylabel("Mechanical core score")
    axd.set_ylim(-1.15, 1.45)
    axd.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    axd.text(0.02, 0.96, "Spearman ρ = 0.94\nD14 vs D0: P = 0.057",
             transform=axd.transAxes, fontsize=8, va="top")
    style_axis(axd)
    panel_label(axd, "d")

    fig.text(0.02, 0.004,
             "CD, collecting duct; highest mechanical core score among "
             "epithelial compartments in sham kidneys.",
             fontsize=7.5, va="bottom", ha="left")
    return fig, []


# ---------------------------------------------------------------------------
# 布局自检：所有文本不得越出画布；Fig1 条目文本不得越出所在方框
# ---------------------------------------------------------------------------

def layout_report(fig, name: str, pairs):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.get_size_inches() * DPI
    issues = []
    all_axes = []
    for ax in fig.axes:
        all_axes.append(ax)
        all_axes.extend(getattr(ax, "child_axes", []) or [])
    for ax in all_axes:
        texts = list(ax.texts)
        texts += list(ax.get_xticklabels()) + list(ax.get_yticklabels())
        for t in texts:
            if not t.get_visible():
                continue
            bb = t.get_window_extent(renderer=r)
            if bb.x0 < -1 or bb.y0 < -1 or bb.x1 > W + 1 or bb.y1 > H + 1:
                issues.append(
                    f"outside-canvas: '{t.get_text()[:30]}' "
                    f"({bb.x0:.0f},{bb.y0:.0f},{bb.x1:.0f},{bb.y1:.0f})")
        # 相邻刻度标签两两重叠检测
        for labels, axis in ((ax.get_yticklabels(), "y"),
                             (ax.get_xticklabels(), "x")):
            vis = [t for t in labels if t.get_visible() and t.get_text().strip()]
            vis.sort(key=lambda t: t.get_window_extent(renderer=r).y0)
            for a, b in zip(vis, vis[1:]):
                ba = a.get_window_extent(renderer=r)
                bb = b.get_window_extent(renderer=r)
                if (ba.x1 > bb.x0 + 1 and bb.x1 > ba.x0 + 1
                        and ba.y1 > bb.y0 + 1 and bb.y1 > ba.y0 + 1):
                    issues.append(
                        f"tick-overlap({axis}): '{a.get_text()[:20]}' vs "
                        f"'{b.get_text()[:20]}'")
        # 面板内注释文本两两重叠检测
        anns = [t for t in ax.texts if t.get_visible() and t.get_text().strip()]
        for i in range(len(anns)):
            for j in range(i + 1, len(anns)):
                ba = anns[i].get_window_extent(renderer=r)
                bb = anns[j].get_window_extent(renderer=r)
                if (ba.x1 > bb.x0 + 1 and bb.x1 > ba.x0 + 1
                        and ba.y1 > bb.y0 + 1 and bb.y1 > ba.y0 + 1):
                    issues.append(
                        f"text-overlap: '{anns[i].get_text()[:20]}' vs "
                        f"'{anns[j].get_text()[:20]}'")
    for (t, rect, ax) in pairs:
        bb = t.get_window_extent(renderer=r)
        x0, y0 = ax.transData.transform((rect[0], rect[1]))
        x1, y1 = ax.transData.transform((rect[0] + rect[2],
                                         rect[1] + rect[3]))
        if (bb.x0 < x0 - 1 or bb.y0 < y0 - 1
                or bb.x1 > x1 + 1 or bb.y1 > y1 + 1):
            issues.append(f"outside-box: '{t.get_text()[:30]}'")
    if issues:
        log(f"{name}: {len(issues)} layout issue(s)")
        for s in issues[:12]:
            log(f"    {s}")
    else:
        log(f"{name}: layout OK (all texts inside canvas/boxes)")
    return issues


def save_figure(fig, name: str, force: bool):
    png_path = os.path.join(PREVDIR, f"{name}.png")
    tif_path = os.path.join(OUTDIR, f"{name}.tif")
    for p in (png_path, tif_path):
        if os.path.exists(p) and not force:
            raise SystemExit(f"refusing to overwrite: {p} (use --force)")
    fig.savefig(png_path, format="png", dpi=DPI, facecolor="white")
    plt.close(fig)
    with Image.open(png_path) as im:
        rgb = np.asarray(im.convert("RGB"))
    tifffile.imwrite(tif_path, rgb, photometric="rgb", planarconfig="contig",
                     resolution=(DPI, DPI), resolutionunit="inch",
                     compression="lzw")
    with Image.open(png_path) as im:
        pw, ph = im.size
    log(f"{name}: {pw}x{ph} px PNG + LZW TIFF")
    return tif_path


def verify_tiff(path: str) -> dict:
    with Image.open(path) as im:
        info = {"size": im.size, "mode": im.mode,
                "dpi": im.info.get("dpi"), "bytes": os.path.getsize(path)}
        tag = im.tag_v2.get(259)
    return {"info": info, "tag259(LZW=5)": tag}


def main() -> int:
    ap = argparse.ArgumentParser(description="PLOS 组图重绘")
    ap.add_argument("--force", action="store_true", help="允许覆盖已有输出")
    args = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(PREVDIR, exist_ok=True)
    setup_fonts()

    builders = [("Fig1", make_fig1), ("Fig2", make_fig2),
                ("Fig3", make_fig3), ("Fig4", make_fig4)]
    total_issues = 0
    for name, fn in builders:
        fig, pairs = fn()
        issues = layout_report(fig, name, pairs)
        total_issues += len(issues)
        save_figure(fig, name, args.force)
    if total_issues:
        log(f"WARNING: {total_issues} layout issue(s) total")
    for i in range(1, 5):
        p = os.path.join(OUTDIR, f"Fig{i}.tif")
        log(f"  verify Fig{i}: {verify_tiff(p)}")
    log("DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

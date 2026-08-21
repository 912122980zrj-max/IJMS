#!/usr/bin/env python
# 57_rebuild_fig5.py —— 仅重建 Fig5（scRNA/空间/时间），不触碰其它图
# 复刻 50_figures_ijms_v2.py 的 Fig5 构建路径，用于修复面板 c 轴标签重叠。
import importlib.util
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sps
import tifffile
from PIL import Image

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"E:/sheng xin/ObstructiveNephropathy_MRG"
PROC = os.path.join(ROOT, "data", "processed")
OUTDIR = os.path.join(ROOT, "submission", "ijms", "figures")
PREVDIR = os.path.join(OUTDIR, "preview")
DPI = 300


def setup_fonts():
    from matplotlib import font_manager
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
    uuo = pd_read_uuo()
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


def pd_read_uuo():
    import pandas as pd
    return pd.read_csv(os.path.join(PROC, "uuo_meta.csv"))


def save_fig(fig, name):
    png_path = os.path.join(PREVDIR, f"{name}.png")
    pdf_path = os.path.join(OUTDIR, f"{name}.pdf")
    tif_path = os.path.join(OUTDIR, f"{name}.tif")
    fig.savefig(pdf_path, format="pdf", dpi=DPI, facecolor="white")
    fig.savefig(png_path, format="png", dpi=DPI, facecolor="white")
    with Image.open(png_path) as im:
        rgb = np.asarray(im.convert("RGB"))
    tifffile.imwrite(tif_path, rgb, photometric="rgb", planarconfig="contig",
                     resolution=(DPI, DPI), resolutionunit="inch",
                     compression="lzw")
    plt.close(fig)
    with Image.open(png_path) as im:
        print(f"Fig5: {im.size[0]}x{im.size[1]} px -> {tif_path}")


def main():
    setup_fonts()
    m = load_fig31()
    m.FIG_W_IN = 7.0
    m.panel_label = lambda ax, t: ax.text(0.0, 1.02, t,
                                          transform=ax.transAxes,
                                          fontsize=9, fontweight="bold",
                                          va="bottom", ha="left", color="black")
    fig5, _ = m.make_fig4()
    fix_fig4_ticks(fig5)
    annotate_fig5_stats(fig5)
    save_fig(fig5, "Fig5")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

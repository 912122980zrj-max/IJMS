#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""43_figures_ijms.py —— 将既有 Fig1–Fig4 按 IJMS 规格重排并重绘

相对 31_figures_plos.py 的改动（只读复用其构建函数，不改原脚本）：
  - 图宽 7.0 in（约 178 mm，MDPI 双栏版心）；
  - 正文字号下调至 7 pt（部分轴标题/图例 6.5–8 pt），面板字母 9 pt；
  - 修复 Fig1 工作流方框标题与首行条目间距过近的问题（标题与首行间距
    由 0.60 增至 0.78 数据单位）；
  - 输出：PDF（矢量）+ 300 dpi PNG 预览 + 300 dpi LZW TIFF。

用法:
    python scripts/43_figures_ijms.py [--force]
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
import tifffile
from matplotlib.figure import Figure
from PIL import Image

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"E:/sheng xin/ObstructiveNephropathy_MRG"
MOD = os.path.join(ROOT, "scripts", "31_figures_plos.py")
OUTDIR = os.path.join(ROOT, "submission", "ijms", "figures")
PREVDIR = os.path.join(OUTDIR, "preview")


def load_fig31():
    spec = importlib.util.spec_from_file_location("fig31_ijms", MOD)
    m = importlib.util.module_from_spec(spec)
    sys.modules["fig31_ijms"] = m
    spec.loader.exec_module(m)
    return m


def fixed_workflow_axis(ax):
    """与 31_figures_plos.workflow_axis 相同的三栏工作流，
    但标题与首行条目之间留出更多间距。"""
    import matplotlib.patches as mp
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    boxes = [
        (0.25, 1.2, 2.7, 8.6, "#EAF2F8", "#0072B2", "Public transcriptomes", [
            ["Human biopsies", "GSE115857 (n = 86)"],
            ["Mouse scRNA-seq", "GSE175412 (3 kidneys)"],
            ["Human spatial CosMx", "GSE282059 (523,855 cells)"],
            ["Mouse UUO time course", "GSE118339 (n = 15)"],
        ]),
        (3.65, 1.2, 2.7, 8.6, "#FFF6E6", "#E69F00", "Mechanical framework", [
            ["2,536 mechanical genes", "(3 GO categories)"],
            ["ssGSEA mechanical score"],
            ["WGCNA β = 4: 38 modules", "100-gene core (ME31)"],
            ["LASSO ∩ random forest", "NDNF, PCDHB7, RRAGB"],
        ]),
        (7.05, 1.2, 2.7, 8.6, "#EDF7F1", "#009E73", "Outcomes", [
            ["IgAN vs controls", "P = 7.01 × 10$^{-10}$"],
            ["ROC AUC 0.853", "C-index 0.848"],
            ["Subtypes C1/C2", "P = 0.0105"],
            ["CD hub: spatial ρ = 0.586", "time course ρ = 0.94"],
        ]),
    ]
    pairs = []
    for x, y, w, h, face, edge, title, items in boxes:
        ax.add_patch(mp.FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.15",
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
        ax.add_patch(mp.FancyArrowPatch((xa, 5.0), (xa + 0.68, 5.0),
                                        arrowstyle="-|>", mutation_scale=13,
                                        color="#333333", lw=1.0))
    return pairs


def fixed_panel_label(ax, text: str) -> None:
    ax.text(0.0, 1.02, text, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom", ha="left", color="black")


def fix_fig4_ticks(fig) -> None:
    """将 Fig4 FOV 散点子面板的刻度减为 3 个，避免 7 pt 下标签重叠。"""
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


_ORIG_ADD_GRIDSPEC = Figure.add_gridspec


def patch_left_margin(delta: float):
    """给 gridspec 的 left 边距加偏移。

    fig.subplots_adjust 对 add_gridspec(left=...) 创建的子图不生效，
    因此临时包装 Figure.add_gridspec，仅在本进程内为后续图形生效。
    """
    def wrapped(self, nrows=1, ncols=1, **kwargs):
        if "left" in kwargs and "right" in kwargs:
            kwargs["left"] = kwargs["left"] + delta
        return _ORIG_ADD_GRIDSPEC(self, nrows, ncols, **kwargs)

    Figure.add_gridspec = wrapped


def save_fig(fig, name: str, force: bool) -> None:
    png_path = os.path.join(PREVDIR, f"{name}.png")
    pdf_path = os.path.join(OUTDIR, f"{name}.pdf")
    tif_path = os.path.join(OUTDIR, f"{name}.tif")
    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(PREVDIR, exist_ok=True)
    for p in (png_path, pdf_path, tif_path):
        if os.path.exists(p) and not force:
            raise SystemExit(f"refusing to overwrite: {p}")
    fig.savefig(pdf_path, format="pdf", dpi=300, facecolor="white")
    fig.savefig(png_path, format="png", dpi=300, facecolor="white")
    with Image.open(png_path) as im:
        rgb = np.asarray(im.convert("RGB"))
    tifffile.imwrite(tif_path, rgb, photometric="rgb", planarconfig="contig",
                     resolution=(300, 300), resolutionunit="inch",
                     compression="lzw")
    plt.close(fig)
    with Image.open(png_path) as im:
        print(f"[figures] {name}: {im.size[0]}x{im.size[1]} px "
              f"({im.size[0] / 300 * 25.4:.1f} mm @ 300 dpi)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    m = load_fig31()
    m.FIG_W_IN = 7.0
    plt.rcParams.update({
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    m.panel_label = fixed_panel_label
    m.workflow_axis = fixed_workflow_axis

    builders = [("Fig1", m.make_fig1), ("Fig2", m.make_fig2),
                ("Fig3", m.make_fig3), ("Fig4", m.make_fig4)]
    total = 0
    for name, fn in builders:
        if name == "Fig2":
            patch_left_margin(0.020)
        elif name == "Fig4":
            patch_left_margin(0.020)
        else:
            patch_left_margin(0.0)
        fig, pairs = fn()
        if name == "Fig4":
            fix_fig4_ticks(fig)
        issues = m.layout_report(fig, name, pairs)
        total += len(issues)
        save_fig(fig, name, args.force)
    print(f"[figures] TOTAL layout issues: {total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python
# 59_rebuild_fig137.py —— 只重建需要修复的图（Fig1/3/6 由 50 脚本，Fig7 由 42 脚本）
# 不触碰 Fig2 / Fig4 / Fig5，避免覆盖已定稿版本。
import importlib.util
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from PIL import Image

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"E:/sheng xin/ObstructiveNephropathy_MRG"
OUTDIR = os.path.join(ROOT, "submission", "ijms", "figures")
PREVDIR = os.path.join(OUTDIR, "preview")
DPI = 300
CB = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"]
BLUE, ORANGE = CB[0], CB[1]


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


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
    w, h = Image.open(png_path).size
    plt.close(fig)
    print(f"{name}: {w}x{h} px -> {tif_path}")


def main():
    m = _load(os.path.join(ROOT, "scripts", "50_figures_ijms_v2.py"), "fj50")
    m.setup_fonts()
    # Fig1, Fig3, Fig6
    for name, fn in [("Fig1", m.make_fig1), ("Fig3", m.make_fig3),
                     ("Fig6", m.make_fig6)]:
        fig, _ = fn()
        save_fig(fig, name)

    # Fig7 from 42_qpcr_ijms.py
    q = _load(os.path.join(ROOT, "scripts", "42_qpcr_ijms.py"), "fp42")
    q.setup_fonts()
    long = pd.read_csv(os.path.join(ROOT, "results", "qpcr_analysis", "ijms",
                                    "qpcr_ijms_long.csv"))
    stats = pd.read_csv(os.path.join(ROOT, "results", "qpcr_analysis", "ijms",
                                     "qpcr_ijms_stats.csv"))
    fig = q.make_fig5(long, stats)
    save_fig(fig, "Fig7")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

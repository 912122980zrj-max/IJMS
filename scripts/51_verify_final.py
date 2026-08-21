#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""51_verify_final.py —— 最终核验 Fig1–Fig7 与稿件 PDF"""

import os

import fitz
import numpy as np
from PIL import Image


BASE = r"E:/sheng xin/ObstructiveNephropathy_MRG/submission/ijms"


def main() -> int:
    figdir = os.path.join(BASE, "figures")
    for i in range(1, 8):
        png = os.path.join(figdir, "preview", f"Fig{i}.png")
        tif = os.path.join(figdir, f"Fig{i}.tif")
        pdf = os.path.join(figdir, f"Fig{i}.pdf")
        ok = all(os.path.exists(p) for p in (png, tif, pdf))
        im = Image.open(png)
        t = Image.open(tif)
        a = np.asarray(im.convert("L"))
        dark = a < 100
        edge = int(dark[:, :2].sum()) + int(dark[:, -2:].sum()) \
            + int(dark[:2].sum()) + int(dark[-2:].sum())
        print(f"Fig{i} {im.size[0]}x{im.size[1]} "
              f"dpi={im.info.get('dpi')} tif={t.mode}/LZW{t.tag_v2.get(259)} "
              f"edge_dark={edge} {'OK' if ok else 'MISSING'}")
    d = fitz.open(os.path.join(BASE, "manuscript",
                              "manuscript_ijms_draft_preview.pdf"))
    print("PDF pages:", d.page_count)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

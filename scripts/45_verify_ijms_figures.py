#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""45_verify_ijms_figures.py —— 只读核验 ijms/figures 下 Fig1–Fig5 规格"""

import os

import numpy as np
from PIL import Image


BASE = r"E:/sheng xin/ObstructiveNephropathy_MRG/submission/ijms/figures"


def ink(a: np.ndarray) -> float:
    return float((np.any(a < 245, axis=2)).mean())


def main() -> int:
    for i in range(1, 6):
        png = os.path.join(BASE, "preview", f"Fig{i}.png")
        tif = os.path.join(BASE, f"Fig{i}.tif")
        pdf = os.path.join(BASE, f"Fig{i}.pdf")
        im = Image.open(png)
        t = Image.open(tif)
        arr = np.asarray(im.convert("RGB"))
        edges = {
            "T": ink(arr[:3]), "B": ink(arr[-3:]),
            "L": ink(arr[:, :3]), "R": ink(arr[:, -3:]),
        }
        ok = (im.size[0] >= 1000 and im.size[1] >= 1000
              and im.info.get("dpi") == (300.0, 300.0)
              and t.mode == "RGB" and t.tag_v2.get(259) == 5)
        print(f"Fig{i}: {im.size[0]}x{im.size[1]} px "
              f"dpi={im.info.get('dpi')} | tif {t.size} {t.mode} "
              f"dpi={t.info.get('dpi')} LZW={t.tag_v2.get(259)} | "
              f"edges={ {k: round(100 * v, 2) for k, v in edges.items()} } | "
              f"png={os.path.getsize(png)} tif={os.path.getsize(tif)} "
              f"pdf={os.path.getsize(pdf)} | {'PASS' if ok else 'FAIL'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

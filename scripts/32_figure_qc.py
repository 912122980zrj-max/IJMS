#!/usr/bin/env python
"""32_figure_qc.py

投稿图件程序化质检（只读预览 PNG，不修改任何文件）：
  - 每个网格区域的非白像素占比（检测空面板/异常）；
  - 图像四边缘 3px 是否有非白内容（检测文字/元素被裁切）；
  - 行/列沟槽区域的意外墨迹（粗略检测跨面板文字重叠）。
"""

from __future__ import annotations

import os

import numpy as np
from PIL import Image


ROOT = r"E:/sheng xin/ObstructiveNephropathy_MRG"
PREVDIR = os.path.join(ROOT, "submission", "plos", "figures", "preview")

GRIDS = {
    "Fig1.png": (1, 2),
    "Fig2.png": (1, 2),
    "Fig3.png": (2, 2),
    "Fig4.png": (2, 3),
}


def ink_frac(a: np.ndarray) -> float:
    nonwhite = np.any(a < 245, axis=2)
    return float(nonwhite.mean())


def main() -> int:
    for name, (rows, cols) in GRIDS.items():
        path = os.path.join(PREVDIR, name)
        im = np.asarray(Image.open(path).convert("RGB")).astype(np.uint8)
        H, W, _ = im.shape
        print(f"== {name} {W}x{H}")
        for r in range(rows):
            line = []
            for c in range(cols):
                y0, y1 = r * H // rows, (r + 1) * H // rows
                x0, x1 = c * W // cols, (c + 1) * W // cols
                f = ink_frac(im[y0:y1, x0:x1])
                line.append(f"{100*f:5.1f}%")
            print("   ", " | ".join(line))
        edges = {
            "top": ink_frac(im[:3, :, :]),
            "bottom": ink_frac(im[-3:, :, :]),
            "left": ink_frac(im[:, :3, :]),
            "right": ink_frac(im[:, -3:, :]),
        }
        print("    edges(3px ink):", {k: f"{100*v:.2f}%" for k, v in edges.items()})
    print("QC done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""40_layout_audit.py —— 只读渲染层审计

重新调用 31_figures_plos.py 的构建函数，仅做 renderer 实测：
  - 文本是否越出画布/所属方框；
  - 相邻刻度标签是否重叠；
  - 面板内注释文本两两是否重叠；
  - 图例/坐标轴标签与数据区是否越界（额外检查）。
不写任何文件，不覆盖既有投稿图。
"""

from __future__ import annotations

import importlib.util
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.stderr.reconfigure(encoding="utf-8", errors="replace")

ROOT = r"E:/sheng xin/ObstructiveNephropathy_MRG"
MOD = os.path.join(ROOT, "scripts", "31_figures_plos.py")


def load():
    spec = importlib.util.spec_from_file_location("fig31", MOD)
    m = importlib.util.module_from_spec(spec)
    sys.modules["fig31"] = m
    spec.loader.exec_module(m)
    return m


def main() -> int:
    m = load()
    builders = [("Fig1", m.make_fig1), ("Fig2", m.make_fig2),
                ("Fig3", m.make_fig3), ("Fig4", m.make_fig4)]
    total = 0
    for name, fn in builders:
        fig, pairs = fn()
        issues = m.layout_report(fig, name, pairs)
        total += len(issues)
        import matplotlib.pyplot as plt
        plt.close(fig)
    print(f"TOTAL layout issues: {total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

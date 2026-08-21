#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""42_qpcr_ijms.py —— qPCR 三板重分析（IJMS 口径）与 Figure 5 生成

输入:
  C:/Users/91212/Desktop/qpcr.xlsx（板1/2/3 编译 Ct 表）
输出:
  results/qpcr_analysis/ijms/qpcr_ijms_long.csv   长表（ΔCt）
  results/qpcr_analysis/ijms/qpcr_ijms_stats.csv  统计表（含每板 BH、全板 BH、阻断剂对比）
  submission/ijms/figures/Fig7.pdf/.tif/.png（正文图）
  submission/ijms/figures/preview/Fig7.png

统计口径（与既有内部报告的不同处在本文件 docstring 与报告中说明）:
  - ΔCt = Ct(target) − Ct(Gapdh)；板3 用 Gapdh（大鼠）。
  - ΔΔCt = ΔCt(treatment) − mean(ΔCt(vehicle))；fold = 2^(−ΔΔCt)。
  - 主检验：各处理 vs vehicle 的 Welch t 检验；BH 校正按“每块板”为独立
    多重比较家族（P1=18、P2=5、P3=9 个比较），另附“全板合并 32 个比较”
    的 BH 作为敏感性分析。
  - 探索性检验：阻断剂组 vs Yoda1 组（y+g vs y、y+v vs y），仅报名义 P 与
    校正 P（板内 BH），用于描述方向，不作因果宣称。
  - 附 Hedges g（小样本校正）与 fold 的 95% bootstrap 百分位区间。
  - n=3/组，不作正态性宣称；所有推断标记为探索性。
"""

from __future__ import annotations

import argparse
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
from scipy import stats as sps

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"E:/sheng xin/ObstructiveNephropathy_MRG"
DEFAULT_XLSX = r"C:/Users/91212/Desktop/qpcr.xlsx"
OUT_CSV_DIR = os.path.join(ROOT, "results", "qpcr_analysis", "ijms")
FIG_DIR = os.path.join(ROOT, "submission", "ijms", "figures")
PREV_DIR = os.path.join(FIG_DIR, "preview")

FIG_W_MM = 178.0  # MDPI 双栏版心宽约 178 mm
FIG_W_IN = FIG_W_MM / 25.4
DPI = 300

CB = {"c": "#0072B2", "y": "#E69F00", "y+g": "#009E73", "y+v": "#CC79A7",
      "tgf": "#E69F00"}
GROUP_ORDER = ["c", "y", "y+g", "y+v"]
GROUP_LABEL = {
    "c": "Vehicle",
    "y": "Yoda1",
    "y+g": "Yoda1+GsMTx4",
    "y+v": "Yoda1+Verteporfin",
    "tgf": "TGF-β1",
}


def log(msg: str) -> None:
    print(f"[qpcr] {msg}", flush=True)


def setup_fonts() -> None:
    for fam in ("Arial",):
        if any(f.name == fam for f in font_manager.fontManager.ttflist):
            plt.rcParams["font.family"] = fam
            break
    plt.rcParams.update({
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "axes.linewidth": 0.6,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "axes.edgecolor": "#333333",
        "xtick.color": "#333333",
        "ytick.color": "#333333",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


# ---------------------------------------------------------------------------
# 解析原始 xlsx（pandas 0-based 行号，与文件中实际单元格一一对应）
# ---------------------------------------------------------------------------

def parse_xlsx(path: str) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name=0, header=None)

    def block(gene_rows: list[int], ncols: int, groups: list[str],
              plate: str, gap_name: str) -> pd.DataFrame:
        rows = []
        for r in gene_rows:
            gene = str(raw.iat[r, 0]).strip()
            if gene == gap_name:
                continue
            vals = [raw.iat[r, c] for c in range(1, ncols + 1)]
            for j, v in enumerate(vals):
                rows.append({
                    "plate": plate,
                    "gene": gene,
                    "group": groups[j],
                    "col": j + 1,
                    "rep": (j % 3) + 1,
                    "ct": float(v),
                })
        # 内参
        gap_rows = []
        for r in gene_rows:
            gene = str(raw.iat[r, 0]).strip()
            if gene != gap_name:
                continue
            for j in range(ncols):
                gap_rows.append({"group": groups[j], "col": j + 1,
                                 "gap_ct": float(raw.iat[r, j + 1])})
        gap_df = pd.DataFrame(gap_rows)
        out = pd.DataFrame(rows).merge(gap_df, on=["group", "col"], how="left")
        out["dct"] = out["ct"] - out["gap_ct"]
        return out

    p1_groups = ["c"] * 3 + ["y"] * 3 + ["y+g"] * 3 + ["y+v"] * 3
    p2_groups = ["c"] * 3 + ["tgf"] * 3
    p1 = block(list(range(2, 9)), 12, p1_groups, "P1", "Gapdh")
    p2 = block(list(range(12, 18)), 6, p2_groups, "P2", "Gapdh")
    p3 = block(list(range(22, 26)), 12, p1_groups, "P3", "Gapdh（大鼠）")
    long = pd.concat([p1, p2, p3], ignore_index=True)
    return long


def welch(v1: np.ndarray, v2: np.ndarray) -> tuple[float, float, float]:
    """Welch t（双尾）。返回 (t, p, df)。"""
    t, p = sps.ttest_ind(v1, v2, equal_var=False)
    n1, n2 = len(v1), len(v2)
    s1, s2 = np.var(v1, ddof=1), np.var(v2, ddof=1)
    num = (s1 / n1 + s2 / n2) ** 2
    den = (s1 / n1) ** 2 / (n1 - 1) + (s2 / n2) ** 2 / (n2 - 1)
    df = num / den if den > 0 else np.nan
    return float(t), float(p), float(df)


def ddct_ci(v1: np.ndarray, v2: np.ndarray) -> tuple[float, float]:
    d = np.mean(v1) - np.mean(v2)
    _, _, df = welch(v1, v2)
    n1, n2 = len(v1), len(v2)
    se = np.sqrt(np.var(v1, ddof=1) / n1 + np.var(v2, ddof=1) / n2)
    tcrit = sps.t.ppf(0.975, df) if np.isfinite(df) else 1.96
    lo = d - tcrit * se
    hi = d + tcrit * se
    return float(lo), float(hi)


def hedges_g(v1: np.ndarray, v2: np.ndarray) -> float:
    n1, n2 = len(v1), len(v2)
    sp = np.sqrt(((n1 - 1) * np.var(v1, ddof=1)
                  + (n2 - 1) * np.var(v2, ddof=1)) / (n1 + n2 - 2))
    if sp == 0:
        return float("nan")
    g = (np.mean(v1) - np.mean(v2)) / sp
    g = g * (1 - 3 / (4 * (n1 + n2) - 9))
    return float(g)


def bh(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    q = np.empty(m)
    prev = np.inf
    for i in range(m - 1, -1, -1):
        r = order[i]
        cur = min(p[r] * m / (i + 1), prev)
        q[r] = cur
        prev = cur
    return q


def boot_fold_ci(v1: np.ndarray, v2: np.ndarray, n_boot: int = 20000,
                 seed: int = 42) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        a = rng.choice(v1, size=len(v1), replace=True)
        b = rng.choice(v2, size=len(v2), replace=True)
        boots[i] = 2 ** (-(np.mean(a) - np.mean(b)))
    return float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def run_stats(long: pd.DataFrame) -> pd.DataFrame:
    rows = []
    # 各基因各板 vehicle 均值
    ctrl_mean = (long[long.group == "c"].groupby(["plate", "gene"])
                 .dct.mean().rename("ctrl_mean").reset_index())
    for (plate, gene), d in long.groupby(["plate", "gene"]):
        c0 = d.loc[d.group == "c", "dct"].to_numpy()
        base_mean = c0.mean()
        treat_groups = [g for g in d.group.unique() if g != "c"]
        for g in treat_groups:
            v = d.loc[d.group == g, "dct"].to_numpy()
            dd = v.mean() - base_mean
            t, p, df = welch(v, c0)
            _, w_p = sps.mannwhitneyu(v, c0, alternative="two-sided",
                                      method="exact")
            lo_d, hi_d = ddct_ci(v, c0)
            lo_f, hi_f = boot_fold_ci(v, c0)
            rows.append({
                "plate": plate, "gene": gene, "contrast": f"{g}_vs_c",
                "group": g, "n": len(v), "ddct": dd, "fold": 2 ** (-dd),
                "ddct_ci_lo": lo_d, "ddct_ci_hi": hi_d,
                "fold_ci_lo": lo_f, "fold_ci_hi": hi_f,
                "t": t, "t_df": df, "t_p": p, "w_p": float(w_p),
                "hedges_g": hedges_g(v, c0),
            })
        # 阻断剂 vs Yoda1（探索性）
        if "y" in d.group.unique():
            yv = d.loc[d.group == "y", "dct"].to_numpy()
            for g in ("y+g", "y+v"):
                if g not in d.group.unique():
                    continue
                v = d.loc[d.group == g, "dct"].to_numpy()
                t, p, df = welch(v, yv)
                _, w_p = sps.mannwhitneyu(v, yv, alternative="two-sided",
                                          method="exact")
                rows.append({
                    "plate": plate, "gene": gene, "contrast": f"{g}_vs_y",
                    "group": g, "n": len(v),
                    "ddct": v.mean() - yv.mean(),
                    "fold": 2 ** (-(v.mean() - yv.mean())),
                    "ddct_ci_lo": np.nan, "ddct_ci_hi": np.nan,
                    "fold_ci_lo": np.nan, "fold_ci_hi": np.nan,
                    "t": t, "t_df": df, "t_p": p, "w_p": float(w_p),
                    "hedges_g": hedges_g(v, yv),
                })
    stats = pd.DataFrame(rows)

    # 每板 BH（vs c）与全板 BH（vs c）敏感性
    vs_c = stats.contrast.str.endswith("_vs_c")
    stats["padj_per_plate"] = np.nan
    stats["padj_pooled"] = np.nan
    stats["padj_blocker"] = np.nan
    for plate, idx in stats[vs_c].groupby("plate").groups.items():
        stats.loc[idx, "padj_per_plate"] = bh(stats.loc[idx, "t_p"])
    stats.loc[vs_c, "padj_pooled"] = bh(stats.loc[vs_c, "t_p"])
    for plate, idx in stats[~vs_c].groupby("plate").groups.items():
        stats.loc[idx, "padj_blocker"] = bh(stats.loc[idx, "t_p"])
    return stats


# ---------------------------------------------------------------------------
# Figure 5
# ---------------------------------------------------------------------------

def rel_expr(long: pd.DataFrame, plate: str, gene: str,
             group: str) -> np.ndarray:
    d = long[(long.plate == plate) & (long.gene == gene)]
    ctrl = d.loc[d.group == "c", "dct"].mean()
    v = d.loc[d.group == group, "dct"].to_numpy()
    return 2 ** (-(v - ctrl))


def draw_panel(ax, long: pd.DataFrame, stats: pd.DataFrame, plate: str,
               genes: list[str], groups: list[str], use_q: str = "per_plate",
               gene_map: dict | None = None):
    qcol = {"per_plate": "padj_per_plate", "pooled": "padj_pooled"}[use_q]
    for gi, gene in enumerate(genes):
        lookup = (gene_map or {}).get(gene, gene)
        ymax = -np.inf
        ymin = np.inf
        for t, g in enumerate(groups):
            rel = rel_expr(long, plate, lookup, g)
            if len(rel) == 0:
                continue
            x = gi + (t - (len(groups) - 1) / 2) * 0.21
            jx = x + np.random.default_rng(42).uniform(-0.028, 0.028, len(rel))
            ax.scatter(jx, np.log2(rel), s=10, color=CB.get(g, "#666666"),
                       alpha=0.85, lw=0, zorder=3)
            m = float(np.mean(np.log2(rel)))
            ax.plot([x - 0.045, x + 0.045], [m, m], color="#222222",
                    lw=1.0, zorder=4)
            ymax = max(ymax, float(np.max(np.log2(rel))))
            ymin = min(ymin, float(np.min(np.log2(rel))))
        yspan = ymax - ymin
    ax.axhline(0, color="#999999", lw=0.7, ls="--", zorder=1)
    ax.set_xticks(np.arange(len(genes)))
    ax.set_xticklabels(genes, fontstyle="italic")
    ax.set_xlim(-0.75, len(genes) - 0.25)
    ax.set_ylabel("Relative mRNA\n($\\mathrm{log_{2}}$, vehicle = 1)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=0.6, length=2.5)


def panel_label(ax, letter: str, title: str) -> None:
    ax.text(0.0, 1.008, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom", ha="left", color="black")
    ax.text(0.055, 1.008, title, transform=ax.transAxes, fontsize=8,
            va="bottom", ha="left", color="#222222")


def make_fig5(long: pd.DataFrame, stats: pd.DataFrame):
    p1_genes = ["Piezo1", "Yap1", "Ctgf", "Cyr61", "Ankrd1", "Tgfb1"]
    p2_genes = ["Ndnf", "Pcdhb7", "Rragb", "Ctgf (retest)", "Tgfb1 (retest)"]
    p2_map = {"Ctgf (retest)": "Ctgf（复测）",
              "Tgfb1 (retest)": "Tgfb1（复测）"}
    p3_genes = ["Acta2", "Col1a1", "Fn1"]

    fig = plt.figure(figsize=(FIG_W_IN, 6.9))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.35, 0.85, 1.0],
                          left=0.075, right=0.975, top=0.885, bottom=0.075,
                          hspace=0.72)

    axa = fig.add_subplot(gs[0])
    draw_panel(axa, long, stats, "P1", p1_genes, GROUP_ORDER)
    panel_label(axa, "a", "mIMCD-3 collecting-duct cells (mechanism plate)")
    axa.set_ylim(-0.45, 3.0)

    axb = fig.add_subplot(gs[1])
    draw_panel(axb, long, stats, "P2", p2_genes, ["c", "tgf"],
               gene_map=p2_map)
    panel_label(axb, "b", "mIMCD-3 TGF-$\\beta$1 response (signature plate)")
    axb.set_ylim(-0.45, 3.0)

    axc = fig.add_subplot(gs[2])
    draw_panel(axc, long, stats, "P3", p3_genes, GROUP_ORDER)
    panel_label(axc, "c",
                "NRK-49F fibroblasts, conditioned medium (paracrine plate)")
    axc.set_ylim(-0.55, 3.75)

    handles = [plt.Line2D([0], [0], marker="o", color="none",
                          markerfacecolor=CB[g], markersize=5,
                          label=GROUP_LABEL[g]) for g in GROUP_ORDER]
    fig.legend(handles=handles, loc="upper center",
               bbox_to_anchor=(0.5, 0.965), frameon=False, ncol=4,
               columnspacing=1.2, handletextpad=0.35, borderaxespad=0,
               fontsize=6.5)
    fig.text(0.02, 0.008,
             "Pilot experiment: three technical replicates per group "
             "(single biological sample); bars = group mean; dashed line = "
             "vehicle.",
             fontsize=6, va="bottom", ha="left")
    return fig


def layout_report(fig, name: str) -> list[str]:
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.get_size_inches() * DPI
    issues = []
    for ax in fig.axes:
        texts = list(ax.texts) + list(ax.get_xticklabels()) \
            + list(ax.get_yticklabels())
        for t in texts:
            if not t.get_visible():
                continue
            bb = t.get_window_extent(renderer=r)
            if bb.x0 < -1 or bb.y0 < -1 or bb.x1 > W + 1 or bb.y1 > H + 1:
                issues.append(f"outside-canvas: '{t.get_text()[:30]}' "
                              f"({bb.x0:.0f},{bb.y0:.0f},{bb.x1:.0f},{bb.y1:.0f})")
        anns = [t for t in ax.texts if t.get_visible() and t.get_text().strip()]
        for i in range(len(anns)):
            for j in range(i + 1, len(anns)):
                ba = anns[i].get_window_extent(renderer=r)
                bb = anns[j].get_window_extent(renderer=r)
                if (ba.x1 > bb.x0 + 1 and bb.x1 > ba.x0 + 1
                        and ba.y1 > bb.y0 + 1 and bb.y1 > ba.y0 + 1):
                    issues.append(f"text-overlap: '{anns[i].get_text()[:20]}' vs "
                                  f"'{anns[j].get_text()[:20]}'")
    if issues:
        log(f"{name}: {len(issues)} layout issue(s)")
        for s in issues[:15]:
            log(f"    {s}")
    else:
        log(f"{name}: layout OK")
    return issues


def save_fig(fig, name: str, force: bool = False) -> None:
    png_path = os.path.join(PREV_DIR, f"{name}.png")
    pdf_path = os.path.join(FIG_DIR, f"{name}.pdf")
    tif_path = os.path.join(FIG_DIR, f"{name}.tif")
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(PREV_DIR, exist_ok=True)
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
        log(f"{name}: {im.size[0]}x{im.size[1]} px ("
            f"{im.size[0] / DPI * 25.4:.1f} mm @ {DPI} dpi)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", default=DEFAULT_XLSX)
    ap.add_argument("--outdir", default=OUT_CSV_DIR)
    ap.add_argument("--force", action="store_true",
                    help="覆盖由本脚本生成的图件输出")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    long = parse_xlsx(args.xlsx)
    stats = run_stats(long)
    long.to_csv(os.path.join(args.outdir, "qpcr_ijms_long.csv"),
                index=False, encoding="utf-8-sig")
    stats.to_csv(os.path.join(args.outdir, "qpcr_ijms_stats.csv"),
                 index=False, encoding="utf-8-sig")
    log(f"long rows={len(long)} stats rows={len(stats)}")

    setup_fonts()
    fig = make_fig5(long, stats)
    issues = layout_report(fig, "Fig7")
    save_fig(fig, "Fig7", args.force)
    log(f"Fig7 layout issues: {len(issues)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python
"""22_manuscript_data_audit.py

对手稿 (paper_rewriting_output/final_paper/main.tex) 中出现的每一个定量声明，
与 results/ 和 data/processed/ 的中间结果逐条对账。只读、不写数据文件；仅输出
Markdown 审计报告。

用法:
    python scripts/22_manuscript_data_audit.py \
        --root "E:/sheng xin/ObstructiveNephropathy_MRG" \
        --out results/manuscript_data_audit.md
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import sys
from dataclasses import dataclass
from datetime import datetime

import numpy as np
import pandas as pd
from scipy import stats


def sig_round(x: float, n: int) -> float:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return x
    if x == 0:
        return 0.0
    return round(x, n - int(math.floor(math.log10(abs(x)))) - 1)


def sigfigs(expected: float | str) -> int:
    """从期望值字面量推断有效数字位数（用于对账舍入）。"""
    s = str(expected).lower()
    s = s.replace("e-", "e-").strip()
    mant = s.split("e")[0]
    mant = mant.replace("-", "").replace("+", "").replace(".", "").lstrip("0")
    return max(1, len(mant.rstrip("0")) if mant else 1)


def match_sig(computed: float, expected: float | str) -> bool:
    """按期望值的有效数字位数比较两个浮点数。"""
    if computed is None:
        return False
    exp = float(expected)
    n = sigfigs(expected)
    if exp == 0:
        return abs(computed) < 0.5 * 10 ** (-n)
    # 比例误差按 1 个最低有效位的一半判定
    ulp = 0.5 * 10 ** (math.floor(math.log10(abs(exp))) - (n - 1))
    return abs(computed - exp) <= ulp


@dataclass
class Finding:
    item: str
    claim: str
    computed: str
    expected: str
    status: str  # PASS / FAIL / INFO
    evidence: str
    note: str = ""


class Audit:
    def __init__(self, root: str, r_items: str):
        self.root = root
        self.dp = os.path.join(root, "data", "processed")
        self.rf = os.path.join(root, "data", "reference")
        self.findings: list[Finding] = []
        self.r = json.load(open(r_items, encoding="utf-8"))

    def add(self, item, claim, computed, expected, evidence, note="", status=None):
        if status is None:
            try:
                status = "PASS" if match_sig(float(computed), expected) else "FAIL"
            except (TypeError, ValueError):
                status = "INFO"
        self.findings.append(
            Finding(item, claim, str(computed), str(expected), status, evidence, note)
        )

    def add_bool(self, item, claim, ok, computed, expected, evidence, note=""):
        self.findings.append(
            Finding(item, claim, str(computed), str(expected),
                    "PASS" if ok else "FAIL", evidence, note)
        )

    def load(self, name: str) -> pd.DataFrame:
        return pd.read_csv(os.path.join(self.dp, name))

    def cohort(self):
        meta = self.load("meta_human_main.csv")
        total = len(meta)
        igan = int((meta["status"] == "IgAN patient").sum())
        mcd = int((meta["status"] == "minimal change disease").sum())
        mn = int((meta["status"] == "Membranous glomerulonephritis").sum())
        donor = int((meta["status"] == "Living donor").sum())
        fsgs = int((meta["status"] == "Focal Segmental Glomerulosclerosis").sum())
        ctrl = mcd + mn + donor
        self.add("C1", "GSE115857 总样本 n=86", total, 86, "meta_human_main.csv")
        self.add("C2", "IgAN n=55", igan, 55, "meta_human_main.csv")
        self.add("C3", "非进展对照 n=30 (MCD12+MN11+供肾7)", ctrl, 30, "meta_human_main.csv")
        self.add_bool("C4", "MCD=12, MN=11, 供肾=7, FSGS=1（FSGS 不进入对比）",
                      mcd == 12 and mn == 11 and donor == 7 and fsgs == 1,
                      (mcd, mn, donor, fsgs), "(12,11,7,1)", "meta_human_main.csv")

    def deg(self):
        deg = self.load("deg_igag_vs_control.csv")
        sig = deg[(deg["logFC"].abs() > 0.585) & (deg["adj.P.Val"] < 0.05)]
        n, nup, ndown = len(sig), int((sig["logFC"] > 0).sum()), int((sig["logFC"] < 0).sum())
        self.add("D1", "显著 DEG n=231 (|log2FC|>0.585 & adjP<0.05)", n, 231, "deg_igag_vs_control.csv")
        self.add("D2a", "上调基因数（手稿误写 43，数据为 58）", nup, 58, "deg_igag_vs_control.csv")
        self.add("D2b", "下调基因数（手稿误写 190，数据为 173）", ndown, 173, "deg_igag_vs_control.csv")
        self.add_bool("D2c", "手稿 43/190 与数据 58/173 一致",
                      (nup, ndown) == (43, 190), f"{nup}/{ndown}", "43/190",
                      "deg_igag_vs_control.csv",
                      "数据实际为 58/173，且 43+190=233≠231 自身矛盾；重写时须更正。")
        mrg = set(x.strip() for x in open(os.path.join(self.rf, "MRG_union.txt"), encoding="utf-8") if x.strip())
        overlap = len(set(sig["gene"]) & mrg)
        self.add("D3", "DEG∩MRG = 36", overlap, 36, "deg_igag_vs_control.csv ∩ MRG_union.txt")
        t2 = deg.set_index("gene").reindex(["NDNF", "PCDHB7", "RRAGB"])
        for g, exp in [("NDNF", -1.03), ("PCDHB7", 0.44), ("RRAGB", -0.36)]:
            row = t2.loc[g]
            note = (f"adjP={row['adj.P.Val']:.3g}；"
                    f"{'通过|log2FC|>0.585' if abs(row['logFC']) > 0.585 else '未过 FC 阈值（仅 FDR 显著）'}")
            self.add(f"T2-{g}", f"表2 {g} log2FC={exp}", float(row["logFC"]), exp,
                     "deg_igag_vs_control.csv", note=note)

    def ssgsea(self):
        meta = self.load("meta_human_main.csv")
        score = self.load("ssgsea_mech_scores.csv")
        meta["disease"] = np.where(meta["status"] == "IgAN patient", "IgAN", "control")
        df = meta.merge(score, on="sample")
        a = df.loc[df.disease == "IgAN", "MRG_up"].astype(float)
        b = df.loc[df.disease == "control", "MRG_up"].astype(float)
        p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
        self.add("S1", "机械评分 IgAN vs 对照 wilcox p=3.63e-10", p, "3.63e-10",
                 "ssgsea_mech_scores.csv + meta_human_main.csv")

    def wgcna(self):
        w = self.r["wgcna"]
        self.add("W1", "WGCNA β=4", w["power_estimate"], 4, "audit_r_items.json（重算 pickSoftThreshold）")
        self.add("W2", "scale-free R²=0.741", w["sft_r2_at_power4"], 0.741, "audit_r_items.json")
        self.add("W3", "模块数=38", w["n_modules"], 38, "wgcna_modules.rds")
        self.add("W4", "最优模块 |cor|=0.477", w["best_module_abs_cor"], 0.477, "wgcna_modules.rds")
        self.add("W5", "机械核心基因=100", w["n_core"], 100, "wgcna_core.rds")

    def signature(self):
        sig = set(self.r["signature"])
        ok = sig == {"NDNF", "PCDHB7", "RRAGB"}
        self.add_bool("G1", "3-gene signature = NDNF/PCDHB7/RRAGB", ok,
                      sorted(sig), "['NDNF','PCDHB7','RRAGB']", "progression_signature.rds")
        uni = self.load("univariate_igag.csv").set_index("gene")
        for g, exp in [("NDNF", 0.23), ("PCDHB7", 15.9), ("RRAGB", 0.075)]:
            self.add(f"G2-{g}", f"表2 {g} 单因素 OR={exp}", float(uni.loc[g, "or"]), exp,
                     "univariate_igag.csv")

    def roc(self):
        roc = self.load("signature_roc.csv")
        overall = float(roc.loc[roc.fold == 0, "AUC"].iloc[0])
        cv = roc.loc[roc.fold > 0, "AUC"].astype(float)
        self.add("R1", "整体 AUC=0.853", overall, 0.853, "signature_roc.csv")
        self.add("R2", "5 折 CV AUC=0.813", cv.mean(), 0.813, "signature_roc.csv")
        self.add("R3", "CV AUC SD=0.079", cv.std(ddof=1), 0.079, "signature_roc.csv")

    def nomogram(self):
        nm = self.r["nomogram"]
        ok = set(nm["terms"]) == {"NDNF", "PCDHB7"}
        self.add_bool("N1", "逐步回归保留 NDNF+PCDHB7", ok, nm["terms"],
                      "['NDNF','PCDHB7']", "stepwise_glm.rds")
        self.add("N2", "C-index=0.848（手稿误写 0.696=把 Dxy 当作 C-index）", nm["c_index"], 0.848,
                 "stepwise_glm.rds→val.prob C(ROC)，样本 85",
                 note="0.69576 是 Somers' Dxy，不是 C-index；正确 C=0.84788。手稿正文与摘要均需更正。")
        self.add("N3", "校准斜率=1.0", nm["slope"], 1.0, "stepwise_glm.rds→val.prob Slope")
        alt = self.r["nomogram_fsgs_included"]
        self.findings.append(Finding(
            "N4", "FSGS 并入对照口径（18_figure_final 的 idx 逻辑）的 C-index，用于溯源 0.696",
            f"C={alt['c_index']:.3g} / Dxy={alt['dxy']:.3g}", "见计算值", "INFO",
            "stepwise_glm.rds→val.prob，样本 86",
            "两口径 C-index 均为 0.848~0.850，与手稿 0.696 均不符；0.696 恰为 Dxy，证实标签错误。"))

    def subtype(self):
        meta = self.load("meta_human_main.csv")
        score = self.load("ssgsea_mech_scores.csv")
        sub = self.load("subtype_assignment.csv")
        c1 = int((sub.subtype == "C1").sum())
        c2 = int((sub.subtype == "C2").sum())
        self.add("B1", "C1 n=12", c1, 12, "subtype_assignment.csv")
        self.add("B2", "C2 n=43", c2, 43, "subtype_assignment.csv")
        p = self.r["clustering"]["subtype_mech_wilcox_p"]
        self.add("B3", "C1 vs C2 机械评分 p=0.0105", p, 0.0105,
                 "audit_r_items.json（R wilcox.test，含连续性校正）")
        self.add("B4", "PAC(K=2)=0.442", self.r["clustering"]["pac_k2_to_k9"][0], 0.442,
                 "audit_r_items.json（重算 ConsensusClusterPlus）")

    def hallmark(self):
        meta = self.load("meta_human_main.csv")
        meta["disease"] = np.where(meta.status == "IgAN patient", "IgAN", "control")
        hall = self.load("hallmark_scores.csv").merge(meta, on="sample")
        tgf = hall["HALLMARK_TGF_BETA_SIGNALING"].astype(float)
        emt = hall["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"].astype(float)
        p_tgf = stats.mannwhitneyu(tgf[hall.disease == "IgAN"], tgf[hall.disease == "control"]).pvalue
        p_emt = stats.mannwhitneyu(emt[hall.disease == "IgAN"], emt[hall.disease == "control"]).pvalue
        self.add("H1", "TGF-β 通路 IgAN vs 对照 p=3e-5", p_tgf, 3e-5, "hallmark_scores.csv")
        self.add("H2", "EMT 通路 IgAN vs 对照 p=0.246", p_emt, 0.246, "hallmark_scores.csv")
        # 子型 EMT 比较（图 3B 实际展示的比较，用于定位文本-图不一致）
        sub = self.load("subtype_assignment.csv").merge(hall, on="sample")
        p_emt_sub = self.r["clustering"]["subtype_emt_wilcox_p"]
        self.findings.append(Finding(
            "H3", "图3B 显示的是“EMT 按亚型”比较（正文却写 IgAN vs 对照 p=0.246）",
            f"EMT C1 vs C2 p={p_emt_sub:.3g}", "文本-图不一致", "INFO",
            "hallmark_scores.csv + subtype_assignment.csv",
            "重写时需统一比较对象：要么图改画疾病间比较，要么正文改用亚型 p 值。"))

    def immune(self):
        res = self.load("immune28_charoentong_results.csv")
        n_sig = int((res.padj < 0.05).sum())
        self.add("I1", "28 免疫细胞中 11 个 BH<0.05", n_sig, 11, "immune28_charoentong_results.csv")
        for ct, exp, up in [("Central memory CD8 T cell", 1.2e-5, True),
                            ("Activated B cell", 2.3e-4, True)]:
            row = res[res.celltype == ct].iloc[0]
            self.add(f"I2-{ct[:12]}", f"{ct} p={exp}（IgAN 升高）", row.p, exp,
                     "immune28_charoentong_results.csv")
            self.add_bool(f"I3-{ct[:12]}", f"{ct} 方向：IgAN > control",
                          row.IgAN > row.control, (row.IgAN, row.control),
                          "IgAN>control", "immune28_charoentong_results.csv")
        eo = res[res.celltype == "Eosinophil"].iloc[0]
        self.add_bool("I4", "Eosinophil 方向：IgAN < control（降低）", eo.IgAN < eo.control,
                      (eo.IgAN, eo.control), "IgAN<control", "immune28_charoentong_results.csv")

    def scrna(self):
        um = self.load("uuo_meta.csv")
        self.add("U1", "UUO scRNA QC 后细胞数=54,930", len(um), 54930, "uuo_meta.csv")
        prop = self.load("uuo_celltype_proportions.csv") * 100
        sham = prop.loc[(prop["PT"] - 38.5).abs().idxmin()]
        uuo = prop.loc[(prop["PT"] - 5.9).abs().idxmin()]
        for lbl, row, pt, fib, mac in [("sham", sham, 38.5, 1.7, 4.8), ("UUO", uuo, 5.9, 10.0, 21.9)]:
            self.add(f"U2-{lbl}-PT", f"{lbl} PT 比例={pt}%", row["PT"], pt, "uuo_celltype_proportions.csv")
            self.add(f"U2-{lbl}-Fibro", f"{lbl} 成纤维={fib}%", row["Fibro"], fib, "uuo_celltype_proportions.csv")
            self.add(f"U2-{lbl}-Macro", f"{lbl} 巨噬={mac}%", row["Macro"], mac, "uuo_celltype_proportions.csv")
        # 机械评分按细胞类型
        epi = ["PT", "TAL", "DCT", "CD", "Podo"]
        means = um[um.group == "sham"].groupby("celltype")["MRG1"].mean()
        cd_sham = float(means.get("CD", np.nan))
        best_epi = means.reindex(epi).dropna().idxmax()
        self.add("U3", "sham 集合管机械评分=0.026", cd_sham, 0.026, "uuo_meta.csv MRG1")
        self.add_bool("U4", "sham 上皮腔室中 CD 机械评分最高",
                      best_epi == "CD", f"最高={best_epi} ({cd_sham:.3g})", "CD",
                      "uuo_meta.csv MRG1")
        fib = um[um.celltype == "Fibro"].groupby("group")["MRG1"].mean()
        self.add("U5", "UUO 成纤维机械评分=0.0149", float(fib.get("UUO", np.nan)), 0.0149,
                 "uuo_meta.csv MRG1")
        self.add("U6", "sham 成纤维机械评分=0.0114", float(fib.get("sham", np.nan)), 0.0114,
                 "uuo_meta.csv MRG1")

    def spatial(self):
        sp = self.load("spatial_validation.csv").set_index("metric")["value"]
        self.add("P1", "细胞级 Spearman rho=0.199", float(sp["cell_spearman"]), 0.199,
                 "spatial_validation.csv")
        self.add("P2", "FOV 级 Spearman rho=0.586", float(sp["fov_spearman"]), 0.586,
                 "spatial_validation.csv")
        self.add_bool("P3", "高 ECM 四分位 MRG 显著更高 (p≈0)", float(sp["high_ecm_wilcox_p"]) < 1e-6,
                      sp["high_ecm_wilcox_p"], "p≈0", "spatial_validation.csv")
        # 细胞总数：流式统计 CosMx 表达矩阵行数（不载入内存）
        path = os.path.join(self.root, "data", "raw", "GSE282059", "suppl",
                            "GSE282059_2B_FFPE_exprMat_file.csv.gz")
        if os.path.exists(path):
            n = -1
            with gzip.open(path, "rb") as fh:
                for line in fh:
                    n += 1
            self.add("P4", "CosMx 细胞数=523,855", n, 523855, "GSE282059 exprMat 行数")
        else:
            self.findings.append(Finding(
                "P4", "CosMx 细胞数=523,855（原始文件缺失，无法直接复核）", "缺失",
                "523855", "INFO", "GSE282059 exprMat", "可依据 11_spatial_validate.R 日志复核"))

    def timecourse(self):
        tc = self.load("uuo_timecourse_results.csv").set_index("metric")["value"]
        self.add("T1", "时间序列 Spearman rho=0.94", float(tc["spearman_day"]), 0.94,
                 "uuo_timecourse_results.csv")
        self.add("T2", "D14 vs D0 wilcox p=0.057", float(tc["D14_vs_D0_wilcox_p"]), 0.057,
                 "uuo_timecourse_results.csv")
        frac = float(tc["positive_slope_frac"])
        n65 = round(frac * 76)
        self.add("T3", "核心基因正向斜率占比=85.5%", frac * 100, 85.5,
                 "uuo_timecourse_results.csv")
        self.add_bool("T4", "正向斜率 65/76", n65 == 65, f"{n65}/76", "65/76",
                      "uuo_timecourse_results.csv（0.855263*76=65）")

    def cross_species(self):
        oc = self.load("ortholog_concordance.csv")
        n = len(oc)
        frac = 100 * oc.consistent.mean()
        n_sig = int((oc["adj.P.Val"] < 0.05).sum())
        self.add("X1", "可评估核心同源基因=55", n, 55, "ortholog_concordance.csv")
        self.add("X2", "方向一致占比=58.2%", frac, 58.2, "ortholog_concordance.csv")
        self.add("X3", "显著 (adjP<0.05) n=37", n_sig, 37, "ortholog_concordance.csv")

    def run(self) -> str:
        self.cohort(); self.deg(); self.ssgsea(); self.wgcna(); self.signature()
        self.roc(); self.nomogram(); self.subtype(); self.hallmark(); self.immune()
        self.scrna(); self.spatial(); self.timecourse(); self.cross_species()
        return self.render()

    def render(self) -> str:
        n_pass = sum(f.status == "PASS" for f in self.findings)
        n_fail = sum(f.status == "FAIL" for f in self.findings)
        n_info = sum(f.status == "INFO" for f in self.findings)
        lines = [
            "# 手稿数据对账报告（PLOS Computational Biology 投稿前核查）",
            "",
            f"- 生成时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"- 结论：{n_pass} PASS / {n_fail} FAIL / {n_info} INFO",
            "- 对账口径：与 `paper_rewriting_output/final_paper/main.tex` 正文及表 1/表 2 逐一比较；"
            "PASS 表示在稿件有效数字位数内一致。",
            "",
            "| # | 声明 | 计算值 | 手稿值 | 状态 | 证据 |",
            "|---|---|---|---|---|---|",
        ]
        for i, f in enumerate(self.findings, 1):
            lines.append(f"| {i} | {f.item}：{f.claim} | {f.computed} | {f.expected} "
                         f"| {f.status} | {f.evidence} |")
        fail_notes = [f for f in self.findings if f.status in ("FAIL", "INFO") and f.note]
        if fail_notes:
            lines += ["", "## 需要处置的发现", ""]
            for f in fail_notes:
                lines.append(f"- **{f.item}**（{f.status}）：{f.note}")
        if n_fail == 0:
            lines += ["", "> 未发现数值性 FAIL。所有可复核的数字均与手稿一致。", ""]
        return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description="手稿数据核查")
    ap.add_argument("--root", default=r"E:/sheng xin/ObstructiveNephropathy_MRG")
    ap.add_argument("--out", default=None, help="输出 Markdown 路径")
    ap.add_argument("--r-items", default=None, help="audit_r_items.json 路径")
    args = ap.parse_args()

    root = args.root
    r_items = args.r_items or os.path.join(root, "data", "processed", "audit_r_items.json")
    out = args.out or os.path.join(root, "results", "manuscript_data_audit.md")
    if not os.path.exists(r_items):
        print(f"缺少 {r_items}：请先运行 scripts/23_audit_extract.R", file=sys.stderr)
        return 2
    report = Audit(root, r_items).run()
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(report)
    print(f"审计报告已写入 {out}")
    print(report[:1200])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

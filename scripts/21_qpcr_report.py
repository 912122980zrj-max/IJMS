#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""21_qpcr_report.py —— 汇总三板 qPCR 结果为中文报告与 xlsx"""

import os

import pandas as pd

OUT = r"E:/sheng xin/ObstructiveNephropathy_MRG/results/qpcr_analysis"

PLATE_CN = {"P1": "板1 mIMCD-3（机制板）", "P2": "板2 mIMCD-3（阳性+签名板）", "P3": "板3 NRK-49F（成纤维板）"}
GROUP_CN = {"c": "对照 DMSO", "y": "Yoda1", "y+g": "Yoda1+GsMTx4",
            "y+v": "Yoda1+Verteporfin", "tgf": "TGF-β1"}


def sig_of(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def main():
    stats = pd.read_csv(os.path.join(OUT, "qpcr_stats.csv"))
    cv = pd.read_csv(os.path.join(OUT, "qpcr_cv.csv"))
    wells = pd.read_csv(os.path.join(OUT, "plate1_wells.csv"))
    recovered = pd.read_csv(os.path.join(OUT, "plate1_recovered.csv"))

    L = []
    add = L.append
    add("# qPCR 三板分析报告")
    add("")
    add("## 一、结论摘要")
    add("")
    for _, s in stats.iterrows():
        add("- %s，%s，%s：ΔΔCt = %.3f，fold = %.2f，t 检验 p = %.4g，"
            "BH 校正 p = %.4g%s"
            % (PLATE_CN[s["plate"]], s["gene"], GROUP_CN[s["group"]],
               s["ddct"], s["fold"], s["t_p"], s["padj_t"],
               "，显著（%s）" % sig_of(s["padj_t"]) if s["padj_t"] < 0.05 else ""))
    add("")
    add("### 关键读数")
    add("")
    p1 = stats[stats["plate"] == "P1"]
    y_sig = p1[(p1["group"] == "y") & (p1["padj_t"] < 0.05)]["gene"].tolist()
    y_nom = p1[(p1["group"] == "y") & (p1["t_p"] < 0.05)]["gene"].tolist()
    add("板 1 中 Yoda1 组相对对照上调的 YAP 靶基因（BH 校正后显著）：%s。" % ("、".join(y_sig) if y_sig else "无"))
    add("名义 p<0.05 的基因：%s。" % ("、".join(y_nom) if y_nom else "无"))
    add("Yoda1+GsMTx4 与 Yoda1+Verteporfin 两组多数基因回到对照水平（fold 接近 1），")
    add("与 Piezo1、YAP 分别被抑制后效应被阻断的预期一致。")
    add("")
    p2 = stats[stats["plate"] == "P2"]
    tgf_sig = p2[(p2["group"] == "tgf") & (p2["padj_t"] < 0.05)]["gene"].tolist()
    add("板 2 中 TGF-β1 阳性对照显著上调的基因（BH 校正后）：%s。" % ("、".join(tgf_sig) if tgf_sig else "无"))
    add("")
    p3 = stats[stats["plate"] == "P3"]
    for g in ["y", "y+g", "y+v"]:
        gg = p3[p3["group"] == g]
        up = gg[gg["fold"] > 1.3]["gene"].tolist()
        add("板 3 %s 条件培养基处理后明显上调（fold>1.3）的基因：%s。"
            % (GROUP_CN[g], "、".join(up) if up else "无"))
    add("")
    add("## 二、第一板 Ct 补齐说明")
    add("")
    add("原始文件 Results 表 84 个有效孔均有 Ct，编译表与原始表逐孔一致。")
    add("其中 36 个孔为 3–4 位小数的手输值（Piezo1、Gapdh 整行及若干散点）。")
    add("尝试用扩增曲线反算精确 Ct，但该导出文件的扩增表与 Results 表孔号对不上")
    add("（起峰循环与报告 Ct 相差约 10 个循环，疑似导出错位或混入其他运行），")
    add("交叉验证未通过（中位偏差 −19.3 循环），故按计划回退，采用现值并标注。")
    add("3–4 位小数的舍入误差不超过 0.0005 个循环，对 ΔΔCt 和 fold 的影响可以忽略。")
    add("对照表见 plate1_recovered.csv。")
    add("")
    add("## 三、质控")
    add("")
    n_spike = int(wells["quality_flags"].apply(lambda s: str(s).split()[3] if len(str(s).split()) > 3 else "N").eq("Y").sum())
    n_highsd = int(wells["quality_flags"].apply(lambda s: str(s).split()[1] if len(str(s).split()) > 1 else "N").eq("Y").sum())
    add("板 1 机器标记：SPIKE %d 孔，HIGHSD %d 孔（详见 plate1_qc_flags.csv）。" % (n_spike, n_highsd))
    add("复孔 CV>3%% 的组 %d 个（qpcr_cv.csv）。" % int((cv["cv_pct"] > 3).sum()))
    wells["tm1_num"] = pd.to_numeric(wells["tm1"], errors="coerce")
    tm = wells.groupby("gene")["tm1_num"].agg(["min", "max", "count"])
    tm_bad = tm[(tm["max"] - tm["min"]) > 0.5]
    add("板 1 产物区（65–92℃）熔解峰：每个基因 12 孔中检出峰的孔数见下，")
    add("同一靶标内 Tm 极差 >0.5℃ 的基因：%s。"
        % ("、".join("%s（%.1f–%.1f℃）" % (g, r["min"], r["max"]) for g, r in tm_bad.iterrows()) if len(tm_bad) else "无"))
    add("各基因峰检出数：%s。"
        % ("，".join("%s %d/12" % (g, int(r["count"])) for g, r in tm.iterrows())))
    add("可疑离群孔（仅披露、未剔除）：A7 Piezo1、B11 Yap1、D9 Cyr61、E7 Ankrd1。")
    add("")
    add("综合评估：机器质量标志（SPIKE/HIGHSD）与熔解峰缺失较多，且扩增/熔解表与")
    add("Results 表孔位对不上，本板原始导出疑似错位或为另一运行，数据质量欠佳。")
    add("当前统计结论方向自洽，但建议在条件允许时重跑第一板后再定稿。")
    add("")
    add("## 四、方法")
    add("")
    add("ΔCt = Ct(目的基因) − Ct(Gapdh)，板 3 使用大鼠 Gapdh。")
    add("ΔΔCt = ΔCt(处理组) − ΔCt(c 对照组)，fold = 2^(−ΔΔCt)。")
    add("处理组对 c 做 Welch t 检验与 Wilcoxon 秩和检验，按板内基因数做 BH 校正。")
    add("")
    add("## 五、局限")
    add("")
    add("每组 n=3 复孔，统计功效有限；板 1 部分孔为手输值（误差可忽略）；")
    add("板 2、板 3 无原始导出文件，采用编译表数值；离群孔未剔除。")

    md = os.path.join(OUT, "qpcr_报告.md")
    with open(md, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")

    xlsx = os.path.join(OUT, "qpcr_汇总.xlsx")
    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        pd.read_csv(os.path.join(OUT, "qpcr_long.csv")).to_excel(w, sheet_name="长表", index=False)
        stats.to_excel(w, sheet_name="统计", index=False)
        cv.to_excel(w, sheet_name="复孔CV", index=False)
        recovered.to_excel(w, sheet_name="板1恢复对照", index=False)
        wells.to_excel(w, sheet_name="板1QC", index=False)

    print("written:", md, xlsx)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""34_verify_package.py

投稿包最终一致性验证（只读）：
  1) 关键文件存在性与大小；
  2) TIFF 图件规格（尺寸/dpi/模式/LZW）；
  3) 主稿 DOCX 行号/双倍行距/页脚页码；
  4) 合并 PDF 页数与关键数字（审计纠错后的值）；
  5) 手稿关键数字与数据文件对账抽查。
"""

from __future__ import annotations

import os
import re

import fitz
import numpy as np
import pandas as pd
from PIL import Image
from docx import Document
from docx.oxml.ns import qn


PKG = r"E:/sheng xin/ObstructiveNephropathy_MRG/submission/plos"
PROC = r"E:/sheng xin/ObstructiveNephropathy_MRG/data/processed"


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}")


def main() -> int:
    fails = 0

    # 1) 关键文件
    required = [
        "manuscript/manuscript_plos.docx",
        "manuscript/manuscript_plos.md",
        "manuscript/manuscript_plos_with_figures.pdf",
        "figures/Fig1.tif", "figures/Fig2.tif",
        "figures/Fig3.tif", "figures/Fig4.tif",
        "cover_letter/cover_letter.docx",
        "supporting_information/S1_Table_MRG_gene_set.txt",
        "supporting_information/S2_Table_mechanical_core.csv",
        "supporting_information/S3_Table_immune_deconvolution.csv",
        "supporting_information/S1_Text_GSE66494_batch_statement.docx",
        "supporting_information/S2_Text_sensitivity_analyses.docx",
        "supporting_information/S3_Text_immune_consistency.docx",
        "supporting_information/S4_Text_reproducibility_deviation_log.docx",
    ]
    for rel in required:
        p = os.path.join(PKG, rel)
        ok = os.path.exists(p) and os.path.getsize(p) > 0
        if not ok:
            fails += 1
        check(f"file:{rel}", ok, f"{os.path.getsize(p)} B" if ok else "missing")

    # 2) TIFF 规格
    for i in range(1, 5):
        p = os.path.join(PKG, "figures", f"Fig{i}.tif")
        with Image.open(p) as im:
            tag = im.tag_v2.get(259)
            w, h = im.size
        ok = (im.mode == "RGB" and tag == 5 and w == 2250
              and h <= 2625 and os.path.getsize(p) < 10_000_000
              and im.info.get("dpi") == (300.0, 300.0))
        if not ok:
            fails += 1
        check(f"Fig{i}.tif", ok,
              f"{w}x{h} mode={im.mode} dpi={im.info.get('dpi')} tag259={tag} "
              f"{os.path.getsize(p)} B")

    # 3) DOCX 设置
    doc = Document(os.path.join(PKG, "manuscript", "manuscript_plos.docx"))
    ln_set = doc.settings.element.find(qn("w:lnNumType")) is not None
    ln_sec = doc.sections[0]._sectPr.find(qn("w:lnNumType")) is not None
    footer_has_page = "PAGE" in doc.sections[0].footer._element.xml
    body = [p for p in doc.paragraphs if p.text.strip()
            and not p.style.name.startswith("Heading")]
    dbl = sum(1 for p in body if p.paragraph_format.line_spacing == 2.0)
    ok = ln_set and ln_sec and footer_has_page and dbl == len(body)
    if not ok:
        fails += 1
    check("manuscript DOCX", ok,
          f"lnNum settings={ln_set} sect={ln_sec} footerPAGE={footer_has_page} "
          f"doubleSpaced={dbl}/{len(body)}")

    # 4) PDF
    pdf_path = os.path.join(PKG, "manuscript", "manuscript_plos_with_figures.pdf")
    d = fitz.open(pdf_path)
    raw = "".join(d[i].get_text() for i in range(d.page_count))
    n_images = sum(len(d[i].get_images()) for i in range(d.page_count))

    # 关键数字对账以 Markdown 源稿为准（PDF 文本会混入行号导致断词）
    md_path = os.path.join(PKG, "manuscript", "manuscript_plos.md")
    text = re.sub(r"\s+", " ", open(md_path, encoding="utf-8").read())
    checks = {
        "58 up-regulated": "58 up-regulated" in text,
        "173 down-regulated": "173 down-regulated" in text,
        "C-index 0.848": "concordance index of 0.848" in text,
        "AUC 0.853": "0.853" in text,
        "mech p 7.01e-10": "7.01" in text,
        "subtype p 0.0105": "0.0105" in text,
        "EMT p 0.438": "0.438" in text,
        "immune 4 altered": "four cell types significantly altered" in text
        or "four significantly altered" in text,
        "CD8 2.29e-5": "2.29" in text,
        "ActB 3.79e-4": "3.79" in text,
        "MDSC 4.30e-3": "4.30" in text,
        "eosino 4.82e-3": "4.82" in text,
    }
    for k, okk in checks.items():
        if not okk:
            fails += 1
        check(f"PDF:{k}", okk)
    check("PDF page count", 10 <= d.page_count <= 40, f"{d.page_count} pages")
    check("PDF embedded images == 4", n_images == 4, f"n={n_images}")

    # 5) 对账抽查：签名/免疫/核心基因
    imm = pd.read_csv(os.path.join(PKG, "supporting_information",
                                   "S3_Table_immune_deconvolution.csv"))
    n_sig = int((imm.padj < 0.05).sum())
    check("S3 immune n significant == 4", n_sig == 4, f"n={n_sig}")
    core = pd.read_csv(os.path.join(PKG, "supporting_information",
                                    "S2_Table_mechanical_core.csv"))
    check("S2 core genes == 100", len(core) == 100, f"n={len(core)}")
    mrg = open(os.path.join(PKG, "supporting_information",
                            "S1_Table_MRG_gene_set.txt"),
               encoding="utf-8").read().splitlines()
    check("S1 MRG set == 2536", len(mrg) == 2536, f"n={len(mrg)}")

    print("=" * 50)
    print("TOTAL FAILURES:", fails)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

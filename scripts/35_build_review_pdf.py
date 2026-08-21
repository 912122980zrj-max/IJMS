#!/usr/bin/env python
"""35_build_review_pdf.py

生成单文件审阅 PDF：
  1) 文本页来自 manuscript_plos.docx（不含图，先由 Word COM 导出为 text_only.pdf）；
  2) 在文末追加 2 个图版页，用 PyMuPDF 以原始 2250 px TIFF 直接嵌入
     （6.5 in 宽 -> 约 346 dpi，避免 Word 导出时把图压到 1299 px/200 dpi）。

用法：
    python scripts/35_build_review_pdf.py --text text_only.pdf --out manuscript_plos_with_figures.pdf [--force]
"""

from __future__ import annotations

import argparse
import os

import fitz


ROOT = r"E:/sheng xin/ObstructiveNephropathy_MRG"
FIGDIR = os.path.join(ROOT, "submission", "plos", "figures")


def main() -> int:
    ap = argparse.ArgumentParser(description="组装审阅 PDF")
    ap.add_argument("--text", required=True)
    ap.add_argument("--out",
                    default=os.path.join(ROOT, "submission", "plos", "manuscript",
                                         "manuscript_plos_with_figures.pdf"))
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if os.path.exists(args.out) and not args.force:
        raise SystemExit(f"refusing to overwrite: {args.out}")
    doc = fitz.open(args.text)
    page_w, page_h = doc[0].rect.width, doc[0].rect.height
    img_w = 6.5 * 72  # 6.5 in
    left = (page_w - img_w) / 2

    figs = ["Fig1.tif", "Fig2.tif", "Fig3.tif", "Fig4.tif"]
    placements = [[figs[0], figs[1]], [figs[2], figs[3]]]
    for group in placements:
        page = doc.new_page(width=page_w, height=page_h)
        y = 54  # 0.75 in 上边距
        for name in group:
            path = os.path.join(FIGDIR, name)
            with fitz.open(path) as tif:
                pix = tif[0].get_pixmap(dpi=300)
            aspect = pix.height / pix.width
            img_h = img_w * aspect
            # 若超页则等比缩到可用高度
            avail = page_h - y - 54
            if img_h > avail / 2:
                img_h = avail / 2
                w2 = img_h / aspect
                left2 = (page_w - w2) / 2
                rect = fitz.Rect(left2, y, left2 + w2, y + img_h)
            else:
                rect = fitz.Rect(left, y, left + img_w, y + img_h)
            page.insert_text((54, y - 12), name.replace(".tif", ""),
                             fontsize=11, fontname="helv",
                             color=(0, 0, 0))
            page.insert_image(rect, pixmap=pix)
            y += img_h + 36
    doc.save(args.out, garbage=4, deflate=True)
    print(f"written -> {args.out} ({doc.page_count} pages)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

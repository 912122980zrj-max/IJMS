#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""19_qpcr_recover.py —— 第一板原始导出 Ct 恢复（扩增曲线阈值插值反算 + 交叉验证）

输入: C:/Users/91212/Desktop/2026151447.xls（QuantStudio 导出）
输出: results/qpcr_analysis/plate1_recovered.csv
  列: well, well_position, gene, ct_current, ct_recalc, delta,
      candidate_rounded, threshold, tm1, quality_flags
逻辑:
  1. Results 表定位 "Well Position" 表头，取 A1-G12 共 84 孔。
  2. Amplification Data 表定位 "Well/Cycle/Delta Rn" 表头，重建逐孔曲线。
  3. 每孔用 Ct Threshold（缺省 0.1）在 Delta Rn 曲线上线性插值反算 Ct。
  4. 用小数位 <=4 识别疑似被覆盖（手输）的孔，输出对照表。
  5. 用小数位 >4 的孔做交叉验证，报告偏差中位数/95 分位。
"""

import csv
import os
import sys
from collections import defaultdict

import xlrd

RAW = r"C:/Users/91212/Desktop/2026151447.xls"
OUT = r"E:/sheng xin/ObstructiveNephropathy_MRG/results/qpcr_analysis"
OUT_CSV = os.path.join(OUT, "plate1_recovered.csv")
OUT_WELLS = os.path.join(OUT, "plate1_wells.csv")

FLAG_COLS = ["BADROX", "HIGHSD", "NOISE", "SPIKE", "EXPFAIL", "PRFDROP", "MTP"]


def header_row(sheet, *needles):
    for r in range(sheet.nrows):
        cells = [str(sheet.cell_value(r, c)).strip() for c in range(sheet.ncols)]
        if all(n in cells for n in needles):
            return r, {name: i for i, name in enumerate(cells)}
    return None, {}


def fraction_digits(x):
    s = repr(float(x))
    if "e" in s or "E" in s:
        return 99
    if "." not in s:
        return 0
    return len(s.split(".")[1])


def main():
    wb = xlrd.open_workbook(RAW)

    res = wb.sheet_by_name("Results")
    hdr, col = header_row(res, "Well Position")
    if hdr is None:
        sys.exit("Results 表未找到 Well Position 表头")
    wells = []
    for r in range(hdr + 1, res.nrows):
        well = res.cell_value(r, col["Well"])
        if not isinstance(well, float) or well > 96:
            continue
        pos = str(res.cell_value(r, col["Well Position"])).strip()
        gene = str(res.cell_value(r, col["Target Name"])).strip()
        ct = res.cell_value(r, col["CT"])
        thr = res.cell_value(r, col["Ct Threshold"])
        tm = res.cell_value(r, col["Tm1"])
        flags = " ".join(
            str(res.cell_value(r, col[f])).strip()
            for f in FLAG_COLS
            if f in col and str(res.cell_value(r, col[f])).strip()
        )
        wells.append(dict(well=int(well), pos=pos, gene=gene, ct=ct,
                          thr=thr, tm=tm, flags=flags))

    amp = wb.sheet_by_name("Amplification Data")
    ahdr, acol = header_row(amp, "Cycle", "Delta Rn")
    if ahdr is None:
        sys.exit("Amplification Data 未找到曲线表头")
    curves = defaultdict(list)
    for r in range(ahdr + 1, amp.nrows):
        w = amp.cell_value(r, acol["Well"])
        cyc = amp.cell_value(r, acol["Cycle"])
        drn = amp.cell_value(r, acol["Delta Rn"])
        if isinstance(w, float) and isinstance(cyc, float):
            curves[int(w)].append((cyc, drn))

    melt = wb.sheet_by_name("Melt Curve Raw Data")
    mhdr, mcol = header_row(melt, "Temperature", "Derivative")
    tm_by_well = {}
    if mhdr is not None:
        for r in range(mhdr + 1, melt.nrows):
            w = melt.cell_value(r, mcol["Well"])
            temp = melt.cell_value(r, mcol["Temperature"])
            deriv = melt.cell_value(r, mcol["Derivative"])
            if isinstance(w, float) and isinstance(temp, float) and isinstance(deriv, float):
                key = int(w)
                if key not in tm_by_well:
                    tm_by_well[key] = []
                tm_by_well[key].append((temp, deriv))

    def recompute(curve, thr):
        prev = None
        for cyc, drn in curve:
            if not isinstance(drn, (int, float)) or not isinstance(cyc, (int, float)):
                continue
            if drn >= thr and prev is not None:
                c0, d0 = prev
                c1, d1 = cyc, drn
                if d1 == d0:
                    return c1
                return c0 + (thr - d0) / (d1 - d0) * (c1 - c0)
            prev = (cyc, drn)
        return None

    default_thr = 0.1
    rows = []
    for w in wells:
        thr = w["thr"] if isinstance(w["thr"], (int, float)) and w["thr"] > 0 else default_thr
        curve = curves.get(w["well"], [])
        ct_calc = recompute(curve, thr) if curve else None
        ct_cur = w["ct"] if isinstance(w["ct"], (int, float)) else None
        delta = None
        if ct_cur is not None and ct_calc is not None:
            delta = ct_calc - ct_cur
        cand = 0
        if ct_cur is not None and fraction_digits(ct_cur) <= 4:
            cand = 1
        rows.append(dict(
            well=w["well"], pos=w["pos"], gene=w["gene"],
            ct_current=round(ct_cur, 6) if ct_cur is not None else "",
            ct_recalc=round(ct_calc, 6) if ct_calc is not None else "",
            delta=round(delta, 6) if delta is not None else "",
            candidate_rounded=cand,
            threshold=round(thr, 4),
            tm1=round(w["tm"], 3) if isinstance(w["tm"], (int, float)) else "",
            quality_flags=w["flags"],
        ))

    # 交叉验证：小数位 >4 的孔（未受损）
    # Tm：产物区（65-92℃）导数最大值处温度；全区间最大值另存
    for r in rows:
        pts = tm_by_well.get(r["well"], [])
        if pts:
            prod = [(t, d) for t, d in pts if 65 <= t <= 92]
            r["tm1"] = round(max(prod, key=lambda p: p[1])[0], 3) if prod else ""
            r["tm_full"] = round(max(pts, key=lambda p: p[1])[0], 3)
        else:
            r.setdefault("tm_full", "")

    val = [r["delta"] for r in rows
           if r["candidate_rounded"] == 0 and r["delta"] != ""]
    if val:
        val_sorted = sorted(val)
        med = val_sorted[len(val_sorted) // 2]
        p95 = val_sorted[int(len(val_sorted) * 0.95)] if val_sorted else None
        p_lt_01 = sum(1 for d in val if abs(d) < 0.1) / len(val)
    else:
        med = p95 = p_lt_01 = None
    recovery_ok = (med is not None and abs(med) < 0.05 and
                   (p_lt_01 is None or p_lt_01 >= 0.95))

    with open(OUT_WELLS, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "well", "pos", "gene", "ct", "tm1", "tm_full", "quality_flags",
            "candidate_rounded"])
        writer.writeheader()
        writer.writerows([
            dict(well=r["well"], pos=r["pos"], gene=r["gene"],
                 ct=r["ct_current"], tm1=r["tm1"], tm_full=r["tm_full"],
                 quality_flags=r["quality_flags"],
                 candidate_rounded=r["candidate_rounded"])
            for r in rows
        ])

    os.makedirs(OUT, exist_ok=True)
    with open(OUT_CSV, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print("wells:", len(rows))
    print("candidate_rounded:", sum(r["candidate_rounded"] for r in rows))
    print("validation n:", len(val),
          "median_delta:", round(med, 4) if med is not None else "NA",
          "p95_abs:", round(max(abs(d) for d in val), 4) if val else "NA",
          "frac_abs<0.1:", round(p_lt_01, 3) if p_lt_01 is not None else "NA")
    print("recovery_status:", "VALIDATED" if recovery_ok else "FAILED (fallback to current values)")
    print("written:", OUT_CSV, OUT_WELLS)


if __name__ == "__main__":
    main()

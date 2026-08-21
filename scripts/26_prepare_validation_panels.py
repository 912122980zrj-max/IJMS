#!/usr/bin/env python
"""26_prepare_validation_panels.py

为组图 4C/4D 制备底层数据（只读原始数据，结果写 data/processed/panel_export/）：
  spatial    —— 流式两遍处理 CosMx 表达矩阵（与 11_spatial_validate.R 的 z() 口径一致），
                输出 FOV 级 MRG/ECM 聚合与 20,000 个细胞的随机子样本；
  timecourse —— 读取 GSE118339 15 个 TPM 文件，Ensembl REST 映射 ID→symbol，
                按核心基因小鼠同源集计算逐样本机械核心评分。

用法:
    python scripts/26_prepare_validation_panels.py spatial
    python scripts/26_prepare_validation_panels.py timecourse
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
import time

import numpy as np
import pandas as pd
import requests


ROOT = r"E:/sheng xin/ObstructiveNephropathy_MRG"
OUT = os.path.join(ROOT, "data", "processed", "panel_export")


def log(msg: str) -> None:
    print(f"[prepare-validation] {msg}", flush=True)


def load_core_genes() -> set[str]:
    w = pd.read_csv(os.path.join(ROOT, "data", "processed", "wgcna_core_gs_mm.csv"))
    return set(w.loc[w.core, "gene"])


def spatial() -> None:
    core = load_core_genes()
    ecm = {"COL1A1", "COL3A1", "DCN", "FN1", "ACTA2", "COL4A1", "COL6A1", "VCAN"}
    path = os.path.join(ROOT, "data", "raw", "GSE282059", "suppl",
                        "GSE282059_2B_FFPE_exprMat_file.csv.gz")
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split(",")
    cols = ["fov", "cell_ID"] + [c for c in header[2:] if c in core or c in ecm]
    gene_cols = cols[2:]
    mrg_cols = [c for c in gene_cols if c in core]
    ecm_cols = [c for c in gene_cols if c in ecm]
    log(f"panel genes: MRG={len(mrg_cols)} ECM={len(ecm_cols)}")

    # Pass 1: mean / sd of log1p per gene column
    n = np.zeros(len(gene_cols))
    sums = np.zeros(len(gene_cols))
    sumsq = np.zeros(len(gene_cols))
    chunks = pd.read_csv(path, usecols=cols, chunksize=200_000)
    for chunk in chunks:
        m = np.log1p(chunk[gene_cols].to_numpy(dtype=np.float64))
        n += (~np.isnan(m)).sum(axis=0)
        m = np.nan_to_num(m)
        sums += m.sum(axis=0)
        sumsq += (m * m).sum(axis=0)
    mean = sums / np.maximum(n, 1)
    # 样本方差（ddof=1），与 R 的 scale()/sd() 口径一致
    var = (sumsq - n * mean * mean) / np.maximum(n - 1, 1)
    sd = np.sqrt(np.clip(var, 0, None))
    keep = sd > 0
    log(f"columns kept (sd>0): {int(keep.sum())}/{len(gene_cols)}")

    # Pass 2: per-cell scores + FOV aggregation + fixed random sample
    rng = np.random.default_rng(42)
    fov_agg: dict[str, list[float]] = {}
    sample_rows: list[tuple[float, float]] = []
    chunks = pd.read_csv(path, usecols=cols, chunksize=200_000)
    for chunk in chunks:
        m = np.log1p(chunk[gene_cols].to_numpy(dtype=np.float64))
        m = np.nan_to_num(m)
        z = (m - mean) / sd
        mrg_idx = np.array([gene_cols.index(c) for c in mrg_cols])
        ecm_idx = np.array([gene_cols.index(c) for c in ecm_cols])
        mrg = z[:, mrg_idx[keep[mrg_idx]]].mean(axis=1)
        ecmv = z[:, ecm_idx[keep[ecm_idx]]].mean(axis=1)
        fovs = chunk["fov"].to_numpy()
        for fov, a, b in zip(fovs, mrg, ecmv):
            fov_agg.setdefault(fov, []).append((float(a), float(b)))
        if len(sample_rows) < 20_000:
            sample_rows.extend(zip(mrg.tolist(), ecmv.tolist()))
    fov_df = pd.DataFrame(
        [(f, np.mean([r[0] for r in v]), np.mean([r[1] for r in v]), len(v))
         for f, v in fov_agg.items()],
        columns=["fov", "MRG", "ECM", "n"],
    )
    fov_df.to_csv(os.path.join(OUT, "spatial_fov.csv"), index=False)
    if len(sample_rows) > 20_000:
        sample_rows = rng.choice(np.array(sample_rows), size=20_000, replace=False)
    pd.DataFrame(sample_rows, columns=["MRG", "ECM"]).to_csv(
        os.path.join(OUT, "spatial_cells_sample.csv"), index=False)
    log(f"spatial done: {len(fov_df)} FOVs, cells={fov_df.n.sum()}")


def timecourse() -> None:
    indir = os.path.join(ROOT, "data", "raw", "GSE118339")
    files = sorted(f for f in os.listdir(indir) if f.endswith(".tpm.txt.gz"))
    frames = {}
    for f in files:
        sample = re.match(r"GSM\d+_(.+?)\.tpm\.txt\.gz$", f).group(1)
        d = pd.read_csv(os.path.join(indir, f), usecols=[0, 1])
        frames[sample] = d.set_index(d.columns[0]).iloc[:, 0]
    M = pd.DataFrame(frames)
    log(f"timecourse matrix: {M.shape[0]} genes x {M.shape[1]} samples")

    ids = list(M.index)
    mapping: dict[str, str] = {}
    base = "https://rest.ensembl.org/lookup/id"
    for i in range(0, len(ids), 1000):
        batch = ids[i:i + 1000]
        for attempt in range(4):
            try:
                r = requests.post(base, headers={"Content-Type": "application/json"},
                                  json={"ids": batch}, timeout=60)
                if r.status_code == 200:
                    for k, v in r.json().items():
                        if v and "display_name" in v:
                            mapping[k] = v["display_name"]
                    break
            except requests.RequestException:
                time.sleep(2 * (attempt + 1))
        time.sleep(0.2)
    log(f"ensembl mapping: {len(mapping)}/{len(ids)} ids")

    orth = pd.read_csv(os.path.join(ROOT, "data", "processed", "ortholog_map.csv"))
    orth = orth[orth["mmusculus_homolog_orthology_type"] == "ortholog_one2one"]
    core = load_core_genes()
    mouse_core = set(orth.loc[orth.hgnc_symbol.isin(core),
                              "mmusculus_homolog_associated_gene_name"].dropna())
    keep = M.index.isin([k for k, v in mapping.items() if v in mouse_core])
    M = M.loc[keep]
    log(f"core mouse genes present: {M.shape[0]}")

    v = np.log2(M.to_numpy(dtype=np.float64) + 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        zmat = (v - v.mean(axis=1, keepdims=True)) / v.std(axis=1, ddof=1, keepdims=True)
    score = np.nanmean(zmat, axis=0)
    day_map = {s: 0 for s in M.columns if s.startswith("normal")}
    for k, d in [("day3", 3), ("day7", 7), ("day14", 14)]:
        day_map.update({s: d for s in M.columns if s.startswith(k)})
    out = pd.DataFrame({"sample": M.columns,
                        "day": [day_map[s] for s in M.columns],
                        "score": score})
    out.to_csv(os.path.join(OUT, "timecourse_samples.csv"), index=False)
    log(f"timecourse done: {len(out)} samples")


def main() -> int:
    ap = argparse.ArgumentParser(description="制备空间/时间序列面板底层数据")
    ap.add_argument("task", choices=["spatial", "timecourse"])
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    if args.task == "spatial":
        spatial()
    else:
        timecourse()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

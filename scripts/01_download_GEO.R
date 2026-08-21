#!/usr/bin/env Rscript
# 01_download_GEO.R —— 选题三数据下载
# 推荐改用 tools/download_geo.py（断点续传/重试/幂等，双击 download_all.bat 即可）；
# 本脚本为 GEOquery 备选实现，保留用于程序化下游读取。
# 数据清单（NCBI GEO esummary 核验，2026-08-17）：
#   人主队列：GSE66494（CKD 肾小管间质，n=61，GPL6480 Agilent 4x44K）
#   小鼠 UUO bulk：GSE299417（sham+UUO RNA-seq，n=24，GPL24247）
#   小鼠 UUO/纤维化 scRNA：GSE175412（n=3）、GSE198621（n=11）、GSE269062/269063（TRPC6，scRNA+空间）
#   人肾空间：GSE282059（CosMx 6,000-plex）
#   可选参考：GSE241634（CD248 机械转导 UUO，n=6）

suppressPackageStartupMessages({
  library(optparse)
  library(GEOquery)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "01_download_GEO.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--acc"), type = "character",
              default = "GSE66494,GSE299417,GSE175412,GSE198621,GSE269062,GSE282059",
              help = "逗号分隔 GSE 编号"),
  make_option(c("--outdir"), type = "character", default = DATA_RAW, help = "原始数据目录"),
  make_option(c("--gsm"), action = "store_true", default = FALSE,
              help = "同时下载补充文件（counts/scRNA 矩阵）")
)
opt <- parse_args(OptionParser(option_list = option_list))
accs <- trimws(strsplit(opt$acc, ",")[[1]])

main <- function() {
  ensure_dirs()
  for (a in accs) {
    log_msg("download", a)
    dir.create(file.path(opt$outdir, a), recursive = TRUE, showWarnings = FALSE)
    gse <- getGEO(a, GSEMatrix = TRUE, destdir = file.path(opt$outdir, a), getGPL = TRUE)
    log_msg("download", paste0(a, " OK: ", length(gse), " platform(s)"))
    if (opt$gsm) {
      getGEOSuppFiles(a, makeDirectory = FALSE, baseDir = file.path(opt$outdir, a))
      log_msg("download", paste0(a, " supplementary files OK"))
    }
  }
  log_msg("download", "DONE")
}

if (!interactive()) main()

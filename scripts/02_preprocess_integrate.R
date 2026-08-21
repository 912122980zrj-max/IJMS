#!/usr/bin/env Rscript
# 02_preprocess_integrate.R —— 三个队列预处理（真实实现）
# 数据现实（见 results/deviation_log.md D1/D2）：
#   人主队列 GSE115857（IgAN 病变分级 G1/G2/G3 轴）：Illumina HT-12，已 log2，探针 ILMN_ → illuminaHumanv4.db
#   人验证队列 GSE66494（CKD 53 / normal 8，discovery/validation）：Agilent GPL6480，线性强度 → log2 + 分位数归一
#   小鼠 GSE299417（6 组 × 4：Sham/UUO × Veh/早期/晚期干预）：gene_exp.xlsx FPKM → log2(FPKM+1)
# 输出：data/processed/{expr_human_main,expr_human_ckd,expr_mouse}.rds + 对应 meta CSV

suppressPackageStartupMessages({
  library(optparse)
  library(GEOquery)
  library(limma)
  library(readxl)
  library(illuminaHumanv4.db)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "02_preprocess_integrate.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--outdir"), type = "character", default = DATA_PROC, help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

collapse_to_genes <- function(expr, symbols) {
  keep <- !is.na(symbols) & symbols != "" & symbols != "---"
  expr <- expr[keep, , drop = FALSE]
  sym <- symbols[keep]
  ord <- order(sym, -rowMeans(expr))
  first <- !duplicated(sym[ord])
  res <- expr[ord[first], , drop = FALSE]
  rownames(res) <- sym[ord[first]]
  res
}

main <- function() {
  ensure_dirs()

  log_msg("human_main", "GSE115857（IgAN 病变分级）")
  g1 <- getGEO(filename = file.path(DATA_RAW, "GSE115857", "matrix", "GSE115857_series_matrix.txt.gz"),
               getGPL = FALSE)
  e1 <- exprs(g1)
  e1 <- normalizeBetweenArrays(e1, method = "quantile")   # D12：GSE115857 全基因组偏移，分位数归一
  sym1 <- unlist(mget(rownames(e1), illuminaHumanv4SYMBOL, ifnotfound = NA))
  expr_main <- collapse_to_genes(e1, sym1)
  p1 <- pData(g1)
  meta_main <- data.frame(
    sample = p1$geo_accession,
    status = gsub("subject status: ", "", p1[["characteristics_ch1"]]),
    grade = gsub("classification: ", "", p1[["characteristics_ch1.1"]]),
    stringsAsFactors = FALSE
  )
  meta_main$group <- ifelse(meta_main$status == "IgAN patient" &
                              meta_main$grade %in% c("G1", "G2", "G3"),
                            paste0("IgAN_", meta_main$grade), "other")
  safe_write(file.path(opt$outdir, "expr_human_main.rds"), function(p) saveRDS(expr_main, p))
  safe_write(file.path(opt$outdir, "meta_human_main.csv"), function(p) write_csv_utf8(meta_main, p))
  log_msg("human_main", paste0(dim(expr_main)[1], " genes x ", dim(expr_main)[2], " samples"))

  log_msg("human_ckd", "GSE66494（CKD vs normal）")
  g2 <- getGEO(filename = file.path(DATA_RAW, "GSE66494", "matrix", "GSE66494_series_matrix.txt.gz"),
               getGPL = FALSE)
  e2 <- log2(exprs(g2))
  e2 <- normalizeBetweenArrays(e2, method = "quantile")
  annot <- read.delim(gzfile(file.path(DATA_RAW, "GSE66494", "GPL6480.annot.gz")),
                      skip = 28, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  sym2 <- annot[["Gene.symbol"]][match(rownames(e2), annot[["ID"]])]
  expr_ckd <- collapse_to_genes(e2, sym2)
  p2 <- pData(g2)
  meta_ckd <- data.frame(
    sample = p2$geo_accession,
    group = ifelse(grepl("normal", p2[["disease status:ch1"]]), "normal", "CKD"),
    study_set = p2[["study set:ch1"]],
    stringsAsFactors = FALSE
  )
  safe_write(file.path(opt$outdir, "expr_human_ckd.rds"), function(p) saveRDS(expr_ckd, p))
  safe_write(file.path(opt$outdir, "meta_human_ckd.csv"), function(p) write_csv_utf8(meta_ckd, p))
  log_msg("human_ckd", paste0(dim(expr_ckd)[1], " genes x ", dim(expr_ckd)[2], " samples"))

  log_msg("mouse", "GSE299417（gene FPKM → log2）")
  x <- read_excel(file.path(DATA_RAW, "GSE299417", "suppl", "GSE299417_Expression_Gene.xlsx"),
                  sheet = "gene_exp", skip = 8)
  x <- as.data.frame(x)
  samp_cols <- grep("^Sample[0-9]+-", colnames(x), value = TRUE)
  m <- as.matrix(x[, samp_cols])
  storage.mode(m) <- "double"
  gname <- trimws(x$Gene_Name)
  dup <- duplicated(gname) | gname == "" | is.na(gname)
  m <- m[!dup, , drop = FALSE]
  gname <- gname[!dup]
  rownames(m) <- gname
  expr_mouse <- log2(m + 1)
  g3 <- getGEO(filename = file.path(DATA_RAW, "GSE299417", "matrix", "GSE299417_series_matrix.txt.gz"),
               getGPL = FALSE)
  title3 <- as.data.frame(g3)$title
  key3 <- paste0("Sample", sub("^\\s*Sample (\\d+)-\\d+.*$", "\\1", title3))
  group_map <- c("Sample1" = "Sham_Veh", "Sample2" = "UUO_Veh",
                 "Sample3" = "Sham_Early", "Sample4" = "UUO_Early",
                 "Sample5" = "Sham_Late", "Sample6" = "UUO_Late")
  meta_mouse <- data.frame(sample = as.data.frame(g3)$geo_accession,
                           key = key3,
                           group = unname(group_map[key3]),
                           stringsAsFactors = FALSE)
  colnames(expr_mouse) <- meta_mouse$sample[match(colnames(expr_mouse), meta_mouse$key)]
  safe_write(file.path(opt$outdir, "expr_mouse.rds"), function(p) saveRDS(expr_mouse, p))
  safe_write(file.path(opt$outdir, "meta_mouse.csv"), function(p) write_csv_utf8(meta_mouse, p))
  log_msg("mouse", paste0(dim(expr_mouse)[1], " genes x ", dim(expr_mouse)[2], " samples"))

  log_msg("preprocess", "DONE")
}

if (!interactive()) main()

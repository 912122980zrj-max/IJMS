#!/usr/bin/env Rscript
# 04_HVG_ssGSEA.R —— 机械基因 HVG + ssGSEA 打分（真实实现）
# 输入：GSE115857（IgAN 队列）表达 + IgAN_G3 vs G1 的 DEG；
# 步骤：DEG ∩ MRG 分三集合 → 各取高变基因（MS/MM/CF 默认 25/46/31，沿用 IPF 口径）→
#       按 logFC 符号拆上/下调 → GSVA ssGSEA 打分 → G3 vs G1 箱线比较（wilcoxon）。

suppressPackageStartupMessages({
  library(optparse)
  library(GSVA)
  library(ggpubr)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "04_HVG_ssGSEA.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--expr"), type = "character",
              default = file.path(DATA_PROC, "expr_human_main.rds"), help = "主队列表达矩阵（GSE115857）"),
  make_option(c("--meta"), type = "character",
              default = file.path(DATA_PROC, "meta_human_main.csv"), help = "主队列表型"),
  make_option(c("--deg"), type = "character",
              default = file.path(DATA_PROC, "deg_igag_vs_control.csv"), help = "03 输出 IgAN vs 对照 DEG 表"),
  make_option(c("--top-ms"), dest = "top_ms", type = "integer", default = 25, help = "MS HVG 数"),
  make_option(c("--top-mm"), dest = "top_mm", type = "integer", default = 46, help = "MM HVG 数"),
  make_option(c("--top-cf"), dest = "top_cf", type = "integer", default = 31, help = "CF HVG 数"),
  make_option(c("--outdir"), type = "character", default = DATA_PROC, help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

select_hvg <- function(expr, genes, top_n) {
  g <- intersect(genes, rownames(expr))
  if (length(g) == 0L) return(character(0))
  v <- apply(expr[g, , drop = FALSE], 1, stats::var)
  names(sort(v, decreasing = TRUE))[seq_len(min(top_n, length(v)))]
}

main <- function() {
  ensure_dirs()
  expr <- readRDS(opt$expr)
  meta <- read.csv(opt$meta, stringsAsFactors = FALSE)
  tt <- read.csv(opt$deg, stringsAsFactors = FALSE, check.names = FALSE)
  lfc <- setNames(tt$logFC, tt$gene)
  sig_genes <- tt$gene[abs(tt$logFC) > 0.585 & tt$adj.P.Val < 0.05]
  cat_mrg <- list(
    MS = readLines(file.path(DATA_REF, "MRG_MS_union.txt")),
    MM = readLines(file.path(DATA_REF, "MRG_MM_union.txt")),
    CF = readLines(file.path(DATA_REF, "MRG_CF_union.txt"))
  )

  hvg <- list(
    MS = select_hvg(expr, intersect(sig_genes, cat_mrg$MS), opt$top_ms),
    MM = select_hvg(expr, intersect(sig_genes, cat_mrg$MM), opt$top_mm),
    CF = select_hvg(expr, intersect(sig_genes, cat_mrg$CF), opt$top_cf)
  )
  safe_write(file.path(opt$outdir, "hvg_by_category.rds"), function(p) saveRDS(hvg, p))
  log_msg("hvg", paste0("MS/MM/CF = ", paste(vapply(hvg, length, integer(1)), collapse = "/")))

  # ssGSEA 基因集：用全部 MRG 交叠 DEG 按 logFC 方向拆分（避免 HVG 全同向导致空集）
  mrg_all <- readLines(file.path(DATA_REF, "MRG_union.txt"))
  overlap <- intersect(sig_genes, mrg_all)
  sets <- list(
    MRG_up = overlap[!is.na(lfc[overlap]) & lfc[overlap] > 0],
    MRG_down = overlap[!is.na(lfc[overlap]) & lfc[overlap] < 0]
  )
  log_msg("ssgsea", paste0("gene sets: up=", length(sets$MRG_up), " down=", length(sets$MRG_down)))
  score <- gsva(ssgseaParam(exprData = as.matrix(expr), geneSets = sets))
  safe_write(file.path(opt$outdir, "ssgsea_mech_scores.csv"), function(p) {
    out <- as.data.frame(t(score))
    out <- data.frame(sample = rownames(out), out, check.names = FALSE)
    write_csv_utf8(out, p)
  })

  meta$disease <- ifelse(meta$status == "IgAN patient", "IgAN", "control")
  df <- data.frame(group = factor(meta$disease, levels = c("control", "IgAN")),
                   score = as.numeric(score["MRG_up", ]))
  p <- ggboxplot(df, x = "group", y = "score", color = "group",
                 palette = CB_PALETTE, add = "jitter", xlab = NULL,
                 ylab = "ssGSEA score (MRG up)") +
    stat_compare_means(method = "wilcox.test") +
    pub_theme()
  save_pub_fig(p, file.path(FIG_DIR, "fig_mrg_score_igag"))
  log_msg("ssgsea", paste0("wilcox p = ", signif(wilcox.test(score ~ group, data = df)$p.value, 3)))
  log_msg("ssgsea", "DONE")
}

if (!interactive()) main()

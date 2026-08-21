#!/usr/bin/env Rscript
# 28_30ctrl_recompute.R
#
# 投稿前口径修正：主手稿明确"非进展对照 = MCD 12 + MN 11 + 供肾 7（共 30），FSGS 单例不进入对比"，
# 但脚本 04/09/16 的 disease 定义把 FSGS 一并归入对照（31 例），导致正文中的机械评分 p、
# Hallmark p 与 28 免疫细胞 p 值基于 31 例对照。本脚本按 30 例对照口径用 R wilcox.test
# （精确检验，含连续性校正）重算这些统计量，并把免疫结果另存为新文件（不覆盖旧结果，
# 保留溯源）。所有重算只读原始/中间数据，输出到 data/processed/panel_export/。

suppressPackageStartupMessages({
  library(optparse)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "28_30ctrl_recompute.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--outdir"), type = "character",
              default = file.path(DATA_PROC, "panel_export"),
              help = "输出目录（新增文件，不覆盖）")
)
opt <- parse_args(OptionParser(option_list = option_list))

main <- function() {
  ensure_dirs()
  if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

  meta <- read.csv(file.path(DATA_PROC, "meta_human_main.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
  # 与手稿一致：对照 = MCD + MN + 供肾（30），FSGS 排除
  keep <- meta$status %in% c("IgAN patient", "minimal change disease",
                             "Membranous glomerulonephritis", "Living donor")
  meta <- meta[keep, , drop = FALSE]
  meta$disease <- ifelse(meta$status == "IgAN patient", "IgAN", "control")

  score <- read.csv(file.path(DATA_PROC, "ssgsea_mech_scores.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  hall <- read.csv(file.path(DATA_PROC, "hallmark_scores.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
  imm <- read.csv(file.path(DATA_PROC, "immune28_charoentong_scores.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)

  m <- match(meta$sample, score$sample)
  score <- score[m, , drop = FALSE]
  h <- match(meta$sample, hall$sample)
  hall <- hall[h, , drop = FALSE]
  i <- match(meta$sample, imm$sample)
  imm <- imm[i, , drop = FALSE]

  wt <- function(x) {
    wilcox.test(x[meta$disease == "IgAN"], x[meta$disease == "control"])$p.value
  }

  summ <- data.frame(
    metric = c("mech_score_igan_vs_30ctrl_wilcox_p",
               "HALLMARK_TGF_BETA_SIGNALING_p",
               "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_p",
               "n_igan", "n_control"),
    value = c(wt(score$MRG_up),
              wt(as.numeric(hall[["HALLMARK_TGF_BETA_SIGNALING"]])),
              wt(as.numeric(hall[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]])),
              sum(meta$disease == "IgAN"),
              sum(meta$disease == "control"))
  )
  out1 <- file.path(opt$outdir, "30ctrl_core_stats.csv")
  if (file.exists(out1)) stop("Refusing to overwrite: ", out1)
  write.csv(summ, out1, row.names = FALSE, fileEncoding = "UTF-8")
  log_msg("30ctrl", paste0("core stats written -> ", out1))
  log_msg("30ctrl", paste0("mech p = ", signif(summ$value[1], 4),
                           " | TGFB p = ", signif(summ$value[2], 4),
                           " | EMT p = ", signif(summ$value[3], 4)))

  cells28 <- setdiff(colnames(imm), "sample")
  res <- do.call(rbind, lapply(cells28, function(nm) {
    x <- as.numeric(imm[[nm]])
    data.frame(celltype = nm,
               IgAN = mean(x[meta$disease == "IgAN"]),
               control = mean(x[meta$disease == "control"]),
               p = wt(x))
  }))
  res$padj <- p.adjust(res$p, method = "BH")
  res <- res[order(res$p), ]
  out2 <- file.path(opt$outdir, "immune28_30ctrl_results.csv")
  if (file.exists(out2)) stop("Refusing to overwrite: ", out2)
  write.csv(res, out2, row.names = FALSE, fileEncoding = "UTF-8")
  log_msg("30ctrl", paste0("immune28 significant (BH<0.05) n = ",
                           sum(res$padj < 0.05)))
  log_msg("30ctrl", "DONE")
}

if (!interactive()) main()

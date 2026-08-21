#!/usr/bin/env Rscript
# 20_qpcr_analysis.R —— 三板 qPCR 质控 + 2^-ΔΔCt + 组间统计 + 报告
# 输入: Desktop/qpcr.xlsx（三板 Ct）、results/qpcr_analysis/plate1_wells.csv（板1 QC）
# 分组: 板1 每基因 12 孔 = c(1-3)/y(4-6)/y+g(7-9)/y+v(10-12)
#       板2 = c(1-3)/tgf(4-6)；板3 同板1

suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "20_qpcr_analysis.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--xlsx"), type = "character",
              default = "C:/Users/91212/Desktop/qpcr.xlsx", help = "三板 Ct 编译表"),
  make_option(c("--wells"), type = "character",
              default = file.path(RES_DIR, "qpcr_analysis", "plate1_wells.csv"),
              help = "板1 原始孔位表"),
  make_option(c("--outdir"), type = "character",
              default = file.path(RES_DIR, "qpcr_analysis"), help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

main <- function() {
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  x <- suppressWarnings(as.data.frame(read_excel(opt$xlsx, col_names = FALSE)))

  parse_block <- function(gene_rows, ncol_use, groups) {
    out <- list()
    for (i in seq_along(gene_rows)) {
      vals <- as.numeric(x[gene_rows[i], 2:(ncol_use + 1)])
      gname <- trimws(as.character(x[gene_rows[i], 1]))
      for (j in seq_len(ncol_use)) {
        out[[length(out) + 1]] <- data.frame(
          gene = gname, col = j, group = groups[j],
          rep = ((j - 1) %% 3) + 1, ct = vals[j], stringsAsFactors = FALSE)
      }
    }
    do.call(rbind, out)
  }

  p1_groups <- rep(c("c", "c", "c", "y", "y", "y", "y+g", "y+g", "y+g",
                     "y+v", "y+v", "y+v"), 1)
  p3_groups <- p1_groups
  p2_groups <- c("c", "c", "c", "tgf", "tgf", "tgf")

  p1 <- parse_block(3:9, 12, p1_groups); p1$plate <- "P1"
  p2 <- parse_block(13:18, 6, p2_groups); p2$plate <- "P2"
  p3 <- parse_block(23:26, 12, p3_groups); p3$plate <- "P3"
  long <- rbind(p1, p2, p3)

  # 板1 原始 QC 关联（孔位 = 基因行 A-G × 列 1-12）
  qc1 <- read.csv(opt$wells, stringsAsFactors = FALSE)
  row_letter <- c(Piezo1 = "A", Yap1 = "B", Ctgf = "C", Cyr61 = "D",
                  Ankrd1 = "E", Tgfb1 = "F", Gapdh = "G")
  qc1$gene <- names(row_letter)[match(substr(qc1$pos, 1, 1),
                                      unname(row_letter))]
  qc1$col <- as.integer(substr(qc1$pos, 2, nchar(qc1$pos)))
  # 用位序解析 SPIKE/HIGHSD（quality_flags 为 7 列 N/Y 拼接，第2列 HIGHSD，第4列 SPIKE）
  long$flag_spike <- NA; long$flag_highsd <- NA; long$tm1 <- NA_real_
  for (i in seq_len(nrow(long))) {
    if (long$plate[i] != "P1") next
    hit <- qc1[qc1$gene == long$gene[i] & qc1$col == long$col[i], ]
    if (nrow(hit) == 1) {
      fl <- unlist(strsplit(hit$quality_flags[1], " "))[1:7]
      fl[is.na(fl)] <- "N"
      long$flag_highsd[i] <- as.integer(fl[2] == "Y")
      long$flag_spike[i] <- as.integer(fl[4] == "Y")
      long$tm1[i] <- hit$tm1[1]
    }
  }

  # Gapdh 内参
  gap <- do.call(rbind, lapply(split(long, long$plate), function(d) {
    g <- d[grepl("^Gapdh", d$gene), ]
    data.frame(plate = d$plate[1], col = g$col, gap_ct = g$ct)
  }))
  long$gap_ct <- gap$gap_ct[match(paste(long$plate, long$col),
                                  paste(gap$plate, gap$col))]
  long <- long[!grepl("Gapdh", long$gene), ]
  long$dct <- long$ct - long$gap_ct

  # 组统计
  stats <- do.call(rbind, lapply(split(long, paste(long$plate, long$gene)), function(d) {
    ctrl <- mean(d$dct[d$group == "c"])
    do.call(rbind, lapply(setdiff(unique(d$group), "c"), function(g) {
      v <- d$dct[d$group == g]
      dd <- mean(v) - ctrl
      tt <- t.test(d$dct[d$group == g], d$dct[d$group == "c"])
      ww <- wilcox.test(d$dct[d$group == g], d$dct[d$group == "c"])
      data.frame(plate = d$plate[1], gene = d$gene[1], group = g,
                 n = length(v), dct_mean = mean(v), dct_sd = sd(v),
                 ddct = dd, fold = 2^(-dd),
                 t_p = tt$p.value, w_p = ww$p.value, stringsAsFactors = FALSE)
    }))
  }))
  stats$padj_t <- p.adjust(stats$t_p, method = "BH")
  stats$padj_w <- p.adjust(stats$w_p, method = "BH")
  stats$sig <- ifelse(stats$padj_t < 0.001, "***",
                      ifelse(stats$padj_t < 0.01, "**",
                             ifelse(stats$padj_t < 0.05, "*", "ns")))

  # 复孔 CV
  cv <- do.call(rbind, lapply(split(long, paste(long$plate, long$gene, long$group)),
                              function(d) {
                                data.frame(plate = d$plate[1], gene = d$gene[1],
                                           group = d$group[1], n = nrow(d),
                                           ct_mean = mean(d$ct), ct_sd = sd(d$ct),
                                           cv_pct = 100 * sd(d$ct) / mean(d$ct),
                                           stringsAsFactors = FALSE)
                              }))

  # 写文件
  write.csv(long, file.path(opt$outdir, "qpcr_long.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(stats, file.path(opt$outdir, "qpcr_stats.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(cv, file.path(opt$outdir, "qpcr_cv.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(qc1, file.path(opt$outdir, "plate1_qc_flags.csv"), row.names = FALSE, fileEncoding = "UTF-8")

  # 图
  make_box <- function(plate_df, outname, ylab = "dCt (gene - Gapdh)") {
    plate_df$gene_lab <- gsub("[^\\x20-\\x7E]", "", plate_df$gene)
    p <- ggplot(plate_df, aes(group, dct, color = group)) +
      geom_boxplot(outlier.size = 1) +
      geom_jitter(width = 0.12, size = 1) +
      facet_wrap(~gene_lab, scales = "free_y") +
      scale_color_manual(values = CB_PALETTE, guide = "none") +
      labs(x = NULL, y = ylab) + pub_theme()
    save_pub_fig(p, file.path(opt$outdir, outname))
  }
  make_box(long[long$plate == "P1", ], "fig_p1_dct_box")
  make_box(long[long$plate == "P2", ], "fig_p2_dct_box")
  make_box(long[long$plate == "P3", ], "fig_p3_dct_box")

  cat("DONE\n")
}

if (!interactive()) main()

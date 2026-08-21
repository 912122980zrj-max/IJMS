#!/usr/bin/env Rscript
# 03_DEG_MRG_sets.R —— 差异表达 + MRG 交集（真实实现）
# 对比设计（results/deviation_log.md D1/D2）：
#   人主队列 GSE115857：IgAN_G3（重病变） vs IgAN_G1（轻病变）——组织学严重度轴
#   人验证队列 GSE66494：CKD vs normal
#   小鼠 GSE299417：UUO_Veh vs Sham_Veh（4v4）
# 阈值沿用 IPF 论文：|log2FC|>0.585 且 adj.p<0.05；MRG 集合 = data/reference/MRG_union.txt（2536，用户提供）

suppressPackageStartupMessages({
  library(optparse)
  library(limma)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "03_DEG_MRG_sets.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--outdir"), type = "character", default = DATA_PROC, help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

run_limma <- function(expr, group, levels_use = NULL) {
  group <- factor(group)
  if (!is.null(levels_use)) group <- factor(group, levels = levels_use)
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  fit <- lmFit(expr, design)
  cm <- makeContrasts(contrasts = paste0(levels(group)[1], "-", levels(group)[2]), levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cm))
  topTable(fit2, coef = 1, number = Inf, sort.by = "none")
}

save_deg <- function(tt, path, mrg) {
  out <- data.frame(gene = rownames(tt), tt, check.names = FALSE)
  safe_write(path, function(p) write_csv_utf8(out, p))
  sig <- out[abs(out$logFC) > 0.585 & out$adj.P.Val < 0.05, ]
  intersect(sig$gene, mrg)
}

main <- function() {
  ensure_dirs()
  mrg <- readLines(file.path(DATA_REF, "MRG_union.txt"))
  cat_mrg <- list(
    MM = readLines(file.path(DATA_REF, "MRG_MM_union.txt")),
    CF = readLines(file.path(DATA_REF, "MRG_CF_union.txt")),
    MS = readLines(file.path(DATA_REF, "MRG_MS_union.txt"))
  )

  log_msg("igag", "GSE115857：IgAN_G3 vs IgAN_G1")
  e1 <- readRDS(file.path(opt$outdir, "expr_human_main.rds"))
  m1 <- read.csv(file.path(opt$outdir, "meta_human_main.csv"), stringsAsFactors = FALSE)
  idx <- m1$group %in% c("IgAN_G3", "IgAN_G1")
  tt1 <- run_limma(e1[, idx, drop = FALSE], m1$group[idx], levels_use = c("IgAN_G3", "IgAN_G1"))
  inter1 <- save_deg(tt1, file.path(opt$outdir, "deg_igag_G3vG1.csv"), mrg)
  safe_write(file.path(opt$outdir, "deg_mrg_intersect_igag_g3g1.rds"), function(p) saveRDS(inter1, p))
  log_msg("igag", paste0("DEG n=", nrow(tt1[abs(tt1$logFC) > 0.585 & tt1$adj.P.Val < 0.05, ]),
                         " | MRG overlap=", length(inter1)))

  log_msg("igag_main", "GSE115857：IgAN（进展性） vs 非进展性对照（MCD/MN/活体供肾）")
  m1b <- m1
  m1b$disease <- ifelse(m1b$status == "IgAN patient", "IgAN",
                        ifelse(m1b$status %in% c("minimal change disease",
                                                 "Membranous glomerulonephritis",
                                                 "Living donor"), "control", NA))
  idx1b <- !is.na(m1b$disease)
  tt1b <- run_limma(e1[, idx1b, drop = FALSE], m1b$disease[idx1b],
                    levels_use = c("IgAN", "control"))
  inter1b <- save_deg(tt1b, file.path(opt$outdir, "deg_igag_vs_control.csv"), mrg)
  safe_write(file.path(opt$outdir, "deg_mrg_intersect_igag.rds"),
             function(p) saveRDS(inter1b, p))   # 覆盖：04/下游以该对比为主线
  log_msg("igag_main", paste0("IgAN n=", sum(m1b$disease == "IgAN", na.rm = TRUE),
                              " | control n=", sum(m1b$disease == "control", na.rm = TRUE),
                              " | DEG n=", nrow(tt1b[abs(tt1b$logFC) > 0.585 & tt1b$adj.P.Val < 0.05, ]),
                              " | MRG overlap=", length(inter1b)))

  log_msg("ckd", "GSE66494：CKD vs normal")
  e2 <- readRDS(file.path(opt$outdir, "expr_human_ckd.rds"))
  m2 <- read.csv(file.path(opt$outdir, "meta_human_ckd.csv"), stringsAsFactors = FALSE)
  tt2 <- run_limma(e2, m2$group, levels_use = c("CKD", "normal"))
  inter2 <- save_deg(tt2, file.path(opt$outdir, "deg_ckd.csv"), mrg)
  safe_write(file.path(opt$outdir, "deg_mrg_intersect_ckd.rds"), function(p) saveRDS(inter2, p))
  log_msg("ckd", paste0("DEG n=", nrow(tt2[abs(tt2$logFC) > 0.585 & tt2$adj.P.Val < 0.05, ]),
                        " | MRG overlap=", length(inter2)))

  log_msg("mouse", "GSE299417：UUO_Veh vs Sham_Veh")
  e3 <- readRDS(file.path(opt$outdir, "expr_mouse.rds"))
  m3 <- read.csv(file.path(opt$outdir, "meta_mouse.csv"), stringsAsFactors = FALSE)
  idx3 <- m3$group %in% c("UUO_Veh", "Sham_Veh")
  tt3 <- run_limma(e3[, idx3, drop = FALSE], m3$group[idx3], levels_use = c("UUO_Veh", "Sham_Veh"))
  out3 <- data.frame(gene = rownames(tt3), tt3, check.names = FALSE)
  safe_write(file.path(opt$outdir, "deg_mouse_uuo.csv"), function(p) write_csv_utf8(out3, p))
  log_msg("mouse", paste0("DEG n=", nrow(out3[abs(out3$logFC) > 0.585 & out3$adj.P.Val < 0.05, ])))

  log_msg("deg", "DONE")
}

if (!interactive()) main()

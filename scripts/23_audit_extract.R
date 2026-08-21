#!/usr/bin/env Rscript
# 23_audit_extract.R —— 手稿数据核查的 R 侧提取器（只读，不改动任何结果文件）
# 目的：把仅存于 RDS/重算才能得到的数字导出为 JSON，供 22_manuscript_data_audit.py 对账：
#   WGCNA 软阈值 R2 / 模块数 / 最优模块 |cor| / 核心基因数 / 3-gene signature /
#   nomogram C-index 与校准斜率 / ConsensusClusterPlus PAC 与 C1-C2 样本数。

suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
  library(rms)
  library(MASS)
  library(ConsensusClusterPlus)
  library(jsonlite)
})

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
this_script <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else getwd()
project_root <- dirname(dirname(normalizePath(this_script)))

option_list <- list(
  make_option(c("--root"), type = "character", default = project_root,
              help = "项目根目录（默认取脚本上一级）"),
  make_option(c("--out"), type = "character",
              default = file.path(project_root, "data", "processed", "audit_r_items.json"),
              help = "输出 JSON 路径")
)
opt <- parse_args(OptionParser(option_list = option_list))

pac_score <- function(cm) {
  x <- cm[lower.tri(cm)]
  mean(x > 0.1 & x < 0.9)
}

main <- function() {
  root <- normalizePath(opt$root)
  dp <- file.path(root, "data", "processed")

  expr <- readRDS(file.path(dp, "expr_human_main.rds"))
  meta <- read.csv(file.path(dp, "meta_human_main.csv"), stringsAsFactors = FALSE)
  keep <- meta$status == "IgAN patient"
  datExpr <- t(expr[, keep, drop = FALSE])
  score <- read.csv(file.path(dp, "ssgsea_mech_scores.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  trait_by_pos <- score$MRG_up[keep]
  trait_by_name <- score$MRG_up[match(meta$sample[keep], score$sample)]
  stopifnot(identical(as.numeric(trait_by_pos), as.numeric(trait_by_name)))
  trait <- as.numeric(trait_by_name)

  # WGCNA soft threshold scan (recomputed, mirrors 05)
  set.seed(42)
  sft <- pickSoftThreshold(datExpr, powerVector = 1:10, networkType = "signed")
  sft_r2 <- setNames(sft$fitIndices$SFT.R.sq, sft$fitIndices$Power)
  beta <- ifelse(is.na(sft$powerEstimate), 4, sft$powerEstimate)

  # Modules from the saved object (no re-run of blockwiseModules)
  net <- readRDS(file.path(dp, "wgcna_modules.rds"))
  n_modules <- length(unique(net$colors))
  me_cor <- abs(as.numeric(WGCNA::cor(net$MEs, trait, use = "pairwise.complete.obs")))
  best_idx <- which.max(me_cor)

  core <- readRDS(file.path(dp, "wgcna_core.rds"))
  sig <- readRDS(file.path(dp, "progression_signature.rds"))

  # Nomogram / calibration, mirrors 07（主口径：排除 FSGS 的 4 类状态）
  # 另算 FSGS-并入对照组口径（18_figure_final 的 idx 逻辑）用于定位 0.696 的来源。
  nomo_run <- function(meta_sub) {
    m <- meta[meta$status %in% meta_sub, , drop = FALSE]
    d <- data.frame(IgAN = as.integer(m$status == "IgAN patient"),
                    t(expr[intersect(sig, rownames(expr)), m$sample, drop = FALSE]))
    f <- glm(IgAN ~ ., data = d, family = binomial())
    s <- MASS::stepAIC(f, direction = "both", trace = 0)
    assign("dd", datadist(d), envir = .GlobalEnv)
    options(datadist = "dd")
    lf <- lrm(IgAN ~ ., data = d, x = TRUE, y = TRUE)
    pr <- predict(lf, type = "fitted")
    cc <- val.prob(pr, d$IgAN, pl = FALSE)
    list(n = nrow(d), terms = names(coef(s))[-1],
         dxy = unname(cc["Dxy"]), c_index = unname(cc["C (ROC)"]),
         slope = unname(cc["Slope"]))
  }
  nomo_main <- nomo_run(c("IgAN patient", "minimal change disease",
                          "Membranous glomerulonephritis", "Living donor"))
  nomo_fsgs_in <- nomo_run(c("IgAN patient", "minimal change disease",
                             "Membranous glomerulonephritis", "Living donor",
                             "Focal Segmental Glomerulosclerosis"))

  # Consensus clustering PAC, mirrors 08 (3 signature genes x 55 IgAN samples)
  dat_cc <- as.matrix(expr[intersect(sig, rownames(expr)), keep, drop = FALSE])
  set.seed(42)
  res <- ConsensusClusterPlus(dat_cc, maxK = 9, reps = 1000, pItem = 0.8,
                              pFeature = 0.8, clusterAlg = "hc",
                              distance = "pearson", seed = 42)
  pacs <- sapply(res[2:9], function(x) pac_score(x$consensusMatrix))
  cls <- res[[2]]$consensusClass

  # 亚型比较用 R 的 wilcox.test 原口径（含连续性校正），与 08/18 脚本一致
  cls_for_p <- cls
  sub_score <- score$MRG_up[match(names(cls_for_p), score$sample)]
  sub_wilcox <- wilcox.test(sub_score ~ factor(paste0("C", cls_for_p)))$p.value
  emt <- as.numeric(read.csv(file.path(dp, "hallmark_scores.csv"),
                             check.names = FALSE)[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]])
  emt_name <- read.csv(file.path(dp, "hallmark_scores.csv"),
                       check.names = FALSE)$sample
  names(emt) <- emt_name
  sub_emt <- wilcox.test(emt[names(cls_for_p)] ~ factor(paste0("C", cls_for_p)))$p.value

  out <- list(
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "Asia/Shanghai"),
    wgcna = list(
      n_igan = sum(keep),
      power_estimate = beta,
      sft_r2_at_power4 = unname(sft_r2[["4"]]),
      sft_r2 = as.list(unname(sft_r2)),
      n_modules = n_modules,
      best_module = colnames(net$MEs)[best_idx],
      best_module_abs_cor = me_cor[best_idx],
      n_core = length(core)
    ),
    signature = as.list(sig),
    nomogram = nomo_main,
    nomogram_fsgs_included = nomo_fsgs_in,
    clustering = list(
      pac_k2_to_k9 = as.list(unname(pacs)),
      c1_n = sum(cls == 1),
      c2_n = sum(cls == 2),
      subtype_mech_wilcox_p = unname(sub_wilcox),
      subtype_emt_wilcox_p = unname(sub_emt)
    )
  )
  write_json(out, opt$out, pretty = TRUE, auto_unbox = TRUE, digits = 8)
  message("[audit-extract] DONE -> ", opt$out)
}

if (!interactive()) main()

#!/usr/bin/env Rscript
# 47_benchmark_enrichment.R —— 对标模板补做：GO/KEGG 富集 + 面板数据导出
#
# 输出:
#   results/benchmark/go_bp.csv / go_cc.csv / go_mf.csv / kegg.csv
#   results/benchmark/core_expr.csv      (100 核心基因 x 86 样本)
#   results/benchmark/sig_expr.csv       (NDNF/PCDHB7/RRAGB x 86 样本)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(DOSE)
})

ROOT <- "E:/sheng xin/ObstructiveNephropathy_MRG"
PROC <- file.path(ROOT, "data", "processed")
OUT <- file.path(ROOT, "results", "benchmark")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

core <- readRDS(file.path(PROC, "wgcna_core.rds"))
cat("core genes:", length(core), "\n")

write.csv(data.frame(gene = core, stringsAsFactors = FALSE),
          file.path(OUT, "core_gene_list.csv"), row.names = FALSE)

ids <- bitr(core, fromType = "SYMBOL", toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)
cat("mapped to ENTREZID:", nrow(ids), "\n")

run_enrich <- function(fun, ...) {
  tryCatch(fun(...), error = function(e) {
    message("enrichment failed: ", conditionMessage(e))
    NULL
  })
}

go_bp <- run_enrich(enrichGO, gene = ids$ENTREZID, OrgDb = org.Hs.eg.db,
                    ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                    readable = TRUE)
go_cc <- run_enrich(enrichGO, gene = ids$ENTREZID, OrgDb = org.Hs.eg.db,
                    ont = "CC", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                    readable = TRUE)
go_mf <- run_enrich(enrichGO, gene = ids$ENTREZID, OrgDb = org.Hs.eg.db,
                    ont = "MF", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                    readable = TRUE)
kegg <- run_enrich(enrichKEGG, gene = ids$ENTREZID, organism = "hsa",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2)

save_result <- function(obj, path) {
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    write.csv(as.data.frame(obj), path, row.names = FALSE)
    cat(basename(path), "rows:", nrow(as.data.frame(obj)), "\n")
  } else {
    write.csv(data.frame(note = "no significant terms"), path,
              row.names = FALSE)
    cat(basename(path), ": no significant terms\n")
  }
}

save_result(go_bp, file.path(OUT, "go_bp.csv"))
save_result(go_cc, file.path(OUT, "go_cc.csv"))
save_result(go_mf, file.path(OUT, "go_mf.csv"))
save_result(kegg, file.path(OUT, "kegg.csv"))

# 核心基因仅名义显著：另存宽阈值结果供如实展示
go_bp_nom <- run_enrich(enrichGO, gene = ids$ENTREZID, OrgDb = org.Hs.eg.db,
                        ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 1,
                        readable = TRUE)
kegg_nom <- run_enrich(enrichKEGG, gene = ids$ENTREZID, organism = "hsa",
                       pvalueCutoff = 0.05, qvalueCutoff = 1)
save_result(go_bp_nom, file.path(OUT, "go_bp_core_nominal.csv"))
save_result(kegg_nom, file.path(OUT, "kegg_core_nominal.csv"))

# ---- 差异基因富集（231 DEG 与 36 MRG∩DEG） ----
deg <- read.csv(file.path(PROC, "deg_igag_vs_control.csv"))
deg_sig <- deg[abs(deg$logFC) > 0.585 & deg$adj.P.Val < 0.05 & !is.na(deg$adj.P.Val), ]
cat("DEG n:", nrow(deg_sig),
    "(up:", sum(deg_sig$logFC > 0), "down:", sum(deg_sig$logFC < 0), ")\n")
mrg <- readLines(file.path(ROOT, "submission", "ijms", "supporting_information",
                           "S1_Table_MRG_gene_set.txt"))
mrg_deg <- deg_sig[deg_sig$gene %in% mrg, ]
cat("MRG ∩ DEG n:", nrow(mrg_deg), "\n")

enrich_list <- function(genes) {
  ids2 <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
  list(
    bp = run_enrich(enrichGO, gene = ids2$ENTREZID, OrgDb = org.Hs.eg.db,
                    ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                    readable = TRUE),
    kegg = run_enrich(enrichKEGG, gene = ids2$ENTREZID, organism = "hsa",
                      pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  )
}

res_deg <- enrich_list(deg_sig$gene)
save_result(res_deg$bp, file.path(OUT, "go_bp_deg.csv"))
save_result(res_deg$kegg, file.path(OUT, "kegg_deg.csv"))

res_mrg <- enrich_list(mrg_deg$gene)
save_result(res_mrg$bp, file.path(OUT, "go_bp_mrg_deg.csv"))
save_result(res_mrg$kegg, file.path(OUT, "kegg_mrg_deg.csv"))

# ---- 导出表达矩阵 ----
expr <- readRDS(file.path(PROC, "expr_human_main.rds"))
core_avail <- intersect(core, rownames(expr))
sig <- intersect(c("NDNF", "PCDHB7", "RRAGB"), rownames(expr))
cat("core available in matrix:", length(core_avail), "\n")
cat("signature available:", paste(sig, collapse = ", "), "\n")

write.csv(data.frame(gene = rownames(expr[core_avail, , drop = FALSE]),
                     expr[core_avail, , drop = FALSE], check.names = FALSE),
          file.path(OUT, "core_expr.csv"), row.names = FALSE)
write.csv(data.frame(gene = rownames(expr[sig, , drop = FALSE]),
                     expr[sig, , drop = FALSE], check.names = FALSE),
          file.path(OUT, "sig_expr.csv"), row.names = FALSE)

cat("DONE\n")

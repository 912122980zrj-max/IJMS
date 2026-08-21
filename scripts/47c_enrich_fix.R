#!/usr/bin/env Rscript
# 单独重算 231 DEG 与 100 核心基因的 GO BP/KEGG，避免 47 中多线程/顺序问题
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(DOSE)
})

ROOT <- "E:/sheng xin/ObstructiveNephropathy_MRG"
PROC <- file.path(ROOT, "data", "processed")
OUT <- file.path(ROOT, "results", "benchmark")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

save_or_note <- function(obj, path) {
  d <- as.data.frame(obj)
  if (!is.null(obj) && nrow(d) > 0) {
    write.csv(d, path, row.names = FALSE)
    cat(basename(path), "rows:", nrow(d), "min padj:", min(d$p.adjust), "\n")
  } else {
    write.csv(data.frame(note = "no terms passed threshold"), path,
              row.names = FALSE)
    cat(basename(path), ": none\n")
  }
}

# --- 231 DEG ---
deg <- read.csv(file.path(PROC, "deg_igag_vs_control.csv"))
deg_sig <- deg[abs(deg$logFC) > 0.585 & deg$adj.P.Val < 0.05
               & !is.na(deg$adj.P.Val), ]
ids_deg <- bitr(deg_sig$gene, fromType = "SYMBOL", toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
cat("DEG mapped:", nrow(ids_deg), "\n")
go_deg <- enrichGO(gene = ids_deg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                   pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
save_or_note(go_deg, file.path(OUT, "go_bp_deg_full.csv"))
k_deg <- enrichKEGG(gene = ids_deg$ENTREZID, organism = "hsa",
                    pvalueCutoff = 1, qvalueCutoff = 1)
save_or_note(k_deg, file.path(OUT, "kegg_deg_full.csv"))

# --- 100 core (名义) ---
core <- readRDS(file.path(PROC, "wgcna_core.rds"))
ids_core <- bitr(core, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
cat("core mapped:", nrow(ids_core), "\n")
go_core <- enrichGO(gene = ids_core$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                    pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
save_or_note(go_core, file.path(OUT, "go_bp_core_full.csv"))
k_core <- enrichKEGG(gene = ids_core$ENTREZID, organism = "hsa",
                     pvalueCutoff = 1, qvalueCutoff = 1)
save_or_note(k_core, file.path(OUT, "kegg_core_full.csv"))
cat("DONE\n")

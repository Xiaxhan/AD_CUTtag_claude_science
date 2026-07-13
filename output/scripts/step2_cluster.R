## Step 2: de-novo RNA clustering + annotation, per region
## Usage: Rscript step2_cluster.R <REGION>   (REGION = HIP or CBL)
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(data.table); library(harmony)})
args <- commandArgs(trailingOnly=TRUE); REGION <- args[1]
stopifnot(REGION %in% c("HIP","CBL"))
set.seed(1)
base <- "/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/output"
cache <- file.path(base,"cache"); dir.create(file.path(base,"cluster"), showWarnings=FALSE)

# Load cached QC libs whose cells' eff_region == REGION (a lib can contribute to either region post-correction)
files <- list.files(cache,"_qc.rds$",full.names=TRUE)
mats <- list(); metas <- list()
for(f in files){
  x <- readRDS(f)
  sel <- x$meta$eff_region==REGION
  if(!any(sel)) next
  m <- x$rna[, sel, drop=FALSE]; md <- x$meta[sel,,drop=FALSE]
  mats[[x$id]] <- m; metas[[x$id]] <- md
  cat("  +", x$id, ":", ncol(m), "cells to", REGION, "\n")
}
# common genes (identical across libs, but be safe)
g <- Reduce(intersect, lapply(mats, rownames))
mats <- lapply(mats, function(m) m[g,,drop=FALSE])
counts <- do.call(cbind, mats)
meta <- rbindlist(metas); meta <- as.data.frame(meta); rownames(meta) <- colnames(counts)
cat(REGION, "merged:", nrow(counts), "genes x", ncol(counts), "cells\n")

so <- CreateSeuratObject(counts=counts, meta.data=meta)
so <- NormalizeData(so, verbose=FALSE)
so <- FindVariableFeatures(so, nfeatures=2000, verbose=FALSE)
so <- ScaleData(so, verbose=FALSE)
so <- RunPCA(so, npcs=30, verbose=FALSE)
# integrate over library (removes batch + mark-technical differences in RNA)
so <- RunHarmony(so, group.by.vars="library", verbose=FALSE)
so <- FindNeighbors(so, reduction="harmony", dims=1:30, verbose=FALSE)
so <- FindClusters(so, resolution=0.6, algorithm=1, verbose=FALSE)  # Louvain (Leiden needs py); res 0.6
so <- RunUMAP(so, reduction="harmony", dims=1:30, verbose=FALSE)
cat("clusters:", nlevels(so$seurat_clusters), "\n"); print(table(so$seurat_clusters))

saveRDS(so, file.path(base,"cluster",paste0(REGION,"_rna_clustered.rds")))
# markers
mk <- FindAllMarkers(so, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.5, verbose=FALSE)
fwrite(mk, file.path(base,"cluster",paste0(REGION,"_cluster_markers_all.csv")))
cat("DONE", REGION, "\n")

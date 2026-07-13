## Step 2b: annotate clusters -> major cell types via AddModuleScore, per region
suppressPackageStartupMessages({library(Seurat); library(data.table); library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE); REGION <- args[1]
base <- "/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/output"
so <- readRDS(file.path(base,"cluster",paste0(REGION,"_rna_clustered.rds")))
feat <- readRDS(file.path(base,"cluster","feature_map.rds"))
sym2ensg <- setNames(feat$ensg, feat$symbol)

# Region-appropriate, specificity-weighted marker panels (symbols)
common <- list(
  ExN   = c("SLC17A7","SATB2","SLC17A6","RORB","FEZF2"),
  InN   = c("GAD1","GAD2","SLC32A1","LHX6","ADARB2"),
  Astro = c("AQP4","GFAP","SLC1A2","GJA1","ALDH1L1"),
  Oligo = c("MOBP","PLP1","MOG","MBP","ST18"),
  OPC   = c("PDGFRA","CSPG4","LHFPL3","OLIG1"),
  Micro = c("CSF1R","P2RY12","CX3CR1","C1QA","DOCK8"),
  Endo  = c("CLDN5","FLT1","PECAM1","VWF"),
  Mural = c("PDGFRB","RGS5","ACTA2","DCN"),          # pericyte/VSMC/fibroblast
  Epend = c("FOXJ1","PIFO","TMEM212","CCDC153")
)
region_specific <- if(REGION=="HIP") list(
  DG    = c("PROX1","C1QL2","DSP"),                  # dentate granule (ExN subtype)
  CA    = c("FIBCD1","NECAB1","CABP7")               # CA pyramidal (ExN subtype)
) else list(
  Granule  = c("GABRA6","CNTN2","NEUROD1","RBFOX3"), # cerebellar granule neurons
  Purkinje = c("PCP2","PCP4","CALB1","ITPR1","CA8"),
  Bergmann = c("GDF10","HEPACAM","AQP4"),            # Bergmann glia (specialized astro)
  UBC_Golgi= c("EOMES","GRM1","LGI2")
)
panel <- c(common, region_specific)
panel_ensg <- lapply(panel, function(v){ e<-sym2ensg[v]; e<-e[!is.na(e)]; e[e %in% rownames(so)] })
panel_ensg <- panel_ensg[sapply(panel_ensg, length)>=2]

so <- AddModuleScore(so, features=panel_ensg, name="ms_", seed=1)
ms_cols <- paste0("ms_", seq_along(panel_ensg))
names(ms_cols) <- names(panel_ensg)
# per-cluster mean module score
cl <- so$seurat_clusters
M <- sapply(ms_cols, function(c) tapply(so[[c]][,1], cl, mean))
colnames(M) <- names(panel_ensg)          # rows=clusters, cols=types
# z-score across clusters within each type, then assign per cluster by max
Mz <- scale(M)
assign_type <- colnames(Mz)[apply(Mz,1,which.max)]
top1 <- apply(Mz,1,max); srt <- t(apply(Mz,1,sort,decreasing=TRUE)); margin <- srt[,1]-srt[,2]
ann <- data.frame(cluster=rownames(M), type=assign_type, top_z=round(top1,2),
                  margin=round(margin,2), n=as.integer(table(cl)[rownames(M)]))
# collapse subtypes to major type for the cross-region major label
major_map <- c(DG="ExN", CA="ExN", Granule="GranuleN", Purkinje="Purkinje",
               Bergmann="Astro", UBC_Golgi="InN", Mural="Vascular", Endo="Vascular", Epend="Ependymal")
ann$major <- ifelse(ann$type %in% names(major_map), major_map[ann$type], ann$type)
ann$subtype <- ann$type
write.csv(ann, file.path(base,"cluster",paste0(REGION,"_cluster_annotation.csv")), row.names=FALSE)
write.csv(round(Mz,3), file.path(base,"cluster",paste0(REGION,"_clusterZ_by_type.csv")))
# attach to object
so$major_type <- ann$major[match(as.character(cl), ann$cluster)]
so$subtype    <- ann$subtype[match(as.character(cl), ann$cluster)]
saveRDS(so, file.path(base,"cluster",paste0(REGION,"_rna_annotated.rds")))
cat("=== ",REGION," annotation ===\n"); print(ann)
cat("\nmajor-type cell counts:\n"); print(table(so$major_type))

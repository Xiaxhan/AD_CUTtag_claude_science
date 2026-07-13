## Step 2c: annotation figures (UMAP + marker dotplot) per region
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(data.table); library(patchwork)})
args <- commandArgs(trailingOnly=TRUE); REGION <- args[1]
base <- "/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/output"
so <- readRDS(file.path(base,"cluster",paste0(REGION,"_rna_annotated.rds")))
feat <- readRDS(file.path(base,"cluster","feature_map.rds"))
sym2ensg <- setNames(feat$ensg, feat$symbol)

pal <- c(ExN="#1f77b4", InN="#17becf", Astro="#2ca02c", Oligo="#8c564b", OPC="#e377c2",
         Micro="#d62728", Vascular="#9467bd", Ependymal="#bcbd22", GranuleN="#ff7f0e",
         Purkinje="#e41a1c")
so$major_type <- factor(so$major_type)

# UMAP by major type
d <- data.table(so@reductions$umap@cell.embeddings); setnames(d, c("UMAP1","UMAP2"))
d$major <- so$major_type; d$mark <- so$mark
set.seed(1); d <- d[sample(.N)]
cen <- d[, .(x=median(UMAP1), y=median(UMAP2)), by=major]
pU <- ggplot(d, aes(UMAP1, UMAP2, color=major)) +
  geom_point(size=0.15, alpha=0.5, stroke=0) +
  scale_color_manual(values=pal, name=NULL) +
  ggrepel::geom_text_repel(data=cen, aes(x,y,label=major), color="black", size=3, seed=1, bg.color="white", bg.r=0.15) +
  guides(color=guide_legend(override.aes=list(size=2,alpha=1))) +
  labs(title=paste0(REGION,": de-novo RNA clusters (n=",ncol(so)," nuclei)")) +
  theme_bw(base_size=9) + theme(panel.grid=element_blank())

# UMAP by mark (technical check: marks should co-mingle after Harmony)
pM <- ggplot(d, aes(UMAP1, UMAP2, color=mark)) +
  geom_point(size=0.15, alpha=0.5, stroke=0) +
  scale_color_manual(values=c(H3K27ac="#D6604D",H3K27me3="#4393C3"), name=NULL) +
  guides(color=guide_legend(override.aes=list(size=2,alpha=1))) +
  labs(title="Integration check (by histone mark)") +
  theme_bw(base_size=9) + theme(panel.grid=element_blank())

# Marker dotplot
dp_sym <- if(REGION=="HIP")
  c("SLC17A7","SATB2","PROX1","FIBCD1","GAD1","GAD2","AQP4","GFAP","MOBP","PLP1",
    "PDGFRA","CSPG4","CSF1R","P2RY12","CLDN5","FLT1","PDGFRB","RGS5","FOXJ1","PIFO") else
  c("GABRA6","CNTN2","PCP2","PCP4","CALB1","SLC17A7","GAD1","GAD2","AQP4","GDF10",
    "MOBP","PLP1","PDGFRA","CSF1R","P2RY12","CLDN5","FLT1","PDGFRB")
dp_ensg <- sym2ensg[dp_sym]; ok <- !is.na(dp_ensg) & dp_ensg %in% rownames(so)
dp_ensg <- dp_ensg[ok]; dp_lab <- dp_sym[ok]
Idents(so) <- factor(so$major_type)
# Extract DotPlot data and rebuild cleanly: cell types on X (readable), markers on Y
dpd <- DotPlot(so, features=dp_ensg)$data
dpd$gene <- factor(setNames(dp_lab, dp_ensg)[as.character(dpd$features.plot)], levels=dp_lab)
dpd$id   <- factor(dpd$id)
pD <- ggplot(dpd, aes(x=id, y=gene, size=pct.exp, color=avg.exp.scaled)) +
  geom_point() +
  scale_size_continuous(range=c(0,6), name="% expressed") +
  scale_color_gradient2(low="#4393C3", mid="grey90", high="#D6604D", midpoint=0, name="scaled\nexpression") +
  labs(x=NULL, y="canonical marker", title=paste0(REGION,": marker expression by major cell type")) +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
        axis.text.y=element_text(size=11, face="italic"),
        legend.title=element_text(size=11), legend.text=element_text(size=10),
        plot.title=element_text(size=13))

# enlarge fonts on the two UMAP panels
bigify <- theme(plot.title=element_text(size=13), legend.text=element_text(size=12),
                axis.title=element_text(size=12), axis.text=element_text(size=10))
pU <- pU + bigify; pM <- pM + bigify

fig <- (pU | pM) / pD + plot_layout(heights=c(1,1.15)) + plot_annotation(tag_levels="a")
ggsave(file.path(base,paste0("celltypes_",REGION,".png")), fig, width=12, height=10.5, dpi=200)
cat("saved celltypes_",REGION,".png\n", sep="")

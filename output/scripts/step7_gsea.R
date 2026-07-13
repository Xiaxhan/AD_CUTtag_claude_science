## Step 7: rank-based GSEA (fgsea) on gene-level AD-vs-CTL rankings — pooled-power enrichment
suppressPackageStartupMessages({library(data.table); library(fgsea)})
base <- "/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/output"
gde <- fread(file.path(base,"DE_genelevel_results.csv.gz"))
# read GMT gene sets
read_gmt <- function(f){
  ln <- readLines(f); sets <- list()
  for(l in ln){ x<-strsplit(l,"\t")[[1]]; if(length(x)>=3) sets[[x[1]]] <- unique(x[3:length(x)][x[3:length(x)]!=""]) }
  sets
}
gmts <- list(Hallmark="/tmp/genesets/Hallmark.gmt", KEGG="/tmp/genesets/KEGG.gmt",
             GO_BP="/tmp/genesets/GO_BP.gmt", WikiPathways="/tmp/genesets/WikiPathways.gmt")
allsets <- do.call(c, lapply(names(gmts), function(n){ s<-read_gmt(gmts[[n]]); names(s)<-paste0(n,":",names(s)); s }))
cat("total gene sets:", length(allsets), "\n")

# rank metric per stratum: signed -log10(P) * sign(logFC); test well-powered strata
strata <- gde[, .N, by=.(stratum,mark)][N>=2000]  # only strata with enough genes tested
res <- list()
set.seed(1)
for(i in seq_len(nrow(strata))){
  st<-strata$stratum[i]; mk<-strata$mark[i]
  d <- gde[stratum==st & mark==mk & !is.na(PValue)]
  d[, stat := sign(logFC) * -log10(PValue)]
  d <- d[is.finite(stat)]
  ranks <- d$stat; names(ranks) <- d$gene
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- sort(ranks)
  fg <- fgsea(pathways=allsets, stats=ranks, minSize=10, maxSize=500, nPermSimple=10000)
  fg <- as.data.table(fg)[, .(pathway, pval, padj, NES, size, stratum=st, mark=mk)]
  res[[paste(st,mk)]] <- fg
  ns <- fg[padj<0.1, .N]
  cat(sprintf("  %s %s: %d sets tested, %d padj<0.1\n", st, mk, nrow(fg), ns))
}
gsea <- rbindlist(res)
gsea[, leadingEdge := NULL] -> tmp   # avoid list-col in fwrite (already dropped in select)
fwrite(gsea, file.path(base,"GSEA_results.csv.gz"))
saveRDS(gsea, file.path(base,"cache","gsea_results.rds"))
cat("\n=== Top enriched pathways (padj<0.1) ===\n")
print(gsea[padj<0.1][order(padj)][1:min(30,.N), .(stratum,mark,pathway=substr(pathway,1,50),NES=round(NES,2),padj=round(padj,3))])
cat("\nTotal significant pathway-stratum hits:", nrow(gsea[padj<0.1]), "\n")

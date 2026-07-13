## Step 5b: gene-level CRE aggregation DE (pooled power) — sum CREs per gene, test AD vs CTL (edgeR QL)
suppressPackageStartupMessages({library(Matrix); library(data.table); library(edgeR)})
base <- "/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/output"
Mac  <- readRDS(file.path(base,"cache","pb_merged_H3K27ac.rds"))
Mme3 <- readRDS(file.path(base,"cache","pb_merged_H3K27me3.rds"))
pk_ac  <- fread(file.path(base,"cache/peakann_H3K27ac.csv.gz"))
pk_me3 <- fread(file.path(base,"cache/peakann_H3K27me3.csv.gz"))

# build peak->gene aggregation matrix (assign each CRE to nearest gene within 100kb of TSS OR genic)
gene_agg <- function(M, pk){
  pk <- pk[!is.na(nearest_gene) & (context %in% c("promoter","genic") | dist_to_tss<=1e5)]
  pk <- pk[peak %in% rownames(M$counts)]
  genes <- sort(unique(pk$nearest_gene))
  gi <- match(pk$nearest_gene, genes); pj <- match(pk$peak, rownames(M$counts))
  A <- sparseMatrix(i=pj, j=gi, x=1, dims=c(nrow(M$counts), length(genes)), dimnames=list(NULL, genes))
  genemat <- t(A) %*% M$counts   # genes x groups (sum of CRE counts per gene)
  genemat
}
gm_ac  <- gene_agg(Mac, pk_ac)
gm_me3 <- gene_agg(Mme3, pk_me3)
cat("gene-level matrices: ac", dim(gm_ac), " me3", dim(gm_me3), "\n")

run_gene_de <- function(genemat, M, mark){
  cm <- copy(M$colmeta); cm[, stratum:=paste(region,major_type,sep=".")]; cm[,valid:=ncells>=25]; cm[,col:=.I]
  strata <- cm[valid==TRUE, .(nAD=sum(dx=="AD"),nCTL=sum(dx=="CTL")), by=stratum][nAD>=3&nCTL>=3, stratum]
  out <- list()
  for(st in strata){
    sub <- cm[stratum==st & valid==TRUE]; y <- genemat[, sub$col, drop=FALSE]
    grp <- factor(sub$dx, levels=c("CTL","AD"))
    d<-DGEList(as.matrix(y),group=grp); keep<-filterByExpr(d,group=grp); d<-d[keep,,keep.lib.sizes=FALSE]; d<-calcNormFactors(d)
    des<-model.matrix(~grp); d<-estimateDisp(d,des); f<-glmQLFit(d,des); q<-glmQLFTest(f,coef=2)
    t<-as.data.table(topTags(q,n=Inf)$table); t$gene<-rownames(topTags(q,n=Inf)$table)
    t[, `:=`(stratum=st, mark=mark)]
    out[[st]]<-t[, .(gene,logFC,logCPM,PValue,FDR,stratum,mark)]
  }
  rbindlist(out)
}
gde <- rbind(run_gene_de(gm_ac,Mac,"H3K27ac"), run_gene_de(gm_me3,Mme3,"H3K27me3"))
fwrite(gde, file.path(base,"DE_genelevel_results.csv.gz"))
saveRDS(list(gde=gde, gm_ac=gm_ac, gm_me3=gm_me3), file.path(base,"cache","genelevel_de.rds"))
# summary
s <- gde[, .(n_genes=.N, sig05=sum(FDR<0.05 & abs(logFC)>=0.5), sig10=sum(FDR<0.1 & abs(logFC)>=0.25)), by=.(stratum,mark)]
cat("\n=== Gene-level DE summary (sig strata) ===\n"); print(s[sig10>0][order(-sig10)])
cat("\nTotal gene-level primary sig:", sum(s$sig05), " relaxed:", sum(s$sig10), "\n")
cat("\nTop gene-level hits:\n")
print(gde[FDR<0.1 & abs(logFC)>=0.25][order(FDR)][1:20, .(stratum,mark,gene,logFC=round(logFC,2),FDR=round(FDR,3))])

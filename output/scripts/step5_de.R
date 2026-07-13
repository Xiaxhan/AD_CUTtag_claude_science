## Step 5: pseudobulk differential CRE AD vs CTL, per region x celltype x mark (edgeR QL)
## Donor = replication unit. Never per-cell. Primary sig FDR<0.05 & |lfc|>=0.5; relaxed FDR<0.1 & |lfc|>=0.25.
suppressPackageStartupMessages({library(Matrix); library(data.table); library(edgeR)})
base <- "/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/output"
outdir <- file.path(base,"DE"); dir.create(outdir, showWarnings=FALSE)
ann_meta <- fread(file.path(base,"cell_annotations.csv.gz"))
dl_region <- readRDS(file.path(base,"cache","donor_library_region.rds"))
# batch = source library per (donor, region, mark)
donlib <- unique(ann_meta[, .(donor, region, mark, library)])

run_de <- function(M, mark_arg){
  mark <- mark_arg
  cm <- copy(M$colmeta); cm[, stratum := paste(region, major_type, sep=".")]
  cm[, valid := ncells>=25]
  cm[, col := .I]
  res_all <- list(); summ <- list()
  strata <- cm[valid==TRUE, .(nAD=sum(dx=="AD"), nCTL=sum(dx=="CTL")), by=stratum][nAD>=3 & nCTL>=3, stratum]
  for(st in strata){
    sub <- cm[stratum==st & valid==TRUE]
    y <- M$counts[, sub$col, drop=FALSE]
    grp <- factor(sub$dx, levels=c("CTL","AD"))
    # batch covariate: library per donor in this region/mark
    reg <- sub$region[1]
    lib <- donlib[donor %in% sub$donor & region==reg & mark==mark_arg, .(donor,library)]
    sub <- merge(sub, unique(lib), by="donor", all.x=TRUE, sort=FALSE)
    sub <- sub[match(colnames(y) , key)]  # realign to y columns
    batch <- factor(sub$library)
    # DGEList + filter + TMM
    d <- DGEList(counts=as.matrix(y), group=grp)
    keep <- filterByExpr(d, group=grp)
    d <- d[keep,,keep.lib.sizes=FALSE]; d <- calcNormFactors(d)
    # design: batch (library) is near-confounded with dx and has up to 4 levels for 8-9 donors,
    # leaving too few residual df and absorbing signal (Decision #08). Include batch ONLY when it
    # adds >=2 residual df beyond ~batch+dx AND does not separate dx; else diagnosis-only.
    # Batch (library) here is near-confounded with diagnosis and has levels with only 1 donor
    # (Decision #08). Adjusting for it manufactures spurious "significant" hits (e.g. GSTM1) that
    # vanish under the diagnosis-only model — a batch-confounding-as-signal artifact. We therefore
    # DO NOT adjust for batch: the primary DE model is diagnosis-only (~grp). Batch structure is
    # reported in QC/records, not used as a covariate given this design.
    use_batch <- FALSE
    design <- model.matrix(~grp)
    d <- estimateDisp(d, design)
    fit <- glmQLFit(d, design)
    qlf <- glmQLFTest(fit, coef=ncol(design))  # last coef = AD vs CTL
    tt <- topTags(qlf, n=Inf)$table
    tt$peak <- rownames(tt)
    tt <- as.data.table(tt)[, .(peak, logFC, logCPM, F, PValue, FDR)]
    tt[, `:=`(stratum=st, mark=mark, region=reg, celltype=sub("^[^.]*\\.","",st),
              nAD=sum(grp=="AD"), nCTL=sum(grp=="CTL"), batch_adj=use_batch, n_tested=nrow(tt))]
    res_all[[st]] <- tt
    sig1 <- tt[FDR<0.05 & abs(logFC)>=0.5]; sig2 <- tt[FDR<0.10 & abs(logFC)>=0.25]
    summ[[st]] <- data.table(stratum=st, mark=mark, region=reg, celltype=sub("^[^.]*\\.","",st),
        nAD=sum(grp=="AD"), nCTL=sum(grp=="CTL"), n_tested=nrow(tt), batch_adj=use_batch,
        n_sig_primary=nrow(sig1), n_up_primary=sum(sig1$logFC>0), n_dn_primary=sum(sig1$logFC<0),
        n_sig_relaxed=nrow(sig2))
    cat(sprintf("  %s %s: tested=%d primary_sig=%d (up %d/dn %d) batch_adj=%s\n",
        mark, st, nrow(tt), nrow(sig1), sum(sig1$logFC>0), sum(sig1$logFC<0), use_batch))
  }
  list(res=rbindlist(res_all), summ=rbindlist(summ))
}
Mac  <- readRDS(file.path(base,"cache","pb_merged_H3K27ac.rds"))
Mme3 <- readRDS(file.path(base,"cache","pb_merged_H3K27me3.rds"))
cat("=== H3K27ac DE ===\n"); DAC  <- run_de(Mac,"H3K27ac")
cat("=== H3K27me3 DE ===\n"); DME3 <- run_de(Mme3,"H3K27me3")
allres <- rbind(DAC$res, DME3$res); allsumm <- rbind(DAC$summ, DME3$summ)
fwrite(allres,  file.path(base,"DE_results_all.csv.gz"))
fwrite(allsumm, file.path(base,"DE_summary.csv"))
saveRDS(list(res=allres, summ=allsumm), file.path(base,"cache","de_results.rds"))
cat("\nDONE DE. Total sig (primary):", sum(allsumm$n_sig_primary), "\n")

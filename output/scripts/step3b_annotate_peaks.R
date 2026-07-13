## Step 3b: annotate consensus peaks to genomic context + nearest gene (hg38, EnsDb v86)
.libPaths(c("/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/output/Rlib", .libPaths()))
suppressPackageStartupMessages({library(GenomicRanges); library(GenomeInfoDb); library(ensembldb); library(EnsDb.Hsapiens.v86); library(data.table)})
base <- "/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/output"; cache <- file.path(base,"cache")
edb <- EnsDb.Hsapiens.v86
gall <- genes(edb); seqlevelsStyle(gall) <- "UCSC"; gall <- keepStandardChromosomes(gall, pruning.mode="coarse")
gpc  <- gall[gall$gene_biotype=="protein_coding"]
# TSS + promoter windows (protein-coding)
tss  <- resize(gpc, width=1, fix="start")                    # strand-aware TSS
prom <- promoters(gpc, upstream=2000, downstream=200)
# gene bodies (protein-coding, for genic vs distal)
body <- gpc
exons_gr <- reduce(exons(edb)); seqlevelsStyle(exons_gr) <- "UCSC"; exons_gr <- keepStandardChromosomes(exons_gr, pruning.mode="coarse")

annotate_peaks <- function(cons){
  names(cons) <- paste0(seqnames(cons),"-",start(cons),"-",end(cons))
  center <- resize(cons, width=1, fix="center")
  dt <- data.table(peak=names(cons), chr=as.character(seqnames(cons)),
                   start=start(cons), end=end(cons), width=width(cons))
  # promoter overlap (peak overlaps any promoter window)
  dt$is_promoter <- overlapsAny(cons, prom)
  dt$is_genic    <- overlapsAny(cons, body)
  dt$is_exonic   <- overlapsAny(cons, exons_gr)
  # nearest gene to peak center + distance
  nn <- distanceToNearest(center, tss, ignore.strand=TRUE)
  dt$nearest_gene_id <- NA_character_; dt$nearest_gene <- NA_character_; dt$dist_to_tss <- NA_integer_
  qi <- queryHits(nn); si <- subjectHits(nn)
  dt$nearest_gene_id[qi] <- gpc$gene_id[si]
  dt$nearest_gene[qi]    <- gpc$symbol[si]
  dt$dist_to_tss[qi]     <- mcols(nn)$distance
  # genomic-context class (priority: promoter > genic(intragenic) > distal)
  dt[, context := fifelse(is_promoter, "promoter",
                    fifelse(is_genic, "genic", "distal"))]
  dt[]
}

for(mk in c("H3K27ac","H3K27me3")){
  cons <- readRDS(file.path(cache, paste0("consensus_",mk,"_raw.rds")))
  a <- annotate_peaks(cons)
  fwrite(a, file.path(cache, paste0("peakann_",mk,".csv.gz")))
  cat(sprintf("\n%s: %d peaks | context: ", mk, nrow(a)))
  print(a[, .N, by=context][order(-N)])
  cat("  median dist_to_tss (distal):", a[context=="distal", median(dist_to_tss)], "bp\n")
}
cat("DONE peak annotation\n")

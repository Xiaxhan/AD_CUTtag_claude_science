## Step 1: QC & valid-cell filtering per library (all 12), with mislabel corrections
suppressPackageStartupMessages({library(Seurat); library(Signac); library(Matrix); library(GenomicRanges); library(data.table)})
inbase <- "/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/input"
cache  <- "/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/output/cache"
dir.create(cache, showWarnings=FALSE, recursive=TRUE)

dx <- c(Br5090="AD",Br5106="AD",Br5120="AD",Br5196="AD",Br5620="AD",
        Br2700="CTL",Br5417="CTL",Br5764="CTL",Br6492="CTL",Br8179="CTL")

libs <- data.frame(
 id=c("CBL_H3K27ac_0506","CBL_H3K27ac_0513_02","CBL_H3K27ac_0513_03",
      "CBL_H3K27me3_0506","CBL_H3K27me3_0513_02","CBL_H3K27me3_0513_03",
      "HIP_H3K27ac_LIB1","HIP_H3K27ac_LIB2","HIP_H3K27ac_LIB3",
      "HIP_H3K27me3_LIB1","HIP_H3K27me3_LIB2","HIP_H3K27me3_LIB3"),
 region=c("CBL","CBL","CBL","CBL","CBL","CBL","HIP","HIP","HIP","HIP","HIP","HIP"),
 mark=c("H3K27ac","H3K27ac","H3K27ac","H3K27me3","H3K27me3","H3K27me3",
        "H3K27ac","H3K27ac","H3K27ac","H3K27me3","H3K27me3","H3K27me3"),
 path=c(
  "CBL/H3K27ac/0506_CBL01/seurat/GAD_0506_CBL_H3K27ac_Br2700_Br5090_Br5120_Br5764.rawSeurat.RDS",
  "CBL/H3K27ac/0513_CBL02/seurat/GAD_CBL_0513BC02_H3K27ac_Br5106_Br5196_Br5417.rawSeurat.RDS",
  "CBL/H3K27ac/0513_CBL03/seurat/GAD_CBL_0513BC03_H3K27ac_Br5620_Br6492_Br8179.rawSeurat.RDS",
  "CBL/H3K27me3/0506_CBL01/seurat/GAD_0506_CBL_H3K27me3_Br2700_Br5090_Br5120_Br5764.rawSeurat.RDS",
  "CBL/H3K27me3/0513_CBL02/seurat/GAD_CBL_0513BC02_H3K27me3_Br5620_Br6492_Br8179.rawSeurat.RDS",
  "CBL/H3K27me3/0513_CBL03/seurat/GAD_CBL_0513BC03_H3K27me3_Br5106_Br5196_Br5417.rawSeurat.RDS",
  "HIP/H3K27ac/2025-06-05-LIB1/seurat/GAD_0605_HIP_H3K27ac_Br5417_Br5196_Br5620.rawSeurat.RDS",
  "HIP/H3K27ac/2025-06-05-LIB2/seurat/GAD_HIP_0605Lib2_H3K27ac_Br8179_Br5106_Br6492.rawSeurat.RDS",
  "HIP/H3K27ac/2025-06-09-LIB3/seurat/GAD_HIP_0609_H3K27ac_Br2700_Br5090_Br5120_Br5764.rawSeurat.RDS",
  "HIP/H3K27me3/2025-06-05-LIB1/seurat/GAD_0605_HIP_H3K27me3_Br5417_Br5196_Br5620.rawSeurat.RDS",
  "HIP/H3K27me3/2025-06-05-LIB2/seurat/GAD_0605LIB2_HIP_H3K27me3_Br8179_Br5106_Br6492.rawSeurat.RDS",
  "HIP/H3K27me3/2025-06-09-LIB3/seurat/GAD_0609_HIP_H3K27me3_Br2700_Br5090_Br5120_Br5764.rawSeurat.RDS"),
 stringsAsFactors=FALSE)

eff_region <- function(donor, region, mark, libid){
  if(donor=="Br5106" && libid=="CBL_H3K27ac_0513_02") return("HIP")
  if(donor=="Br5106" && libid=="CBL_H3K27me3_0513_03") return("HIP")
  if(donor=="Br5090" && libid=="HIP_H3K27ac_LIB3") return("CBL")
  if(donor=="Br5090" && libid=="HIP_H3K27me3_LIB3") return("CBL")
  if(donor=="Br5196" && libid=="CBL_H3K27me3_0513_03") return("DROP")
  region
}

TH <- list(minRNAcount=500, minRNAfeat=300, minPeakcount=300, maxMT=5)

std_chr <- paste0("chr", c(1:22,"X","Y"))
peaknames_to_gr <- function(pn){
  m <- regmatches(pn, regexec("^(.*)-([0-9]+)-([0-9]+)$", pn))
  ok <- lengths(m)==4
  m <- m[ok]
  gr <- GRanges(sapply(m,`[`,2), IRanges(as.integer(sapply(m,`[`,3)), as.integer(sapply(m,`[`,4))))
  gr[as.character(seqnames(gr)) %in% std_chr]
}

qc_rows <- list()
for(i in seq_len(nrow(libs))){
  L <- libs[i,]
  cat(sprintf("[%d/%d] %s ...\n", i, nrow(libs), L$id)); flush.console()
  so <- readRDS(file.path(inbase, L$path))
  DefaultAssay(so) <- "RNA"
  if("regemt" %in% Assays(so)){ so[["regemt"]] <- NULL; gc() }
  md <- so@meta.data
  n_raw <- nrow(md)
  keep_v <- md$info=="valid" & !(md$sample %in% c("multiplet","other")) & md$sample %in% names(dx)
  n_valid <- sum(keep_v, na.rm=TRUE)
  rna <- so[["RNA"]]$counts
  mt <- grep("^MT-", rownames(rna), value=TRUE)
  if(length(mt)){
    pmt <- Matrix::colSums(rna[mt,,drop=FALSE]) / pmax(Matrix::colSums(rna),1) * 100
    names(pmt) <- colnames(rna)
    md$percent.mt <- as.numeric(pmt[rownames(md)])
    mt_ok <- !is.na(md$percent.mt) & md$percent.mt < TH$maxMT
  } else {
    md$percent.mt <- NA_real_       # no MT genes in this reference; filter not applicable
    mt_ok <- TRUE
  }
  pass <- keep_v &
          !is.na(md$nCount_RNA)   & md$nCount_RNA   >= TH$minRNAcount &
          !is.na(md$nFeature_RNA) & md$nFeature_RNA >= TH$minRNAfeat &
          !is.na(md$nCount_peaks) & md$nCount_peaks >= TH$minPeakcount &
          mt_ok
  pass[is.na(pass)] <- FALSE
  bc <- rownames(md)[pass]
  n_pass <- length(bc)
  don <- md[bc,"sample"]
  effr <- mapply(eff_region, don, L$region, L$mark, L$id)
  drop_cortex <- effr=="DROP"
  bc <- bc[!drop_cortex]; don <- don[!drop_cortex]; effr <- effr[!drop_cortex]
  n_after_drop <- length(bc)
  dc <- table(don)
  rna_pass <- rna[, bc, drop=FALSE]
  metacols <- intersect(c("sample","nCount_RNA","nFeature_RNA","nCount_peaks","nFeature_peaks","FIE","FIP","PC.enrichment"), colnames(md))
  slim <- md[bc, metacols, drop=FALSE]
  slim$percent.mt <- md[bc,"percent.mt"]
  slim$donor <- don; slim$eff_region <- effr; slim$mark <- L$mark
  slim$library <- L$id; slim$dx <- dx[don]
  rownames(slim) <- paste0(L$id, "#", bc)
  colnames(rna_pass) <- rownames(slim)
  gr <- peaknames_to_gr(rownames(so[["peaks"]]))
  saveRDS(list(rna=rna_pass, meta=slim, peaks=gr, mark=L$mark, region=L$region, id=L$id),
          file.path(cache, paste0(L$id, "_qc.rds")))
  qc_rows[[L$id]] <- data.frame(library=L$id, region=L$region, mark=L$mark,
     n_raw=n_raw, n_valid_singlet=n_valid, n_qc_pass=n_pass, n_after_cortex_drop=n_after_drop,
     n_peaks_macs2=length(gr),
     donor_counts=paste(sprintf("%s=%d",names(dc),as.integer(dc)), collapse=";"),
     stringsAsFactors=FALSE)
  rm(so, rna, rna_pass, md, slim, gr); gc()
}
qc <- do.call(rbind, qc_rows)
write.csv(qc, "/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/output/qc_summary_table.csv", row.names=FALSE)
cat("DONE\n"); print(qc[,c("library","n_raw","n_valid_singlet","n_qc_pass","n_after_cortex_drop","n_peaks_macs2")])

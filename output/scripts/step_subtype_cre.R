## Per-subtype CRE re-quantification + edgeR DE (HIP astro/micro states)
.libPaths(c("/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/output/Rlib", .libPaths()))
suppressPackageStartupMessages({library(Signac); library(GenomicRanges); library(Matrix); library(data.table); library(edgeR)})
# state maps (corrected: mito-high clusters excluded consistently)
.sm <- readRDS("/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/output/cache/glial_state_maps.rds")
base <- "/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/output"
cache <- file.path(base,"cache")
inbase <- "/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/input"

ga <- readRDS(file.path(cache,"glial_subtype_assignments.rds"))
ga <- ga[!grepl("excl",state) & !is.na(state_collapsed)]   # use corrected collapsed states
ga[, library:=sub("#.*","",cell_id)]
ga[, mark:=ifelse(grepl("H3K27ac",library),"H3K27ac","H3K27me3")]
ga[, sc:=state_collapsed]
ga[, grp:=paste(donor, celltype, sc, sep="|")]
setkey(ga, cell_id)

FRELMAP <- c(
   CBL_H3K27ac_0513_02="CBL/H3K27ac/0513_CBL02", CBL_H3K27me3_0513_03="CBL/H3K27me3/0513_CBL03",
   HIP_H3K27ac_LIB1="HIP/H3K27ac/2025-06-05-LIB1", HIP_H3K27ac_LIB2="HIP/H3K27ac/2025-06-05-LIB2",
   HIP_H3K27ac_LIB3="HIP/H3K27ac/2025-06-09-LIB3", HIP_H3K27me3_LIB1="HIP/H3K27me3/2025-06-05-LIB1",
   HIP_H3K27me3_LIB2="HIP/H3K27me3/2025-06-05-LIB2", HIP_H3K27me3_LIB3="HIP/H3K27me3/2025-06-09-LIB3")
frag_path <- function(id) file.path(inbase, FRELMAP[id], "DNA/fragments.tsv.gz")
consensus <- list(
  H3K27ac  = readRDS(file.path(cache,"consensus_H3K27ac_raw.rds")),
  H3K27me3 = readRDS(file.path(cache,"consensus_H3K27me3_raw.rds")))
for(mk in names(consensus)) names(consensus[[mk]]) <- paste0(seqnames(consensus[[mk]]),"-",start(consensus[[mk]]),"-",end(consensus[[mk]]))

pb_list <- list()
for(id in names(FRELMAP)){
  mk <- ifelse(grepl("H3K27ac",id),"H3K27ac","H3K27me3")
  sub <- ga[library==id]; if(nrow(sub)==0) next
  bc <- sub("^.*#","",sub$cell_id); names(bc) <- sub$cell_id
  cons <- consensus[[mk]]
  fp <- frag_path(id)
  cat(sprintf("[%s] mark=%s  %d glial cells\n", id, mk, length(bc))); flush.console()
  fr <- CreateFragmentObject(path=fp, cells=unname(bc), verbose=FALSE, validate.fragments=FALSE)
  fm <- FeatureMatrix(fragments=fr, features=cons, cells=unname(bc), verbose=FALSE)  # peaks x bc
  # map bc -> global cell_id -> group
  cid <- sub$cell_id[match(colnames(fm), bc)]
  grp <- sub$grp[match(cid, sub$cell_id)]
  ug <- sort(unique(grp)); gi <- match(grp, ug)
  G <- sparseMatrix(i=seq_along(grp), j=gi, x=1, dims=c(length(grp), length(ug)))
  pb <- fm %*% G; colnames(pb) <- paste(ug, mk, sep="|")
  pb_list[[id]] <- pb
  rm(fm, fr, pb, G); gc()
}
# combine across libraries (same peak rows per mark); need per-mark separate matrices
save_pb <- list()
for(mk in c("H3K27ac","H3K27me3")){
  ids <- names(pb_list)[grepl(mk, names(pb_list))]
  mats <- pb_list[ids]
  # union columns (groups), same rows
  allcols <- unique(unlist(lapply(mats, colnames)))
  rn <- rownames(mats[[1]])
  M <- matrix(0, nrow=length(rn), ncol=length(allcols), dimnames=list(rn, allcols))
  for(m in mats){ M[, colnames(m)] <- M[, colnames(m)] + as.matrix(m) }
  save_pb[[mk]] <- M
  cat(mk, "combined pb:", nrow(M), "peaks x", ncol(M), "donor-celltype-state groups\n")
}
saveRDS(save_pb, file.path(cache,"subtype_pb.rds"))
cat("DONE requant\n")

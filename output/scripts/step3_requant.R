## Step 3: re-quantify fragments over consensus peaks + aggregate to pseudobulk
## Per library: FeatureMatrix(consensus peaks x QC-pass cells) -> donor x celltype pseudobulk sums.
suppressPackageStartupMessages({library(Signac); library(GenomicRanges); library(Matrix); library(data.table)})
base <- "/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/output"
cache <- file.path(base,"cache"); fragdir <- file.path(base,"frag_index")
pbdir <- file.path(base,"pseudobulk"); dir.create(pbdir, showWarnings=FALSE)
inbase <- "/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/input"

# cell annotations (major_type per global cell_id)
ann <- fread(file.path(base,"cell_annotations.csv.gz"))
setkey(ann, cell_id)

# fragment paths: indexed me3 ones in fragdir; ac ones use original input .tbi
FRELMAP <- c(
   CBL_H3K27ac_0506="CBL/H3K27ac/0506_CBL01", CBL_H3K27ac_0513_02="CBL/H3K27ac/0513_CBL02",
   CBL_H3K27ac_0513_03="CBL/H3K27ac/0513_CBL03", CBL_H3K27me3_0506="CBL/H3K27me3/0506_CBL01",
   CBL_H3K27me3_0513_02="CBL/H3K27me3/0513_CBL02", CBL_H3K27me3_0513_03="CBL/H3K27me3/0513_CBL03",
   HIP_H3K27ac_LIB1="HIP/H3K27ac/2025-06-05-LIB1", HIP_H3K27ac_LIB2="HIP/H3K27ac/2025-06-05-LIB2",
   HIP_H3K27ac_LIB3="HIP/H3K27ac/2025-06-09-LIB3", HIP_H3K27me3_LIB1="HIP/H3K27me3/2025-06-05-LIB1",
   HIP_H3K27me3_LIB2="HIP/H3K27me3/2025-06-05-LIB2", HIP_H3K27me3_LIB3="HIP/H3K27me3/2025-06-09-LIB3")
frag_path <- function(id) file.path(inbase, FRELMAP[id], "DNA/fragments.tsv.gz")

consensus <- list(
  H3K27ac  = readRDS(file.path(cache,"consensus_H3K27ac_raw.rds")),
  H3K27me3 = readRDS(file.path(cache,"consensus_H3K27me3_raw.rds")))
# name peaks chr-start-end
for(mk in names(consensus)) names(consensus[[mk]]) <- paste0(seqnames(consensus[[mk]]),"-",start(consensus[[mk]]),"-",end(consensus[[mk]]))

qcfiles <- list.files(cache,"_qc.rds$",full.names=TRUE)
peak_totals <- list()   # per mark: accumulate total fragment count per consensus peak
for(f in qcfiles){
  x <- readRDS(f); id <- x$id; mk <- x$mark
  cellids <- colnames(x$rna); bc <- sub("^.*#","",cellids)
  names(cellids) <- bc
  cons <- consensus[[mk]]
  fp <- frag_path(id)
  cat(sprintf("[%s] mark=%s  %d QC cells  frag=%s\n", id, mk, length(bc), basename(fp))); flush.console()
  fr <- CreateFragmentObject(path=fp, cells=bc, verbose=FALSE, validate.fragments=FALSE)
  fm <- FeatureMatrix(fragments=fr, features=cons, cells=bc, verbose=FALSE)  # peaks x cells (raw bc)
  # relabel cols to global cell_id
  colnames(fm) <- cellids[colnames(fm)]
  # peak totals for atlas
  pt <- Matrix::rowSums(fm)
  peak_totals[[mk]] <- if(is.null(peak_totals[[mk]])) pt else peak_totals[[mk]] + pt[names(peak_totals[[mk]])]
  # aggregate to donor x celltype pseudobulk (sum fragments)
  a <- ann[colnames(fm)]
  grp <- paste(a$donor, a$major_type, sep="|")
  # build group-sum via sparse matrix multiply
  ug <- sort(unique(grp)); gi <- match(grp, ug)
  G <- sparseMatrix(i=seq_along(grp), j=gi, x=1, dims=c(length(grp), length(ug)))
  pb <- fm %*% G                      # peaks x groups
  colnames(pb) <- ug
  # metadata per group
  meta <- data.table(group=ug, donor=sub("\\|.*","",ug), major_type=sub(".*\\|","",ug),
                     region=x$region, mark=mk, library=id)
  meta$dx <- ann[match(meta$donor, ann$donor), dx]
  meta$ncells <- as.integer(table(grp)[ug])
  saveRDS(list(pb=pb, meta=meta, id=id, mark=mk, region=x$region), file.path(pbdir, paste0(id,"_pb.rds")))
  cat("   groups:", length(ug), " peaks:", nrow(fm), "\n")
  rm(x, fm, fr, pb, G); gc()
}
saveRDS(peak_totals, file.path(cache,"peak_totals_by_mark.rds"))
cat("DONE requant\n")

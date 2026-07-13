# Methods, in plain English

*A step-by-step account of how the analysis was done, written to be readable without specialist background. Software versions are given the first time each tool is mentioned. hg38 genome build throughout.*

## What the data are
Each donor's brain tissue was processed with **Droplet Paired-Tag**, a method that reads out, from the very same individual nuclei, both which genes are switched on (RNA) and where two chemical "tags" sit on the DNA-packaging proteins (histones). The two tags have opposite meanings: **H3K27ac** marks regulatory switches that are *active*, and **H3K27me3** marks regions that are *shut down*. We had 10 donors (5 with Alzheimer's, 5 controls), two brain regions each (hippocampus, which is hit hard in Alzheimer's, and cerebellum, which usually is not), and both tags — delivered as 12 sequencing "libraries", each mixing 3–4 donors together.

## Step 0 — Fixing mislabeled samples
A file that came with the data flagged five samples as physically mislabeled in the wet lab. We worked out which nuclei each note referred to (by matching the donors present in each library) and, after confirming the choices with the project lead, moved two donors' material to its correct region and discarded one sample that turned out to be from a different brain area entirely (cortex, outside this study). This left a clean set where every donor is counted once per region and mark.

## Step 1 — Keeping only real cells
Raw droplet data contains many "empty" droplets. We kept only droplets confidently assigned to a single donor and carrying enough signal to be a real nucleus (at least 500 RNA molecules, 300 distinct genes, and 300 histone-tag fragments). A common mitochondrial-content filter did not apply here because the gene reference contained no mitochondrial genes. This left **147,228 nuclei**. Software: **Seurat 5.5.1** and **Signac 1.17.1** (R 4.5.3).

## Step 2 — Discovering the cell types from scratch
Rather than copying labels from elsewhere, we let the cell types emerge from the RNA. For each region we combined the nuclei, standardized the measurements, reduced them to their main axes of variation (principal component analysis), and used **Harmony 1.2** to remove technical differences between libraries so that biology, not batch, drove the grouping. We then grouped nuclei (Louvain clustering) and labeled each group by comparing it to textbook marker-gene lists for known brain cell types. Hippocampus gave 8 major cell types; cerebellum gave 9 (dominated, as expected, by granule neurons).

## Step 3 — Building one common ruler for the regulatory elements (the CRE atlas)
The regulatory "peaks" that came with each library were called separately and did not line up across donors (barely a hundred were shared), so they could not be compared. We fixed this by building **one consensus set of regulatory elements per tag** — merging every library's peaks into a common list (859,854 active-tag elements, 567,404 repressive-tag elements) — and then **re-counting every library's raw fragments against that common list** so all donors are measured on the same ruler. Peak calling used **MACS2 2.2**; re-counting used Signac; genomic annotation (which gene each element sits near) used **EnsDb.Hsapiens.v86** with hg38 gene models. Counts were summed within each donor × cell-type × region group ("pseudobulk"). Software for the merge: **GenomicRanges 1.60**.

## Step 4 — Labeling each element's regulatory state
The two tags sit on the *same* spot of the histone protein, so a single molecule cannot carry both at once. When we see both at one element in the pooled data, it means different cells in the population are doing different things — a *transitional/mixed* state, not the classic "poised" state (which would require a different tag we did not measure). We labeled each element, in each cell type, as **active** (active tag only), **repressed** (repressive tag only), **mixed/transitional** (both), or **unclear** (neither), using a 1-CPM signal threshold, then summarized to a genome-wide state. Notably, the Alzheimer's-risk gene APOE came out "cell-type-switched" — active in some cell types and repressed in others — matching what is known about its regulation.

## Step 5 — Testing for Alzheimer's-vs-control differences
This is the core statistical test. We compared AD to control **using each donor as one data point** (never individual cells — that would fake a much larger sample), separately for each cell type, region and tag. The test was **edgeR 4.x** quasi-likelihood, a standard tool for count data. We did *not* correct for sequencing batch, because batch was almost perfectly tangled with diagnosis in this cohort; including it invented false hits that disappeared once removed. Thresholds for calling a difference "significant": strict = 5% false-discovery rate and at least 1.4-fold change; relaxed = 10% and 1.2-fold.

**Result:** at this sample size, essentially nothing reached significance element-by-element (one lone element genome-wide). This is an honest limitation of 8–9 donors against hundreds of thousands of tests, not a data-quality problem.

## Step 6 — Recovering signal by pooling (gene-level and pathway-level)
Because single elements were underpowered, we pooled related signals two ways, as the project lead requested. First, we summed elements up to their nearest gene and re-tested. Second, we ranked all genes by their AD-vs-control trend and asked whether whole **biological pathways** moved coherently, using **fgsea 1.34** (a fast pathway-enrichment method) against gene-set collections from **Enrichr** (Gene Ontology, KEGG, Hallmark, WikiPathways). Pathways revealed clear, coordinated shifts (431 significant) that single elements had hidden.

## Step 7 — Linking regulatory elements to the genes they control
For every element we measured whether its tag signal rose and fell **together with** the nearby gene's expression across all the donor × cell-type groups (Spearman correlation, within 100,000 bases). Active-tag elements correlated *positively* with expression (they turn genes on) and repressive-tag elements *negatively* (they turn genes off) — exactly as they should, which is a strong internal check that the whole pipeline is measuring real regulation.

## Step 8 — Building the model and picking validation targets
We combined everything — where remodeling is strongest, in which direction, which genes are linked, and which are known Alzheimer's-risk genes — into a ranked shortlist of candidate regulatory elements to test experimentally, both overall and separately for each region × cell-type context. We removed genes that are unreliable in small post-mortem cohorts (sex chromosomes, highly copy-number-variable gene families) before ranking. Each shortlisted target names a specific element (its genomic coordinates), its cell type, and the predicted direction of change.

## Step 9 — Glial cell sub-states
For the three main "support" cell types (astrocytes, microglia, oligodendrocytes) we asked whether Alzheimer's shifts the balance of their internal activity states. There were enough cells to do this reliably only in the hippocampus. Microglia split cleanly into resting, activated and inflammatory states, but **the proportions did not differ significantly between Alzheimer's and controls** — telling us the disease signal lives in the fine regulatory tuning *within* each cell type, not in a wholesale change of cell-state mix.

## A note on honesty and reproducibility
Every figure and table is saved as a versioned file; every methodological choice is logged with its rationale, pros and cons in `DECISIONS.md`; and a dated running log is in `NOTEBOOK.md`. Where the data could not support a conclusion (single-element differences, cerebellar microglia numbers, glial-proportion shifts), that is stated plainly rather than glossed over. The strongest findings are pooled and pathway-level, and should be treated as well-supported hypotheses to be confirmed in an independent, larger cohort.

---

# Technical appendix — exact parameters & settings

*For readers who want the precise values. Software versions are given the first time each tool appears. Genome build hg38 throughout; all analysis in R 4.5.3.*

## Software
| Tool | Version | Use |
|---|---|---|
| R | 4.5.3 | analysis environment |
| Seurat | 5.5.1 | single-cell object handling, normalization, clustering |
| Signac | 1.17.1 | chromatin assay, FeatureMatrix requantification, fragment objects |
| Harmony | 1.2 | batch integration over library |
| MACS2 | 2.2 | peak calls (as delivered with the data) |
| GenomicRanges | 1.60 | consensus peak construction (union + reduce) |
| EnsDb.Hsapiens.v86 | 2.99.0 | hg38 gene models / peak annotation |
| edgeR | 4.x | pseudobulk differential testing (quasi-likelihood) |
| fgsea | 1.34 | rank-based pathway enrichment |
| htslib (tabix/bgzip) | 1.23.1 | fragment file indexing |

## Cohort (post-correction)
- HIP H3K27ac n=9 (4 AD / 5 CTL); HIP H3K27me3 n=9 (4/5); CBL H3K27ac n=9 (4/5); CBL H3K27me3 n=8 (3 AD / 5 CTL, weakest arm).
- Mislabel corrections: Br5106 CBL libraries → HIP; Br5090 HIP-LIB3 material → CBL; Br5196 CBL-H3K27me3 (cortex) dropped.

## QC thresholds
- Kept demultiplexed **singlets** only (`info == valid`); dropped multiplets and background droplets.
- nCount_RNA ≥ 500; nFeature_RNA ≥ 300; nCount_peaks ≥ 300.
- percent.mt: **not applied** (no mitochondrial genes in the delivered reference → percent.mt = NA in all libraries).
- Result: **147,228 QC-pass nuclei** across 12 libraries.

## Clustering & annotation
- Per region; both marks' nuclei together (same nuclei in Paired-Tag).
- LogNormalize (scale factor 10⁴) → FindVariableFeatures (vst, 2,000 HVG) → ScaleData → RunPCA (30 PCs) → RunHarmony(group.by.vars = "library") → FindNeighbors (30 Harmony dims) → FindClusters (Louvain, **resolution 0.6**) → RunUMAP.
- Annotation: AddModuleScore over region-appropriate canonical marker panels; z-scored cluster means; argmax → subtype; collapse to major type. HIP 8 major types; CBL 9.

## Consensus CRE matrix
- Per mark: union of all per-library MACS2 peak GRanges → `reduce()` → standard chromosomes only.
- **H3K27ac: 859,854 peaks** (median width 677 bp); **H3K27me3: 567,404 peaks** (median 2,022 bp).
- Requantification: Signac `FeatureMatrix` per library over the mark's consensus, QC-pass barcodes only; columns relabeled to global cell IDs; summed to donor × cell-type × region pseudobulk (**152 H3K27ac / 144 H3K27me3 groups** after region correction).
- Peak annotation: context priority promoter (−2 kb / +200 bp) > genic > distal; nearest protein-coding gene by strand-aware TSS `distanceToNearest`.

## Differential CRE (edgeR)
- Strata: cell-type × region × mark with ≥3 AD and ≥3 CTL donors and ≥25 cells/donor group (**29 analyzable**).
- Pipeline: DGEList on donor pseudobulk → TMM normalization → `filterByExpr` → `estimateDisp` → `glmQLFit` / `glmQLFTest`, design **`~diagnosis`** (batch excluded — near-confounded with diagnosis).
- Significance: **primary** FDR < 0.05 & |log2FC| ≥ 0.5; **relaxed** FDR < 0.1 & |log2FC| ≥ 0.25.
- Outcome: 1 relaxed-significant CRE genome-wide (IMPAD1 distal H3K27ac, HIP.InN, log2FC −1.38, FDR 0.02).

## Pooled-power analyses
- Gene-level DE: CREs summed to nearest gene, same edgeR pipeline.
- GSEA: fgsea on genes ranked by signed −log10(p) × sign(log2FC); **minSize 10, maxSize 500, nPermSimple 10,000**; gene sets = Enrichr GMTs (GO Biological Process 2021, KEGG 2021 Human, MSigDB Hallmark 2020, WikiPathways 2019). Significance padj < 0.1 → **431 enrichments**.

## CRE–gene associations
- Per CRE with nearest gene ≤ 100 kb: Spearman correlation of CRE log-CPM vs gene expression across matched pseudobulk groups; BH FDR. Filter: ≥5 nonzero CRE groups and sd > 0 both sides.
- H3K27ac: 163,756 activating / 3,054 repressive significant links. H3K27me3: 40,449 repressive / 11,934 activating.

## CRE regulatory-state classification
- 616,374 unified loci = reduce(union of both marks' consensus).
- Per common cell type (threshold **1 CPM**): active = ac ≥ 1 & me3 < 1; repressed = me3 ≥ 1 & ac < 1; dual = both ≥ 1; unclear = both < 1.
- Global state = dominant per-cell-type state; active-in-some + repressed-in-others → cell-type-switched.

## Glial subtype sub-analysis (HIP only)
- Per glial major type: subset → NormalizeData → 2,000 HVG → ScaleData → RunPCA (20) → FindNeighbors → FindClusters (**resolution 0.3**) → AddModuleScore on state panels.
- AD-vs-CTL reactive-score shift: per-donor Wilcoxon (n=9). Microglia p=0.73; astrocytes p=0.56; oligodendrocytes p=0.29 (none significant).

## Target prioritization
- Priority score = (−gene_AD_log2FC) × (1 + best_rho) × (1.5 if AD-GWAS gene else 1), per region × cell-type context.
- Filtered sex-chromosome and germline-CNV/hypervariable loci (USP17L, GSTM1/GSTT1, HLA, Ig/OR/LILR/KIR, UGT2B, C4, Y-genes) as post-mortem confounders.

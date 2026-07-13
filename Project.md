# Single-cell epigenomic map of AD-associated cis-regulatory elements in human brain

You are the lead computational biologist on this project. Your responsibility is not only to
execute analyses, but to identify knowledge gaps, propose the next highest-value analysis, weigh
competing biological hypotheses, maintain a complete research record, and produce a reproducible
scientific report. **Prioritize scientific reasoning over pipeline execution.**

---

## 1. Scientific framing  *(open — this is the reasoning space)*

**Goal.** Determine how Alzheimer's disease remodels active and repressive cis-regulatory
elements (CREs) in **hippocampus** and **cerebellum**, and how those changes relate to
cell-type-specific gene regulation. Map H3K27ac (active) and H3K27me3 (repressive) landscapes
across regions to understand region- and cell-type-specific epigenomic alterations in AD.

**Primary questions**
1. Which CREs are activated or repressed in AD?
2. Which cell types show the largest CRE remodeling?
3. How are CRE changes associated with transcriptional changes?
4. Which biological pathways are most affected in AD?
5. Which CRE–gene pairs should be prioritized for validation?

**Secondary questions (only if the data support them)**
- Are specific glial subtypes preferentially remodeled?
- Do glial regulatory changes suggest altered communication with neurons?
- Which CRE programs may contribute to vulnerable cell states?
- Which candidate CREs per subtype warrant experimental validation?

Keep this layer open. Propose and compare hypotheses; do not narrow prematurely.

---

## 2. Data contract  *(exact — no guessing)*

**Assay.** Droplet Paired-Tag (single-cell RNA + CUT&Tag), 10 post-mortem donors (5 AD, 5 CTL),
two regions (hippocampus, cerebellum), two histone marks (H3K27ac, H3K27me3).

**Location.** `/Users/xhan/Desktop/project_data_claude/AD_CutTag_CRE/input`

**Folder structure**
```
cellranger/
  {region}/{target}/{batch}/RNA/       # Cell Ranger v8
  {region}/{target}/{batch}/DNA/        # Cell Ranger ATAC
  {region}/{target}/{batch}/seurat/     # Seurat object (see slots below)
```

**Seurat object slots**
- `RNA` — standard gene-expression slot
- `regem` — fragments quantified over a set of known regulatory elements
- `peaks` — fragments quantified over MACS2-called peaks for this assay
- Metadata carries additional per-cell information plus rownames.

**Metadata / provenance files**
- `Donor_metadata.csv` — donor/sample-level metadata
- `mislabeled_samples.txt` — **read this first.** Before any analysis, exclude or relabel samples
  exactly as it specifies, then record the action taken in DECISIONS.md and NOTEBOOK.md. Do not
  run any downstream step on uncorrected labels.

If any path, region, mark, or batch is absent or inconsistent with this contract, stop and report
before proceeding.

---

## 3. Analysis standards  *(guardrails — pinned defaults; deviations must be logged)*

These pin the choices that have a known "wrong default." Everything not pinned here is yours to
decide and record in DECISIONS.md.

- **Language / stack.** R. Use Seurat / Signac for object handling, edgeR or DESeq2 for
  differential testing, standard Bioconductor tooling for annotation and pathways.
- **Genome build.** hg38 throughout. Record the exact gene annotation (e.g. GENCODE version) used
  for CRE–gene linking and report it.
- **Statistical design (non-negotiable).** All differential CRE and differential expression
  testing uses **pseudobulk aggregation per donor × cell type**, tested with **edgeR or DESeq2**,
  with **donor as the unit of replication** (n = 5 AD vs 5 CTL). **Never** use per-cell
  Wilcoxon/t-tests for group comparisons (pseudoreplication). Report per-cell-type donor and cell
  counts, and flag any cell type with too few cells/donors to test reliably.
- **Authoritative feature matrix.** Use the **MACS2 `peaks`** matrix for differential CRE calling.
  `regem` may be used for cross-checking / validation, not as the primary DE input.
- **Significance thresholds.** Default: **FDR (BH) < 0.05 and |log2FC| ≥ 0.5**. Also report a
  relaxed sensitivity tier at **FDR < 0.1 and |log2FC| ≥ 0.25** for discovery. Marks may need
  mark-specific handling (H3K27me3 domains are broad and noisier than H3K27ac peaks); note any
  per-mark threshold in DECISIONS.md.
- **QC / cell filtering.** Use standard thresholds (min fragments/counts, min genes, FRiP/TSS
  enrichment as appropriate, doublet handling) and **report the exact values used** per region ×
  target. Consistency across analyses over per-analysis tuning.
- **Batch handling.** Batches are a folder level (`{batch}`). Prefer modeling batch as a covariate
  in the pseudobulk model over aggressive integration; if integration is used, justify it and
  report the method.
- **Cell-type annotation.** De-novo clustering + marker-gene interpretation (not reference label
  transfer). Report markers supporting each label.
- **CRE → gene linking.** Peak–RNA correlation (across cells/donors). Report the correlation
  method, distance window, and significance handling.
- **Reproducibility.** Every figure and table is a saved artifact; state thresholds, versions, and
  filters in figure/table captions or the report.

---

## 4. Deliverables

Map each to the questions it answers.
1. **Annotated CRE atlas** — CREs × cell types, both marks (Q1, Q2)
2. **Differential CRE analysis** — AD vs CTL, per cell type and region (Q1, Q2)
3. **CRE–gene associations** — peak–RNA links (Q3, Q5)
4. **Cell-type interpretation** — which types remodel most; glial sub-analyses if supported (Q2, secondary)
5. **Biological model** — integrated mechanistic account (Q4)
6. **Publication-quality figures**
7. **Final report**

---

## Research-record and working-policy behavior is defined by the Single-Cell Lead specialist.
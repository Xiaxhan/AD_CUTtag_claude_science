# Critical evaluation — sharpening the AD CUT&Tag CRE discovery

*Companion to the four-axis CRE–gene prioritization. Scope: what the study's central claim should be, which deliverable best tests it, and where a demanding reviewer will push. Genome build hg38; donor is the unit of replication throughout (5 AD / 5 CTL, 8–9 per arm after label correction).*

---

## 1. The central biological question this study actually answers

The study was framed around five questions, of which Q1 ("which individual CREs are activated or repressed in AD") cannot be answered at this cohort size: single-element differential testing returns **one** CRE genome-wide at the relaxed threshold. That null is a power result, not a data-quality result — dispersions are tight (BCV 0.12–0.18) and it is 8–9 donors against hundreds of thousands of tests.

The question the data **do** answer, and the one that should anchor the paper, is:

> **Does AD remodel the cell-type-resolved *cis*-regulatory landscape in a coordinated, direction-consistent way that converges on neuroimmune and AD-risk gene programs — and can each convergence be reduced to a specific, testable CRE → gene → cell-type hypothesis?**

The evidence for "yes" is not any single element. It is (i) the pooled, direction-concordant pathway signal (431 GSEA enrichments where single-CRE testing found ~none); (ii) the internal validation that H3K27ac links track expression up and H3K27me3 links track it down genome-wide; and (iii) — the piece this analysis adds — that when a CRE's AD-associated H3K27ac change, its regulatory link to a gene, and that gene's AD-associated RNA change are examined *together*, they agree in direction in **54%** of well-linked models and cluster on AD-risk genes in the disease-vulnerable region. The contribution is the move from *cataloguing* differential signal to a *mechanistic, testable* CRE–gene model.

## 2. The most useful deliverable for testing the model

**The four-axis ranked CRE–gene table (`csv/cre_gene_model_ranked.csv`)** — not the pathway list, and not the original prioritized-targets table. Three reasons:

1. **It restores the missing axis.** The project's existing `prioritized_targets_*` tables score targets on chromatin remodeling × CRE–gene link × GWAS weight. Their `gene_log2FC` column is **chromatin signal aggregated to genes**, not RNA. In a Paired-Tag design that measures the transcriptome in the *same* nuclei, leaving the actual transcriptional consequence out of the prioritization is the single largest missed opportunity. This table adds a genuine donor-level RNA edgeR log2FC per stratum (`csv/rna_de_by_stratum.csv.gz`) as an independent axis.
2. **Direction-concordance is a stronger filter than any single axis.** A CRE that loses H3K27ac, is positively linked to its gene, *and* whose gene falls in RNA is a coherent regulatory hypothesis; any one of those alone is not. Concordance is the column a reviewer will most trust and the one that most cleanly nominates an experiment.
3. **Every row is already an experiment.** Each row names a cell type, a mark, a specific CRE interval, a target gene, an expected direction, and whether the gene carries AD genetic risk — i.e. exactly the specification for a CRISPRi/a enhancer perturbation or an allele-specific reporter assay.

## 3. The model to test first: *SORL1* in hippocampal microglia

Selected as the highest-scoring hippocampus AD-GWAS model in which all four axes agree (see `csv/top_AD_hippocampus_models.csv`; visualized in `locus_model_SORL1.png`):

| Axis | Value |
|---|---|
| Epigenomic remodeling | H3K27ac **loss**, log2FC −0.90, p = 0.023 at a genic enhancer **chr11:121,481,988–121,483,897** (+30.7 kb from the *SORL1* TSS) |
| CRE → gene link | **activating**, Spearman ρ = +0.335 (this element tracks *SORL1* expression up) |
| RNA expression change | *SORL1* **down** in AD microglia, log2FC −0.40 — **concordant** with enhancer loss |
| AD relevance | *SORL1* is a top-tier AD gene (endo-lysosomal sorting; both common-variant GWAS and highly-penetrant rare coding variants), and the change is in **microglia**, the cell type where AD heritability concentrates |

The signal is **cell-type-specific**: at neighbouring *SORL1* CREs, every other hippocampal and cerebellar cell type shows a gain or a weaker/non-significant change — only microglia lose this element. The fragment-level pileup (window-normalized to remove sequencing-depth differences) reproduces the depletion the edgeR model reports (−1.05 vs −0.90), so the track and the statistic agree.

**Testable prediction.** CRISPRi silencing (or deletion) of chr11:121,481,988–121,483,897 in a human microglial model (iMGL / HMC3) should lower *SORL1* expression; the AD-associated loss of this enhancer is predicted to reduce microglial *SORL1* and, downstream, endo-lysosomal / Aβ-clearance capacity. An allele-specific or reporter assay across AD-risk genotypes at the locus would test whether genetic risk acts through this element.

## 4. What a demanding reviewer will attack — and the honest answers

1. **"Your single-CRE effects are not genome-wide significant."** Correct, and stated on every claim. The −0.90 at *SORL1* has p = 0.023 but does **not** survive genome-wide FDR (the near-null Q1 result). The models are **hypothesis-generating**: the value is in the *convergence* of four independent axes on the same locus, not in any one FDR. This is the same regime as early single-cell AD transcriptomics — real, cell-type-specific, subtle per-feature.
2. **"The original headline was cerebellar microglia."** The prior single-axis framing rested on CBL microglia H3K27ac immune loss — the **thinnest** arm (only **41 CREs** passed `filterByExpr`), in a region classically spared in AD, with immune-biased Enrichr gene sets that can inflate immune themes. Prioritizing **hippocampus** (AD-vulnerable, well-powered: HIP microglia H3K27ac has 4 AD / 5 CTL donors, 3,812 cells) puts the testable model on the region where the disease biology and the statistical power both live.
3. **"Correlation, not causation."** All links are Spearman correlations across pseudobulk groups; the CRE–gene assignment is nearest-TSS within 100 kb. The deliverable is explicitly a *prioritization* for perturbation, not a causal claim — hence a single, concrete CRE interval per row.
4. **"RNA change is also not significant."** True — *SORL1* RNA log2FC −0.40 has p = 0.22; single-gene RNA DE is as power-limited as the chromatin. It enters the model as a **directional effect-size axis** that must *agree* with the chromatin, not as an independent significant finding. Concordance across two independently under-powered measurements is the signal.
5. **"Confounds."** The diagnosis-only edgeR model is the defensible choice given near-complete batch–diagnosis confounding (documented in DECISIONS.md #08). Sex-chromosome and germline-CNV/hypervariable loci are filtered from the model table. Post-mortem, cross-sectional design — associations only.

## 5. Deliverables produced here

- **`csv/cre_gene_model_ranked.csv`** — 58,714 four-axis CRE–gene models (columns: region, cell type, mark, gene, CRE peak, `cre_log2FC`/`cre_p`, `link_rho`/`link_class`, `rna_log2FC`/`rna_p`, `ad_gwas`, `concordant`, `score`). 31,677 concordant; 191 on AD-GWAS genes.
- **`csv/rna_de_by_stratum.csv.gz`** — the new axis: donor-level RNA edgeR AD-vs-CTL log2FC per region × cell type (17 strata, 186k gene-stratum rows).
- **`model_prioritization_scatter.png`** — the four axes on one plane (chromatin-predicted vs observed RNA), sized by link strength, AD-GWAS HIP genes labeled.
- **`locus_model_SORL1.png`** — the selected model as locus-level chromatin tracks.
- **`csv/top_AD_hippocampus_models.csv`** — the 15 highest-scoring hippocampus AD-GWAS concordant models.

# DECISIONS — AD CUT&Tag CRE project

Format: Decision #NN; Reason; Pros; Cons; Confidence.

---

## Decision #01 — Sample mislabel decode (Donor-Region-Batch → delivered folders)
**Context.** `mislabeled_samples.txt` uses wet-lab nomenclature that does not map 1:1 to the delivered `Region/Mark/Library` folders. Raw Seurat objects carry NO region/batch/target field (region & mark are folder-defined only). Decode was done by matching each entry's donor to its physical library membership (from `barcode_assignments.csv.gz` demultiplexing).
- `Br5106-CBL-Batch02 → hippocampus`  ⇒ Br5106 in CBL/H3K27ac/0513_CBL02 (Batch02=ac library)
- `Br5106-CBL-Batch03 → hippocampus`  ⇒ Br5106 in CBL/H3K27me3/0513_CBL03 (Batch03=me3 library)
- `Br5196-CBL-Batch03 → cortex`       ⇒ Br5196 in CBL/H3K27me3/0513_CBL03
- `Br5090-Hip-Batch02 → cerebellum`   ⇒ mapped to Br5090's delivered HIP LIB3 (only HIP library it is in)
- `Br5090-WM-Batch01 → hippocampus`   ⇒ UNMAPPABLE (no white-matter folder delivered); ignored
**Reason.** Donor library membership is the only anchor tying wet-lab names to delivered files.
**Confidence.** Medium-high for the three CBL entries; low for Br5090 entries (batch/region names not delivered).

## Decision #02 — Handling of mislabeled cells (USER-CONFIRMED)
**Choice.** Reassign, don't blanket-exclude:
- Br5106: two delivered CBL samples (ac 0513_CBL02, me3 0513_CBL03) reassigned to HIP. Br5106 already has genuine HIP LIB2 samples, so Br5106 will have TWO source libraries per mark in HIP → these are POOLED into one donor-level pseudobulk (donor stays the replication unit; no pseudoreplication).
- Br5090: `Hip-Batch02→cerebellum` applied at dissection level → Br5090's delivered HIP LIB3 (both marks) reassigned to CBL. Br5090 already has genuine CBL 0506_CBL01 samples → POOLED per mark in CBL. `WM-Batch01` entry ignored (undeliverable region).
- Br5196: CBL-me3 (0513_CBL03) is cortex → DROPPED (out of scope). Br5196 CBL-ac (0513_CBL02) kept as cerebellum.
**Reason.** User elected to retain data via reassignment rather than discard; cortex is out of the HIP/CBL scope.
**Pros.** Preserves donor n; keeps region labels biologically correct.
**Cons.** Two donors (Br5106, Br5090) become 2-library pools in their reassigned arm → must aggregate at donor level (already required by pseudobulk design). CBL-me3 arm drops to n=8 (3 AD / 5 CTL).
**Confidence.** High for the corrected design table (below); pooling handled by donor×celltype aggregation.

### Corrected per-arm design (donor as replication unit)
| Arm | n | AD | CTL |
|-----|---|----|----|
| HIP H3K27ac  | 9 | 4 (Br5106,Br5120,Br5196,Br5620) | 5 |
| HIP H3K27me3 | 9 | 4 (Br5106,Br5120,Br5196,Br5620) | 5 |
| CBL H3K27ac  | 9 | 4 (Br5090,Br5120,Br5196,Br5620) | 5 |
| CBL H3K27me3 | 8 | 3 (Br5090,Br5120,Br5620)          | 5 |
Pooled (2-library) donor-arms: CBL-ac & CBL-me3 Br5090; HIP-ac & HIP-me3 Br5106.
Full row-level mapping: output/sample_manifest_corrected.csv

## Decision #03 — CRE regulatory-state framework (2-mark, K27-aware) [USER-REQUESTED STEP]
**Context.** Only H3K27ac + H3K27me3 measured. Both modify the SAME residue (H3 Lys27), so they are chemically mutually exclusive on a single nucleosome tail (literature: <0.1% co-occurrence on the same tail). True "bivalency"/"poised" states are defined by H3K4me3+H3K27me3 (promoters) or H3K4me1+H3K27me3 (enhancers) — neither H3K4 mark is available here.
**Framework (per CRE x cell type).**
- Active            = H3K27ac+ / H3K27me3-
- Repressed (PcG)   = H3K27me3+ / H3K27ac-
- Dual/transitional = both marks called at the CRE (population-level) — NOT single-nucleosome bivalency; reflects cell/allele heterogeneity or an active<->repressed transition
- Unclear/quiescent = neither mark above threshold
**Reason.** Faithful to chromatin-state literature while not overclaiming bivalency from 2 antibodies on the same residue.
**Confidence.** High (framework); mark-call thresholds set empirically & logged separately.

## Decision #04 — Authoritative CRE matrix: consensus MACS2 peaks + re-quantify (USER-CONFIRMED)
**Context.** MACS2 `peaks` assay is library-specific (687k-863k peaks/library; only ~106 shared between two libraries) → NOT comparable across donors as delivered. RNA (36,601 genes) and regemt (1,638,269 REs) ARE identical across all libraries. Fragment files present for all 12 (bgzip); 5 H3K27me3 libraries lack .tbi (re-index needed).
**Choice.** Build consensus MACS2 peak set per mark (merge per-library peaks, standard chromosomes), re-quantify each library's fragments over consensus via Signac FeatureMatrix, restricted to QC-pass cells. Use as primary DE input. regemt reserved for cross-check.
**Reason.** Honors the pinned default (data-driven MACS2 peaks primary) while making features cross-comparable (Signac/ArchR standard).
**Pros.** Assay-specific, data-driven CRE definition; cross-donor comparable. **Cons.** Heaviest compute; runs in background.
**Confidence.** High.

## Decision #05 — QC thresholds (applied uniformly across all 12 libraries)
**Choice (per cell, on info==valid demultiplexed singlets; drop multiplet/other/background).**
- nCount_RNA >= 500, nFeature_RNA >= 300 (RNA modality usable for de-novo clustering)
- nCount_peaks >= 300 (CUT&Tag signal present)
- percent.mt < 5 — SPECIFIED but NOT APPLICABLE: this gene reference contains NO mitochondrial (MT-) genes, so percent.mt could not be computed and the filter was not applied in any library (percent.mt set to NA). Only the three thresholds above were operative.
**Reason.** Raw objects' `info==valid` flag is loose (median RNA counts 2-3 = near-empty droplets); distributions are clearly bimodal with a trough ~200-500 and the real-cell mode >1000 for both modalities. Thresholds sit above the trough. Uniform across libraries per contract (consistency > per-analysis tuning).

## Decision #06b — CORRECTION: pseudobulk region label bug (found & fixed Step 5)
**Issue.** First pb merge (168 ac / 161 me3 groups) tagged each pseudobulk group's `region` by its source library's NOMINAL folder region, not the mislabel-corrected assignment. This mis-placed Br5090's HIP-LIB3 cells (should be CBL) as HIP and Br5106's CBL-library cells (should be HIP) as CBL — creating spurious cross-region singleton strata (HIP.GranuleN, HIP.Purkinje, CBL.Ependymal) and a stray Br5090 HIP donor-group.
**Fix.** Rebuilt merge deriving region per (donor × source-library) from the authoritative csv/cell_annotations.csv.gz key (which was always correct). Corrected groups: 152 ac / 144 me3. Br5090 now CBL-only, Br5106 HIP-only. Corrected analyzable strata = 29 (HIP 4AD/5CTL, CBL 3-4AD/4-5CTL). Atlas (Deliverable 1) & CRE-state (Step 4) were regenerated from the corrected merge.
**Impact.** Small (a few hundred misassigned cells among ~150 groups) but removes spurious strata and restores the donor-as-unit design exactly per Decision #02. Confidence: High (verified Br5090/Br5106 single-region post-fix).

## Decision #08 — DE model: diagnosis-only (no batch covariate)
**Decision.** Primary differential-CRE model is edgeR quasi-likelihood `~diagnosis` (AD vs CTL), donor pseudobulk as replication unit, TMM norm, filterByExpr, per region×cell-type×mark stratum. Batch (sequencing library) is NOT included as a covariate. Strata require ≥3 AD & ≥3 CTL donors with ≥25 cells (29 analyzable). Significance: primary FDR<0.05 & |log2FC|≥0.5; relaxed FDR<0.1 & |log2FC|≥0.25.
**Reason.** Project.md prefers batch-as-covariate, but that presumes batch is estimable. Here each donor contributes to exactly one library per region×mark, libraries carry up to 4 levels for 8–9 donors, some levels hold a single donor, and library is near-confounded with diagnosis. Including it leaves almost no residual df AND manufactures spurious hits: a batch-adjusted run produced "significant" GSTM1/USP17L4/TMEM106B CREs (FDR~2e-5) that ALL vanish (minFDR 0.66–0.74, zero sig) under the diagnosis-only model — a textbook batch-confounding-as-signal artifact. Diagnosis-only is the sound choice for this design.
**Pros.** No false positives from an unstable/confounded covariate; honest power; standard for small-n pseudobulk. **Cons.** Cannot remove genuine library-level technical variance (mitigated by TMM + pseudobulk aggregation). Confidence: High.

## Decision #09 — Complementary pooled-power analyses for the null single-CRE result
**Decision.** Single-CRE DE is essentially null (1 CRE relaxed-sig genome-wide: distal H3K27ac near IMPAD1, HIP.InN, down in AD). Per user ruling, keep the honest single-CRE result and ADD pooled-power analyses: (a) gene-level CRE aggregation DE (sum CREs per nearest gene, ~19k genes/stratum), and (b) rank-based GSEA (fgsea) on gene-level signed −log10P rankings against MSigDB Hallmark / KEGG 2021 / GO-BP 2021 / WikiPathways gene sets (from Enrichr GMTs; org.Hs.eg.db/GO.db unavailable in sandbox).
**Reason.** n=8–9 donors vs 70k–440k tests has near-zero single-CRE power after FDR; concordant weak effects are only detectable when pooled. GSEA recovered 431 significant pathway-stratum enrichments where single-CRE found ~none.
**Pros.** Recovers interpretable biology without inflating single-CRE false positives. **Cons.** Gene/pathway-level, not individual-CRE resolution; gene sets are immune-blood-biased (Enrichr). Confidence: Medium-High.

## Decision #10 — Prioritized validation target ranking scheme
**Decision.** Rank candidate CRE-gene validation targets from the strongest, most coherent AD signal (CBL-microglia immune H3K27ac loss). Steps: (1) extract fgsea leading-edge genes from the significant immune gene sets in CBL.Micro H3K27ac (203 genes); (2) for each, take its CBL.Micro gene-level AD log2FC (loss direction) and its best activating H3K27ac CRE-gene link (max rho, from Deliverable 3); (3) priority score = (−gene_AD_log2FC) × (1 + best_rho) × (1.5 if AD-GWAS gene else 1); (4) report top 15 each with the specific CRE coordinate to validate.
**Reason.** Targets should combine a coherent disease signal (immune-enhancer loss), a concrete regulatory link (CRE→gene correlation gives the element to test), and directional AD effect. AD-GWAS weighting prioritizes translatable candidates.
**Pros.** Each target is a testable CRE+gene+cell-type+direction hypothesis. **Cons.** Ranking depends on the immune-loss finding being real (cerebellum-specific, small strata); score weights are heuristic; gene-level effects individually non-significant. Confidence: Medium (hypothesis-generating).

## Decision #11 — Two-axis model framing + artifact gene filtering (user-directed revision)
**Decision.** (a) Frame the immune-enhancer findings cautiously: reduced H3K27ac at inflammatory programs = region-specific immune dampening / altered responsiveness / possible compensation — NOT "silencing" or unqualified suppression. (b) Give hippocampus co-equal weight with cerebellum as an AD-relevant axis (HIP is the more AD-vulnerable region): HIP.ExN & HIP.Oligo carry the 2nd/3rd highest remodeling burden, and HIP.Astro shows reactive-astrocyte H3K27ac GAIN at NF-κB/TLR genes (CHI3L1/CCL2/TLR4). (c) Provide prioritized targets per region×cell-type context, not just one global list. (d) Filter sex-chromosome (Y-genes, XIST) and germline-CNV / hypervariable loci (USP17L, GSTM1/GSTT1, HLA, Ig/OR/LILR/KIR, UGT2B, C4) from target lists — these reflect donor sex imbalance / CNV, not AD regulation.
**Reason.** User feedback (three PIs' worth of caution): the original single-axis "microglial immune silencing, cerebellum" framing over-reached on direction and under-weighted the AD-vulnerable region. Sex/CNV genes are classic confounders in small post-mortem cohorts.
**Pros.** Biologically honest, region-balanced, translatable targets. **Cons.** More diffuse headline (two axes, cautious); some filtered genes (e.g. GSTM1) could carry real signal but can't be distinguished from CNV here. Confidence: High for the caution; Medium for which axis dominates.

## Decision #07 — CRE regulatory-state classification method
**Decision.** Classify CRE states from matched pseudobulk signal, NOT peak-set overlap. Steps: (1) unified CRE loci = GenomicRanges::reduce(union of H3K27ac + H3K27me3 consensus peaks) = 616,374 loci. (2) Per common region×cell-type (17 populations present in BOTH marks with >=50 cells), map each mark's mean pseudobulk CPM onto union loci (max over overlapping mark-peaks). (3) Per locus per cell type: active = ac>=1 CPM & me3<1; repressed = me3>=1 & ac<1; dual = both>=1; unclear = both<1. (4) Global state = dominant per-celltype state, EXCEPT loci active in some cell types AND repressed in others => "celltype_switched". Threshold 1 CPM = median of nonzero within-celltype signal.
**Reason.** H3K27ac & H3K27me3 mark the same residue (H3K27) and are chemically mutually exclusive on one nucleosome (Decision #03); they are also separate antibody libraries here (a nucleus carries one mark), so comparison is only valid at matched region×cell-type pseudobulk. Peak-set-overlap classification inflates "dual" to 68% because reduce() merges broad me3 domains (up to 40kb) with every narrow ac peak inside — an artifact, not biology. Signal-based per-cell-type classification respects mutual exclusivity and yields interpretable states.
**Pros.** Biologically grounded; per-cell-type resolution; separates true co-occurrence (dual, same population = heterogeneity/transitional) from cell-type switching. Result: active 11%, dual 19%, repressed 18%, celltype_switched 15%, unclear 37%. State classes have coherent properties (active/dual promoter-enriched; unclear narrow & promoter-poor).
**Cons.** "Dual" cannot be resolved to single-nucleosome bivalency (no H3K4 marks); 1-CPM threshold is a defensible but not unique cutoff; loci width heterogeneity (me3 broad) means some dual calls reflect promoter-in-domain overlap. Confidence: Medium-High.
**Confidence.** Medium-high; pass rates reported per library in csv/qc_summary_table.csv; will revisit if a cell type is lost.

## Decision #06 — De-novo clustering & annotation strategy
**Choice.** Cluster RNA per REGION (HIP, CBL separately), both histone marks together (same nuclei; Paired-Tag). Pipeline: LogNormalize, 2000 HVGs, PCA(30 PCs), Harmony integration over `library` (removes batch + mark-technical variation), SNN graph, Louvain resolution 0.6, UMAP. Annotate via Seurat AddModuleScore on specificity-weighted canonical marker panels (region-appropriate: HIP adds DG/CA ExN subtypes; CBL adds Granule/Purkinje/Bergmann/UBC), z-score cluster means across clusters, argmax → subtype, collapse to major_type.
**Reason.** Contract mandates de-novo clustering + marker interpretation (not label transfer). Per-region because cell-type repertoires differ (cerebellar granule/Purkinje absent in HIP). Marks clustered together so a single cell-type definition serves both CUT&Tag pseudobulks. Harmony over integration alternatives: contract prefers covariate/integration; Harmony is the standard scalable batch-integration for this cell count on CPU. Louvain (not Leiden) because Leiden needs the python leidenalg backend; Louvain res 0.6 gives well-separated major types.
**Pros.** Recovers expected region-specific biology (CBL ~75% granule; HIP oligo-rich); marks co-mingle post-Harmony (validated). **Cons.** Louvain not Leiden; module-score annotation depends on panel choice (mitigated by marker dotplot validation + confirmatory FindAllMarkers on separable clusters).
**Confidence.** High for major types; medium for rare subtypes (Purkinje, UBC/Golgi, ependymal) given low counts.

## Decision #12 — Glial subtype states: marker + reference-based, mito-QC exclusion
**Choice.** HIP astrocytes & microglia subclustered (res 0.3); states assigned by combining de-novo FindAllMarkers with AddModuleScore against published consortium signatures (DAM Keren-Shaul 2017; ROSMAP microglial states Sun 2023; DAA Habib 2020; SEA-AD MHC-II). Argmax reference score + marker cross-check. Low-quality mitochondrial-high subclusters (>6% MT-encoded reads, ~700 median genes), oligo/neuron-doublet, and lymphocyte-contamination subclusters excluded consistently in BOTH cell types. Fine subclusters collapsed to 2-3 states for power. Final — micro: Homeostatic 2664, Inflammatory 226, PVM/lipid-DAM 153; astro: Homeostatic 4442, Reactive 881, MT_stress 997.
**Reason.** User annotation ⑥ requested combined marker+reference subtyping citing named consortium sources. Mito-high clusters are QC artifacts (low nFeature, high mito) not biological states — folding them into Homeostatic (initial error) or MT_stress was wrong; metallothionein content is negligible except astro cl4 (genuine MT_stress).
**Pros.** Reference grounding makes state labels defensible/comparable to literature; consistent QC exclusion removes a confound. **Cons.** State proportions sensitive to exclusion thresholds; fine states collapsed so intra-state heterogeneity lost.
**Confidence.** High for homeostatic/reactive/activated identity; medium for the minority states (n cells small).

## Decision #13 — Per-state CRE DE filter & pooled GSEA readout
**Choice.** For subtype pseudobulk DE, use CPM-based filter (CPM≥1 in ≥3 samples) not filterByExpr; edgeR QL diagnosis-only, same thresholds as main DE. Report single-CRE result honestly (near-null) and treat pooled per-state GSEA (fgsea on gene-aggregated signed −log10P) as the interpretable readout.
**Reason.** filterByExpr sets its CPM threshold from the median library size; the small minority states (Reactive 1.5e5, Activated 7.6e4, MT_stress 1.4e4 vs Homeostatic 1.2-1.3e6) then retain ~0-1 peaks. CPM filter recovers 10k-400k testable peaks per stratum. State-splitting raises BCV (0.17-0.18 homeostatic → 0.63 MT_stress), so single-CRE power is low by construction; GSEA pools it.
**Pros.** Powered, interpretable state-resolved signal (astro state-switch; micro de-repression). **Cons.** GSEA inherits Enrichr immune-set bias; gene-aggregation (max|stat| per gene) is a lossy summary.
**Confidence.** High that signal is state-partitioned; medium on specific pathway identities given gene-set bias.

## Decision #14 — Report revision: drop two-axis model figure, reframe as region-divergent glial axis
**Choice.** Drop biological_model.png (old Fig 10) from Results per user annotation ③; reframe the "two AD-relevant axes" narrative (§4.5) as a single region-divergent glial immune axis — same inflammatory programs marked in opposite directions in HIP (reactive) astrocytes (gain) vs CBL microglia (loss). Replace weak GSEA_pathways Fig 8 with the per-cell-type enrichment scatter. Keep original prioritized_targets_heatmap (Fig 16) and ADD two per-cell-type target figures.
**Reason.** User judged the model figure and two-axis framing unsupported; the enrichment scatter (annotation ①) directly shows per-cell-type direction. Region-divergence is the more defensible, data-grounded story.
**Pros.** Narrative now matches the figures' evidence; state-resolved subtype analysis is the biological core. **Cons.** Loses the single summarizing schematic (mitigated: scatter + burden + state figures carry it).
**Confidence.** High.

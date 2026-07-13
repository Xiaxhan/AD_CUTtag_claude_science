# NOTEBOOK — AD CUT&Tag CRE

## 2026-07-08 — Session 1: setup & data validation

**Did.** Read Project.md, style.md, Donor_metadata.csv, mislabeled_samples.txt. Confirmed no prior artifacts/records exist (fresh start). Profiled donor demultiplexing across all 12 `barcode_assignments.csv.gz` files. Built `sc-epi` R env (Seurat 5.5.1, Signac 1.17.1, edgeR, DESeq2, limma, fgsea, clusterProfiler, EnsDb.Hsapiens.v86, BSgenome.Hsapiens.UCSC.hg38). Probed one raw Seurat object.

**Observed.**
- Design = clean 10 donors x 2 regions x 2 marks; donors pooled 3-4 per library; 12 libraries. 5 AD / 5 CTL.
- Raw Seurat object (CBL/H3K27me3/0506): assays `regemt` (1.64M feat), `peaks` (687k MACS2 feat), `RNA` (36,601 genes); ~524k RAW droplets/library. Meta has donor (`sample`), `info` (valid/background), QC (TSS proxy via FIE/FIP, nCount_*), but NO cell-type annotation, NO clustering, NO region/batch field.
- Region & mark are defined ONLY by folder path.
- Cell-type annotations existed only in one ancillary ArchR CSV (HIP/H3K27ac/LIB1) — not usable as delivered labels; will do de-novo clustering as mandated.

**Thought.** The mislabel file is the gating issue (data contract says stop & report). Decoded via donor library membership; 3 CBL entries map cleanly, 2 Br5090 entries are ambiguous. Brought both to user.

**Decisions (user-confirmed).** Reassign Br5106 CBL→HIP (pool with its HIP LIB2); reassign Br5090 HIP LIB3→CBL (pool with its CBL 0506); drop Br5196 CBL-me3 (cortex). See DECISIONS #01/#02. Corrected design: HIP arms n=9 (4AD/5CTL), CBL-ac n=9 (4AD/5CTL), CBL-me3 n=8 (3AD/5CTL).

**Next.** QC + valid-singlet filtering per library, then de-novo RNA clustering & annotation.

## 2026-07-08 — Session 1 cont: Step 1 QC & filtering DONE

**Did.** Discovered feature-matrix structure (Decision #04): RNA (36,601) & regemt (1.64M) identical across libs; MACS2 `peaks` library-specific (~700-860k, ~106 shared) → consensus+requant needed (user-confirmed). Ran QC across all 12 libs: dropped regemt on load (mem), kept info==valid demultiplexed singlets, applied thresholds (Decision #05), applied mislabel corrections, cached RNA counts + slim meta + MACS2 peak GRanges per lib.

**Observed.**
- No MT genes in this gene reference → percent.mt filter not applicable (recorded; not a bug).
- 147,228 QC-pass cells total. Per-arm (post-correction): HIP-ac 52,730 (9 donors, 4AD/5CTL); HIP-me3 35,010 (9, 4/5); CBL-ac 31,284 (9, 4/5); CBL-me3 28,204 (8, 3/5).
- Corrections verified: Br5090 pools 2 libs into CBL; Br5106 pools 2 libs into HIP; Br5196 CBL-me3 (cortex) dropped (CBL_H3K27me3_0513_03: 9,408→8,941 in that lib).
- Retained-cell QC distributions sit well above floors (RNA mode ~1.5-2k genes; CUT&Tag ~5-10k frags) — thresholds cleanly separate cells from empty droplets.

**Artifacts.** qc_metrics.png, qc_summary_table.csv, qc_donor_composition.csv. Cache: output/cache/<lib>_qc.rds (12 files).

**Next.** Step 2: de-novo RNA clustering + marker annotation per region.

## 2026-07-08 — Session 1 cont: Step 2 de-novo clustering & annotation DONE

**Did.** Per region (HIP, CBL), merged cached QC-pass RNA counts across contributing libraries (post-correction eff_region), LogNormalize → 2000 HVG → ScaleData → PCA(30) → Harmony(group=library) → SNN → Louvain(res=0.6) → UMAP. Annotated clusters by AddModuleScore over specificity-weighted canonical panels (region-appropriate), z-scored across clusters, argmax → subtype, collapsed to major_type. Both marks clustered TOGETHER per region (same nuclei measured), each cell keeps its mark label for downstream CUT&Tag pseudobulk.

**Observed.**
- Gene matrix uses ENSG rownames (not symbols) — built ENSG↔symbol map from features.tsv.gz (output/cluster/feature_map.rds).
- HIP: 87,740 nuclei → 26 clusters → 8 major types. Oligo 30,738; ExN 21,006 (DG+CA subtypes); Micro 6,557; OPC 6,723; InN 5,713; Astro 9,479; Vascular 3,714; Ependymal 3,810.
- CBL: 59,488 nuclei → 21 clusters → 9 major types. GranuleN 44,521 (~75%, textbook cerebellum); InN 4,479 (incl UBC/Golgi); ExN 3,527; Astro 2,978 (incl Bergmann); Oligo 1,527; Micro 733; OPC 742; Purkinje 488; Vascular 493.
- Integration check: the two histone marks co-mingle within every cell type on UMAP → Harmony removed antibody-technical variation; cell identity (not mark) drives clustering.
- Marker dotplots show clean specificity diagonal in both regions → annotations validated.
- FindAllMarkers (Wilcoxon, one-vs-rest, logfc>=0.5) yielded markers for only 4 HIP clusters (neuronal subclusters too similar); module-score approach used for annotation instead. Confirmatory: C0=MOG/ST18 (Oligo), C1=ALDH1A1/ADGRV1 (Astro), C4=PDGFRA/LHFPL3 (OPC).

**Flag.** CBL Purkinje H3K27me3 = 135 cells across all donors → thin for per-donor pseudobulk; apply min-cell floor in DE.

**Artifacts.** celltypes_HIP.png, celltypes_CBL.png, cell_annotations.csv.gz (147,228 cells; the downstream pseudobulk key), celltype_counts_by_arm.csv, {HIP,CBL}_cluster_annotation.csv. Objects: output/cluster/{HIP,CBL}_rna_annotated.rds.

**Next.** Step 3: build consensus MACS2 peak set per mark + re-quantify fragments → annotated CRE atlas (Deliverable 1).

## 2026-07-08 — Session 1 cont: Step 3 consensus peaks + re-quant + CRE atlas (Deliverable 1) DONE

**Did.**
1. tabix-indexed the 5 unindexed H3K27me3 fragment files (htslib 1.23.1) in place next to source (non-destructive; completes delivered data).
2. Built consensus MACS2 peak set per mark by merging (GenomicRanges::reduce) per-library MACS2 peak coordinates (cached from Step 1), restricted to standard chromosomes: H3K27ac 859,854 peaks (median 677 bp); H3K27me3 567,404 peaks (median 2,022 bp — broader, expected for repressive mark).
3. Re-quantified all 12 libraries over their mark's consensus via Signac FeatureMatrix (QC-pass cells only), aggregated to donor×celltype×region pseudobulk sums. ~2.5 hr on CPU. Per-lib pseudobulk in output/pseudobulk/.
4. Merged per-lib pseudobulk per mark (vectorized sparse group-sum; pooled Br5090 CBL & Br5106 HIP across their 2 source libs). H3K27ac 168 groups, H3K27me3 161 groups.
5. Annotated consensus peaks to genomic context (promoter[-2kb/+200]/genic/distal, priority order) + nearest protein-coding gene + dist-to-TSS via EnsDb.Hsapiens.v86 (hg38, strand-aware TSS).
6. Built annotated CRE atlas per mark: coords + context + nearest gene + total signal + per-region×celltype mean CPM.

**Observed.**
- Context split: H3K27ac 50.7% distal / 46.7% genic / 2.6% promoter; H3K27me3 53.2% / 43.2% / 3.6%. Typical enhancer-rich distribution.
- Known AD genes (APOE,BIN1,TREM2,MAPT,CLU,PICALM) all have associated consensus CREs in both marks.
- Atlas correlation heatmap: primary split by REGION, then by cell class within region → consensus CRE quantification captures real cell-type-specific active-enhancer biology.
- Cross-contamination artifacts flagged as tiny pseudobulk groups: HIP GranuleN (1 group, granule-like; clusters with CBL granule in heatmap = confirmed contamination), HIP Purkinje (6 cells), CBL Ependymal (5 cells). WILL EXCLUDE strata with too few donor-groups/cells from DE.

**Env note.** bioconductor-ensdb.hsapiens.v86 conda post-link download was blocked in sandbox (0-byte tarball); downloaded EnsDb.Hsapiens.v86_2.99.0.tar.gz from bioconductor.org/packages/release and R-CMD-INSTALLed into output/Rlib (user lib; add via .libPaths in annotation scripts).

**Artifacts.** CRE_atlas_{H3K27ac,H3K27me3}.csv.gz (Deliverable 1 tables; >1MB, reference by name), CRE_atlas_heatmap.png, CRE_atlas_context.png. Checkpoints: pb_merged_{H3K27ac,H3K27me3}.rds, peakann_{mark}.csv.gz.

**Next.** Step 4: CRE regulatory-state classification (active/repressed/dual-transitional/unclear; Decision #03).

## 2026-07-08 — Session 1 cont: Step 4 CRE regulatory-state classification DONE

**Did.** Built 616,374 unified CRE loci (reduce of ac+me3 consensus). Classified each locus per cell type from matched pseudobulk CPM (17 common region×celltype populations), then summarized to global state (method = Decision #07). Threshold 1 CPM.

**Observed.**
- Global states: unclear 37%, dual/transitional 19%, repressed 18%, celltype_switched 15%, active 11%.
- Peak-set-OVERLAP classification (rejected) gave 68% dual — artifact of reduce() merging broad me3 domains with narrow ac peaks. Signal-based per-celltype method respects K27 mutual exclusivity.
- State properties coherent: active/dual promoter-enriched (6.4%/9.1%), unclear narrow (639bp) & promoter-poor (1.1%), repressed/dual broad (3.9/5.7kb).
- AD-gene promoter states: APOE = celltype_switched (active in some types, repressed in others — matches known cell-type-specific APOE regulation); PICALM, PSEN1 = active; BIN1, CLU, MAPT, TREM2, APP = dual/transitional.

**Artifacts.** CRE_regulatory_states.csv.gz (per-locus global + per-celltype state, nearest gene; >1MB ref by name), CRE_state_classification.png.

**Next.** Step 5: pseudobulk differential CRE AD vs CTL (edgeR, per cell type×region×mark; exclude thin strata; CBL-me3 n=8 underpowered).

## 2026-07-08 — Session 1 cont: pseudobulk region-label BUG fix + Step 5 DE + pooled-power analyses

**Bug fixed (Decision #06b).** Pseudobulk merge had tagged region by source-library nominal folder, mis-placing reassigned donors (Br5090 HIP-LIB3 cells stayed HIP; Br5106 CBL-lib cells stayed CBL). Rebuilt merge deriving region per (donor×source-library) from cell_annotations.csv.gz. Corrected: 152 ac / 144 me3 groups (was 168/161). Br5090 CBL-only, Br5106 HIP-only. Spurious cross-region singleton strata removed. **Regenerated atlas (Deliverable 1) & CRE-state (Step 4) from corrected merge** — global state fractions essentially unchanged (unclear 36.8, repressed 19.5, dual 19.0, celltype_switched 13.6, active 11.0).

**Step 5 DE (Decision #08).** edgeR QL, ~diagnosis (NO batch covariate — near-confounded, up to 4 levels for 8-9 donors, some 1-donor levels). Donor pseudobulk = replication unit. 29 analyzable strata (≥3 AD & ≥3 CTL, ≥25 cells). HIP 4AD/5CTL; CBL 3-4AD/4-5CTL (CBL-me3 weakest 3AD/5CTL).
- **Single-CRE result is essentially NULL:** 1 CRE relaxed-sig genome-wide (distal H3K27ac near IMPAD1, HIP.InN, log2FC −1.38, FDR 0.02). p-value distributions at/below null; low dispersion (BCV 0.12-0.18 = good data). Cause: n=8-9 vs 70k-440k tests → near-zero single-CRE power after FDR. **Common, defensible outcome for this cohort size.**
- **CAUGHT ARTIFACT:** a batch-adjusted trial run gave "significant" GSTM1/USP17L4/TMEM106B CREs (FDR~2e-5) that ALL vanished (minFDR 0.66-0.74) under diagnosis-only → batch-confounding-as-signal. Excluded; diagnosis-only is the sound model.

**Pooled-power complements (Decision #09, user-approved).**
- Gene-level CRE aggregation DE: still near-null in well-powered strata (only CBL.Purkinje, a tiny 3v3/328-cell stratum, showed hits — likely unstable). Confirms null is genuine, not just a multiple-testing artifact.
- **Rank-based GSEA (fgsea, gene-level rankings vs Hallmark/KEGG/GO-BP/WikiPathways): 431 significant pathway-stratum enrichments (padj<0.1).** Strongest = **cerebellar microglia H3K27ac: 45 significant immune/inflammatory pathways, ALL with reduced enhancer marking in AD** (Interferon-γ response, TYROBP/DAP12 causal network, TNF-α/NF-κB, LPS response, cytokine production, complement, IL-2/STAT5). CBL.Micro also GAINS H3K27ac at synaptic/neuronal genes (likely ambient/identity effect — interpret cautiously).
- **Directionality — verified, corrected:** sign convention NES>0 = gained mark (GSTM1 sanity check logFC +1.41). The immune-enhancer-loss signal is SPECIFIC to CBL microglia. It does NOT replicate in HIP microglia: HIP.Micro H3K27ac has 0 significant immune pathways, and HIP.Micro H3K27me3 has only 2 significant immune pathways both with NEGATIVE NES (reduced repressive mark, i.e. de-repression tendency — NOT convergent silencing). So the finding is cell-type- AND region-specific, not a microglia-wide signature.

**Artifacts.** DE_volcanoes.png (v2), DE_summary.csv (v2), DE_significant_CREs.csv (v2, 1 row), DE_results_suggestive.csv.gz (v2, ref by name), DE_genelevel_results.csv.gz, GSEA_results.csv.gz, GSEA_pathways.png. Atlas & CRE-state CSVs re-versioned from corrected merge.

**Next.** Step 6 CRE-gene associations (peak-RNA correlation); then pathway/cell-type interpretation (Step 7 partly done via GSEA), biological model, glial subtypes, figures+report.

## 2026-07-08 — Session 1 cont: Step 6 CRE-gene associations DONE (Deliverable 3)

**Did.** Built RNA pseudobulk (36,601 genes × 152 donor×celltype×region groups, corrected regions). For each CRE with a nearest gene within 100kb, computed Spearman correlation of CRE signal (log CPM) vs that gene's expression across the shared pseudobulk groups (152 for ac, 144 for me3). Sig FDR<0.05.

**Observed (strong internal validation).**
- H3K27ac: 455,939 pairs tested; 163,756 significant POSITIVE (activating), only 3,054 negative. rho median +0.14, right-shifted. Active mark tracks target expression.
- H3K27me3: 281,494 tested; 40,449 significant NEGATIVE (repressive) vs 11,934 positive. rho median −0.04, left-shifted. Repressive mark anti-correlates with expression.
- AD-gene links coherent: APOE H3K27ac promoter/genic activating (rho +0.26..+0.45, 5 sig CREs) AND a H3K27me3 promoter CRE repressive (rho −0.28) — matches APOE cell-type-switched state (Step 4). APP H3K27me3 CREs repressive (rho −0.22..−0.36, 7 sig CREs), H3K27ac activating. ABCA7 promoter H3K27ac activating.

**Artifacts.** CRE_gene_associations.csv.gz (737,433 links: peak, gene, mark, context, dist_to_tss, rho, fdr, link-class; >1MB ref by name), CRE_gene_associations.png (rho density by mark; directional link counts; APOE enhancer-expression scatter).

**Next.** Step 7 cell-type pathway interpretation (extend GSEA), then biological model + prioritized validation targets, glial subtypes, figures+report.

## 2026-07-08 — Session 1 cont: Step 7 cell-type pathway interpretation DONE (Deliverable 4, Q4)

**Did.** Consolidated the pooled GSEA into a cell-type-resolved remodeling map: per-stratum sig-pathway burden + theme tagging (immune, synaptic, metabolic, development, proteostasis) with net direction (mean NES → mark gain/loss).

**Observed (answers Q4: which cell types remodel most, in what pathways).**
- Remodeling burden ranks: CBL.Micro H3K27ac 148 » HIP.ExN H3K27ac 79 > HIP.Oligo H3K27ac 69 > CBL.GranuleN H3K27me3 32 > CBL.GranuleN H3K27ac 29 > HIP.Astro H3K27ac 28 > HIP.InN 17 > HIP.Micro H3K27me3 12. Tail small.
- Immune theme (69 sig): LOST H3K27ac in CBL.Micro (37 pw, meanNES −1.84), HIP.ExN (−1.83), HIP.Oligo (−1.83); HIP.Micro H3K27me3 = repress_LOSS (−1.67, NOT silencing). CBL.GranuleN uniquely GAINS both ac and me3 at immune sets (mixed/atypical — interpret cautiously).
- Synaptic theme (34): reciprocal — CBL.Micro & HIP.Oligo GAIN H3K27ac at synaptic genes; HIP.Astro & CBL.GranuleN lose. (CBL.Micro synaptic gain likely ambient-neuronal signal, flagged.)
- Themes overall: Other 237, Immune 69, Development 60, Synaptic 34, Metabolic 23, Proteostasis 8.

**Interpretation.** The dominant, most coherent AD signal = loss of active-enhancer (H3K27ac) marking at microglial immune programs, strongest in cerebellum. Neurons (HIP.ExN) and oligodendrocytes also carry substantial H3K27ac remodeling. Repressive-mark (H3K27me3) remodeling is concentrated in CBL granule neurons (gain).

**Artifacts.** celltype_interpretation.png (burden barplot + theme×celltype direction heatmap), celltype_remodeling_burden.csv, celltype_pathway_themes.csv.

**Next.** Biological model + prioritized validation targets (Deliverable 5, Q5), glial subtypes, figures+report.

## 2026-07-08 — Session 1 cont: Step 8 biological model + prioritized targets DONE (Deliverable 5, Q5)

**Model (answers Q5).** Integrating all steps: AD does not produce genome-wide single-CRE differential signal at n=8-9 donors, but pooled analysis reveals coordinated, cell-type-specific CRE remodeling. The dominant, most coherent axis = **loss of H3K27ac (active-enhancer marking) at microglial immune/inflammatory programs, strongest in cerebellar microglia** (IFN-γ response, TYROBP/DAP12 network, TNF-α/NF-κB, complement, cytokine production). This is consistent with a shift of microglia away from a canonical inflammatory-enhancer state in AD cerebellum. Neurons (HIP.ExN) and oligodendrocytes carry additional H3K27ac remodeling; CBL granule neurons show the main H3K27me3 (repressive) remodeling. Caveats: cerebellum-specific (not replicated in HIP microglia); small microglial strata; Enrichr immune-set bias; single cohort.

**Prioritized validation targets (Deliverable 5).** Ranked CBL-microglia immune leading-edge genes (203) by AD effect (gene-level H3K27ac loss) × best activating CRE-gene link × AD-GWAS relevance. Top 15 each carry a SPECIFIC CRE coordinate to test. Top: CD86 (rho 0.57, chr3-122073933-122075566), BATF3, CYTL1, IL6 (AD-GWAS), PTGS2/COX2, NRP1, LTBR, PLA2G4A, MAP3K8, TNFAIP3 (AD-GWAS). All predicted to LOSE enhancer activity in AD cerebellar microglia.

**Artifacts.** biological_model.png (state landscape + remodeling burden + target dotplot), prioritized_validation_targets.csv (top 15 with CRE coords), validation_targets_microglia_immune.csv (full 203-gene ranking).

**Next.** Glial subtype clustering + preferential-remodeling sub-analysis, then figures + REPORT.md.

## 2026-07-08 — Session 1 cont: Step 9 glial subtype sub-analysis DONE

**Cell-number gate.** CBL Micro is genuinely too thin (733 cells, only 2/9 donors >=100 cells). CBL Astro (8/9 donors adequate) and CBL Oligo (7/9 donors adequate, median 119/donor) WOULD support subclustering by the same donor-adequacy criterion — the initial "CBL too thin" framing was over-broad and applied inconsistently (it excluded CBL Oligo on total-count 1,527 while keeping comparable strata). For this pass the subtype analysis was run on HIPPOCAMPUS only (Astro 9,479; Micro 6,557; Oligo 30,738; all 9 donors adequate), because HIP is the AD-vulnerable region of primary interest and the microglial state axis (the main glial finding) is only powered in HIP. CBL Astro/Oligo subclustering is a supportable follow-up, not excluded on data grounds.

**Did.** Subclustered each HIP glial major type (RNA, 2000 HVG, PCA 20, Louvain res 0.3), scored canonical state panels (Micro: homeostatic P2RY12/CX3CR1/TMEM119 vs activated APOE/CD68/SPP1/TREM2 vs inflammatory IL1B/CCL2/CD86; Astro: homeostatic SLC1A2/AQP4 vs reactive GFAP/VIM/CHI3L1; Oligo: mature PLP1/MOG vs stress HSPA1A/FOS/QDPR). Tested AD vs CTL per-donor reactive score (Wilcoxon).

**Observed.**
- Microglia: 9 subclusters resolve a homeostatic core (clusters 0/1, high P2RY12) + distinct activated (cl6) + inflammatory (cl5/7/8, high IL1B/CD86) states.
- Astro: ALL 9 subclusters argmax-labeled homeostatic (no discrete reactive subcluster resolved — reactive genes GFAP/CHI3L1 vary continuously, not as a separate cluster). Oligo: 7 mature + 1 stress subcluster (cl5). So a discrete reactive/stress state was resolved only in oligodendrocytes; astrocyte reactivity is a graded score, not a subcluster.
- **AD vs CTL reactive-state shift NOT significant** in any glial type: Micro p=0.73, Astro p=0.56 (AD median react 0.072 vs CTL −0.131, non-sig trend up), Oligo p=0.29. n=9 donors.

**Interpretation (KEY).** At this cohort size, the AD-associated signal is in the CRE/enhancer LANDSCAPE WITHIN glial cell types, NOT in gross subtype-proportion shifts. Microglial state structure is intact; astrocytes show a non-significant trend toward more reactive score in AD, congruent with (but not proving) the HIP.Astro NF-κB H3K27ac gain. Preferential remodeling of a specific subtype could not be established with confidence — flagged as underpowered, hypothesis-level.

**Artifacts.** glial_subtypes.png (microglia state UMAP + per-donor reactive score AD vs CTL for 3 glial types), glial_subtype_states.csv, glial_state_AD_vs_CTL.csv.

**Next.** Figures compilation + REPORT.md + plain-English methods.

## 2026-07-08 — Session 1 cont: Step 8 REVISED per user feedback (region balance, caution, per-context targets)

**User feedback (4 points), all addressed.**
1. Frame reduced H3K27ac at inflammatory programs more cautiously — as region-specific dampening / compensation / impaired immune responsiveness, NOT "silencing." → Adopted cautious language.
2. Hippocampal remodeling under-discussed; evaluate HIP as a major AD axis (HIP more AD-vulnerable than CBL). → Re-examined: HIP.ExN (79 sig) and HIP.Oligo (69) are the 2nd/3rd highest remodeling burden and ALSO show immune-enhancer LOSS (T-cell/LPS/TYROBP, negative NES). **HIP.Astro GAINS H3K27ac at TNF-α/NF-κB (+1.69)** with leading edge = CHI3L1/YKL-40, CCL2, TLR4/6/7, TNFRSF1A, RIPK1/2, NFKB2 = classic reactive astrogliosis in the vulnerable region. HIP is now a co-equal axis in the model.
3. Per-region×celltype prioritized targets. → Built prioritized_targets_by_context.csv (top 10 per each of 9 remodeling contexts).
4. Multi-layer heatmap of AD/CRE-linked genes × context. → prioritized_targets_heatmap.png (tile=gene log2FC, dot=CRE-gene link size/direction, bold=AD-GWAS).

**Also caught: sex-chromosome & CNV artifacts** in leading edges (USP17L family, Y-genes UTY/NLGN4Y/RPS4Y, GSTM1/GSTT1, HLA, immunoglobulin/OR/LILR/KIR loci, UGT2B, C4). Filtered out before ranking — these reflect donor sex imbalance / germline CNV, not AD regulation. 50+ genes removed.

**Revised model (two-axis, cautious).** AD is associated with cell-type- and region-specific CRE remodeling, not genome-wide single-CRE change. Two AD-relevant axes: (i) HIPPOCAMPUS (AD-vulnerable) — reactive-astrocyte H3K27ac GAIN at NF-κB/TLR inflammatory genes (CHI3L1/CCL2/TLR4), plus H3K27ac loss at immune programs in ExN & oligodendrocytes; (ii) CEREBELLUM — microglial H3K27ac reduction at immune/inflammatory enhancers (IFN-γ, TYROBP, TNF-α), interpreted cautiously as regional immune dampening / altered responsiveness rather than silencing. 21 AD-GWAS genes appear in leading edges across contexts (TREM2, TYROBP, SPI1, INPP5D, PLCG2, CD33, APOE, CR1, PTK2B, ABCA7, CASS4, IL6, TNFAIP3, MEF2C, SORL1, APP, ADAM10, EPHA1...).

**Artifacts (revised).** biological_model.png v2 (both regions, immune direction by celltype), prioritized_targets_heatmap.png, prioritized_targets_by_context.csv, prioritized_validation_targets.csv (CBL-micro top15 retained).

**Next.** Glial subtype clustering + preferential-remodeling sub-analysis, then figures + REPORT.md.

## 2026-07-09 — Report polish (Deliverables 6,7 revision)

**Did.** Rewrote REPORT.md into a full scientific-report structure: Abstract → Introduction (AD non-coding risk, myeloid enhancer enrichment) → Cohort/provenance → detailed Methods (§3, tools+versions+parameters) → Results (§4.0–4.7, all 13 figures embedded by artifact reference) → Prioritized targets → **Discussion** (literature-grounded) → Limitations → Future work → References (19, PubMed-verified DOIs) → Deliverables. Reorganized narrative around the two-axis biological story. Added a technical-parameter appendix to METHODS_plain_english.md.

**Literature (verified via PubMed metadata, real DOIs).** Method refs: Zhu 2021 / Xie 2023 (Paired-Tag), Stuart 2021 (Signac), Korsunsky 2019 (Harmony), Zhang 2008 (MACS), Robinson 2010 (edgeR), Subramanian 2005 (GSEA), Thurman 2021 (multi-subject sc-DE). Biology refs: Mathys 2019/2023 (sc AD atlas), Nott 2019 (PMID 31727856) + Novikova 2021 (AD myeloid enhancer heritability), Deczkowska 2018 + Podleśny-Drabiniok 2020 (DAM), Haure-Mirande 2018 + Audrain 2020 (TYROBP/DAP12), Sun 2022 (NF-κB in AD), Lananna 2020 (CHI3L1/YKL-40), Kong 2025 (spatial glial chromatin AD).

**Auditor fixes applied this revision.** (1) DOI btab337 correctly attributed to Thurman et al. (not Crowell/muscat). (2) Nott 2019 verified (PMID 31727856, Science, doi:10.1126/science.aay0793) before citing.

**Artifacts.** REPORT.md v3 (5e6e38b0), METHODS_plain_english.md v2 (8338a4f3).

---

## 2026-07-11 — Session 2: glial-subtype revision (reference+marker astro/micro) & report rewrite

**Trigger.** User annotations on REPORT.md v3 driving a 7-step revision (plan ec5e848c): figure legibility; reference+marker astro/micro subtyping; per-subtype composition + CRE change; per-cell-type enrichment scatter replacing weak Fig-8 claim; per-cell-type target figures (keep original Fig 13); Results+Discussion rewrite.

**Fig legibility (Step 1).** Rebuilt celltypes_HIP/CBL dotplots (cell types on x, italic markers on y, size=pct.exp, color=avg.exp.scaled — replaced broken DotPlot+coord_flip); bigger fonts. CRE_state_classification.png narrower bars + count+% labels. celltype_interpretation.png numeric bar labels + "NES = GSEA normalized enrichment score" spelled out.

**Subtype annotation (Step 2), marker + reference.** De-novo FindAllMarkers + AddModuleScore vs published signatures (DAM Keren-Shaul 2017; ROSMAP microglial states Sun 2023; DAA Habib 2020; SEA-AD MHC-II). Micro states: Homeostatic 2664, Inflammatory 226 (CCL3/4/TNFAIP3/CD83), PVM/lipid-DAM 153. Astro: Homeostatic 4442, Reactive 881 (CHI3L1/CXCL2/CCL2), MT_stress 997.
**Auditor fix (mito clusters).** micro cl1 (6.1% mito/717 feat) & astro cl3 (7.1%/700) are low-QUALITY not states → excluded consistently in BOTH cell types (metallothionein content negligible except astro cl4 = genuine MT_stress). Also excluded oligo/neuron-doublet & lymphocyte subclusters. Correction STRENGTHENED astro result.

**Composition AD vs CTL (Step 3), Wilcoxon n=9.** All reactive/activated up, homeostatic down (biologically expected). Significant: MT_stress astrocytes expand 0.004→0.104 (p=0.032). Others trend only: astro Homeostatic p=0.11, Reactive p=0.29; micro all p>0.7.

**Per-subtype CRE change (Step 4).** Re-quantified consensus H3K27ac(859854)+H3K27me3(567404) over donor×celltype×state pseudobulks (42 groups/mark; subtype_pb.rds). edgeR QL AD-vs-CTL, CPM≥1-in≥3 filter (filterByExpr fails on small-library states). Single-CRE NULL (1 me3 CRE Astro-Homeostatic FDR 0.003); BCV 0.17-0.18 (homeostatic) → 0.63 (MT_stress) confirming power drop with state-splitting. **Pooled per-state GSEA = key result:** astrocyte STATE-SWITCH — homeostatic LOSE inflammatory H3K27ac (TNF/NF-kB, IL-1, IL-2), reactive GAIN (TNF/NF-kB +1.63, IL-6 +1.71, IL-1 +1.75, Fc-γ, Oncostatin M). Micro: Activated GAIN immune H3K27ac (TNF/NF-kB +1.56, IL-2/STAT5, IL-6); Homeostatic LOSE repressive H3K27me3 at immune loci (IFN-γ −1.73, TNF −1.70, IL-6/JAK/STAT3, IL-7 = de-repression/priming).

**Env repairs during requant.** GenomeInfoDbData 1.2.13 (conda post-link left 0-byte placeholder) downloaded from mghp.osn.xsede.org archive (network access granted), R CMD INSTALL to abs Rlib path. Consensus GRanges seqinfo rebuilt under current GenomeInfoDb (Bioc 3.20; cached objects serialized under newer Bioc where Seqinfo class split into own package). Both fixes logged.

**Per-cell-type enrichment scatter (Step 5).** celltype_enrichment_scatter.png — columns=celltype×region×mark, rows=top+unique terms, size=−log10padj, color=NES. Replaces weak Fig-8 GSEA_pathways claim per annotation ①. CBL.Micro ac coherent immune LOSS; HIP.Astro ac gains; CBL.GranuleN me3 gains.

**Per-cell-type target figures (Step 6).** prioritized_targets_astrocyte.png + _microglia.png (2 panels each: prioritized CRE-linked genes + state-resolved immune direction). Original Fig 13 heatmap KEPT per user.

**Report rewrite (Step 7).** §4.4 grounded in per-cell-type enrichment scatter (Fig 8 now = scatter). §4.5 reframed two-axis → region-divergent glial immune axis; explained "AD-GWAS in leading edges". Dropped old Fig 10 (biological_model) per annotation ③; figures renumbered 1-17 continuous. §4.7 fully replaced with reference+marker subtype analysis (Figs 11-15). §5 targets renumbered + 2 new per-cell-type figures (Figs 16-17). Discussion rewritten around glial-state reactive switch + HIP-vs-CBL divergence + expected-vs-novel. Added refs Keren-Shaul 2017, Sun 2023, Habib 2020. REPORT.md v4.

**Auditor findings this session.** (1) composition "HIP n=9" allegedly not region-filtered → REBUTTED (CBL-named libs are Br5106 reassigned CBL→HIP; region field uniformly HIP). (2) BCV "0.13→0.63" mixed discarded filterByExpr run's low-end → corrected to 0.17-0.18→0.63 in report.

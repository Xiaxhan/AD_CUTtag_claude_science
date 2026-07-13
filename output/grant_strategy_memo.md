# Follow-up grant strategy memo: the glial-state chromatin switch in AD

*Strategic analysis built on `AD_CutTag_CRE/output/REPORT.md`. Framing, not proposal.*

## What the data establish (the floor)
1. **Distributed, not focal signal.** Single-CRE DE is null (1 element genome-wide); biology appears only at pooled pathway level. Every downstream claim is pathway-level and directional, never per-enhancer.
2. **Inflammatory program partitions by glial state.** HIP astrocytes: reactive state gains H3K27ac at NF-kB/IL-1/IL-6, homeostatic loses it. Microglia: activated state gains immune H3K27ac; homeostatic loses H3K27me3 (de-repression). Computed on donor x state pseudobulk -> partly abundance-controlled.
3. **Remodeled programs enriched for AD-GWAS genes** (TREM2, TYROBP, INPP5D, PLCG2, SPI1, CD33, APOE). Ties epigenomic signal to genetic architecture.

## Strongest mechanistic question
Is the homeostatic->reactive glial transition a locally-instructed, ORDERED chromatin switch -- de-repression (H3K27me3 loss) licensing enhancer activation (H3K27ac gain) -- set in motion by proximity to Ab/tau/neuronal injury, with Polycomb release as the rate-limiting, reversible licensing step?

Rationale: pairs the two independent observations (homeostatic-microglial H3K27me3 loss + reactive-state H3K27ac gain) into a sequence. Ordered, causal, druggable (PRC2/EZH2 -> p300/CBP-NF-kB). Literature-supported: PRC2/H3K27me3 maintains homeostatic microglial identity; Ab and cytokines drive H3K27ac/H3K4me3 trained-glia reprogramming.

## Can current data separate local transition from recruitment/migration? NO.
- No spatial coordinates (dissociated nuclei).
- Cross-sectional, single time point -> no direction.
- No lineage/trajectory anchor.
- Proportion tests underpowered (n=9); non-significant != no change. "Switch not abundance" is partly argument-from-absence.
- Defensible version: within-state AD-vs-CTL enhancer intensification is real and distinct from having more reactive cells.
- Indefensible version: "transition is in situ rather than recruited" -- data cannot address at all. This gap justifies the spatial follow-up.

## Spatial design that TESTS the model (not re-describes states)
Co-register in same tissue: glial enhancer/RNA state + Ab (6E10/methoxy-X04) + tau (AT8) + neuronal injury.
1. Dose-response: reactive-switch score vs distance-to-plaque / local tau density (graded gradient = locally instructed).
2. Find the intermediate: primed (H3K27me3-low, H3K27ac-not-yet-gained) microglia predicted in a peripheral ring.
3. Recruitment vs conversion: glial density + proliferation (recruitment) vs continuous spatial pseudotime homeostatic->primed->reactive (conversion). Need both.
4. Clearance vs vulnerability: co-localization of reactive-switch niche with vulnerable neurons. Field data (human plaque-niche ST) show low-Ab/high-glia spots have MORE neurodegeneration -> do not assume protective.

DIFFERENTIATOR: do it in the CHROMATIN modality. RNA-around-plaques is crowded (PIGs, DAM/DAA niches, SERPINA3+ astrocyte). Spatial-CUT&Tag / spatial-ATAC-RNA for H3K27ac/H3K27me3 exist and work in brain but are barely applied to human AD around pathology = white space your own finding is positioned to claim.

## Three-Aim structure (INDEPENDENT; shared input, not serial)

DESIGN PRINCIPLE: all three Aims draw on the SAME already-completed input -- the snPaired-Tag glial-state enhancer signatures (homeostatic/primed/reactive as H3K27me3/H3K27ac programs), the remodeled loci/programs (NF-kB/IL-1/IL-6, CHI3L1/CCL2/TLR-axis, TNF/IL-2/IL-6 module), and the nominated enzymatic axis (PRC2/EZH2 de-repression arm; p300/CBP-NF-kB acetylation arm). No Aim consumes another Aim's RESULT. Convergence is inferential, not operational: each answers a different necessary clause of the central question and is publishable alone.

- **Aim 1 - WHERE: is the switch spatially graded around pathology?** Input: completed signatures + pathology stains. Spatial epigenome-transcriptome co-profiling (H3K27ac + H3K27me3 + RNA) of human AD hippocampus, co-registered with Ab (methoxy-X04/6E10), tau (AT8), neuronal injury. Distance dose-response + primed-intermediate ring; density/proliferation vs spatial pseudotime for recruitment vs in-situ conversion. Prediction: monotonic rise toward pathology with de-repressed periphery; uniform switch falsifies local instruction. Stands alone: needs only existing molecular signatures.
- **Aim 2 - HOW: does local pathology instruct the switch, in what order, via which enzymes?** Input: completed signatures (readout) + PRC2/p300 axis (from data + published PRC2-maintains-microglial-identity biology). iPSC-microglia/astrocytes (+/- xenograft/organoid) challenged with defined Ab species, tau, injured-neuron-conditioned medium; time-resolved H3K27ac/H3K27me3 CUT&Tag for inducibility/ordering/stimulus-specificity. Perturb EZH2/PRC2 and p300/CBP-NF-kB for necessity + reversibility. Prediction: inducible, ordered (H3K27me3 loss precedes H3K27ac gain), blocked by EZH2 stabilization or NF-kB/p300 inhibition. Stands alone: reductionist system + enzymatic handles from existing data/literature, NOT Aim 1's gradient.
- **Aim 3 - SO WHAT: clearance or vulnerability? (two independent legs, either sufficient)** Input: PRC2/EZH2 handle (literature+data, not Aim 2) + completed signatures. (a) Cell-type-targeted perturbation of the de-repression step in vivo (amyloidosis/tauopathy or human xenograft) with spatial readout of amyloid/tau burden, plaque compaction, neuritic dystrophy, neuronal survival; (b) human correlation: reactive-switch-niche proximity vs neuronal loss (cf. low-Ab/high-glia -> more degeneration). Prediction as fork: forcing/blocking the step moves pathology and survival in opposite directions. Stands alone: perturbation target from published biology + own data; correlational leg needs only signatures + tissue.

INDEPENDENCE CHECK: none of Aim 1/2/3 requires a positive or specific result from another; each has its own model system, readout, and falsifiable prediction. A negative anywhere is informative, not fatal (e.g. uniform Aim-1 pattern + inducible Aim-2 switch = cell-intrinsic priming without spatial gradient).

## Evidence ledger
| Claim | Status | Honest read |
|---|---|---|
| Inflammatory program partitions by glial state | Supported | Pooled/pathway-level, state-stratified; 1 CRE sig genome-wide |
| Remodeled programs enriched for AD-GWAS genes | Supported | Leading-edge; strong genetic link |
| "Dominated by immune-enhancer changes" | Partly overinterpreted | Enrichr immune bias; largest burden = CBL-microglial loss in spared region. "Most coherent axis" fair; "dominated by" overweights biased annotation |
| "Switch rather than abundance shift" | Split | Within-state intensification supported; "rather than abundance" overreaches (proportions underpowered) |
| Homeostatic-microglial H3K27me3 loss = priming | Hypothesis | Pooled GSEA; interpretation, not demonstration |
| Ordered switch (de-repression licenses activation) | Hypothesis (best one) | Not testable here; engine of Aim 2 |
| Local Ab/tau drives transitions | Hypothesis | Zero support in current data |
| In situ vs recruitment | Cannot address | Design blind; Aim 1 |
| Clearance vs vulnerability | Cannot address | No neuronal-fate/causal axis; Aim 3 |
| Cerebellum as clean spared control | Caution | Largest burden in CBL microglia, possibly technical; do not lean on it without replication |

## Preliminary-data analyses runnable NOW on existing artifacts
1. Motif enrichment in gained-H3K27ac enhancers (NF-kB/RELA, AP-1, PU.1/SPI1) -> nominate driver TF, make Aim 2 axis concrete.
2. Overlap of de-repressed (H3K27me3-lost, homeostatic microglia) loci with activated (H3K27ac-gained, reactive/activated) loci -> direct in-data support for ordered-switch hypothesis; candidate Figure 1.
3. Per-donor reactive-enhancer score vs donor pathology staging (Braak/CERAD/Thal/Ab), if metadata exists -> within-cohort dose-response, closest to "pathology drives switch" pre-spatial.

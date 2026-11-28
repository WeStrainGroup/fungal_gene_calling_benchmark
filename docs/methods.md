# Benchmark methods

This protocol reuses frozen inputs and native predictions. It does not retrain models, change their native postprocessing, or rerun successful inference.

## Inputs and cohorts

The benchmark contains 17 unique fungal genomes selected from the ANNEVO, Tiberius, Carbon-A and Helixer evaluations. Shared assemblies are counted once. All 17 species, five conditions and four tools have completed prediction and independent validation. Carbon-A's five selected species are officially `seen`; no uniformly unseen-species claim is made. See [accessions and exposure](datasets.md).

Conditions are original assemblies and full-genome random partitions of 1–2 kb, 2–5 kb, 5–10 kb and 10–100 kb. Existing DNA, annotations, genetic codes, seeds and retained maps are frozen. For the five new assemblies, each condition uses an independent `random.Random(species_seed)`. Cut points do not follow genes. All tools receive identical DNA. Main scoring uses records >=1,000 bp; smaller original records are retained and counted in reconstruction checks.

The evaluator checks RefSeq MD5 files, frozen input SHA256 where available, source DNA equality, exact consecutive fragment coverage, each fragment's sequence checksum, and reconstruction of every original CDS from its visible pieces. Native prediction hashes must remain unchanged during evaluation and independent validation.

## Reference qualification and strata

A CDS means the spliced coding region of one transcript. Original reference qualification is unchanged: the record must be in the frozen nuclear input, initial phase must be zero, no original partial annotation, valid introns, ATG start, length divisible by three, a terminal stop under the sample's genetic code, and no internal stop. Ambiguous bases follow the historical rules; no new ambiguity filter is introduced. Explicit stop-codon features are unioned with CDS coordinates to harmonize GFF conventions without double counting them.

Completeness is derived from the original qualified annotation and fragment map, never from predicted codons or native model status:

- Complete: all original coding bases remain in the fragment.
- Truncated: some original coding bases fall outside it; distinguish missing 5' end, missing 3' end and missing both ends in transcriptional orientation.
- Cuts in introns count as truncation when coding exons lie outside the fragment, even if no visible CDS touches its physical edge.
- A fragment may contribute to both strata. No whole-fragment exclusion based on the presence of a truncated gene is applied.

Retain all eligible isoforms for gene matching, with at most one TP per original gene in each fragment. A gene with qualified isoforms in both complete and truncated classes is unscored. If all its qualified visible isoforms are <200 bp, it is unscored. Short isoform extensions outside an eligible gene span are masked; shared regions already represented by an eligible isoform are not converted to background. Truncated genes whose isoforms have different missing-end types remain in aggregate truncated scores, with subtype `mixed_ends`.

## Filters and strand-aware scoring domains

All tools retain their native filters. Apply the following external filters once, before every metric:

1. Input record >=1,000 bp.
2. Visible spliced CDS >=200 bp, counting retained terminal-stop bases and partial boundary bases.
3. Carbon-A confidence >=0.1 or explicit `NA`. A missing confidence attribute is an error, not `NA`.

No model is removed because its predicted start/stop, complete/incomplete status, or phase is wrong. An inconsistent predicted exon-phase chain cannot earn a gene TP, but remains in counts and interval metrics. Raw confidence, native status, inherited native partial flag, gene Parent and phase are read directly from the original GFF; prediction IDs remain unchanged.

Reference territory is the visible CDS span on each record and strand, including internal introns. Same-strand overlapping spans are resolved as connected components. A component containing complete and truncated references, or any unscored reference, is conservatively unscored in full; affected references are removed from scoring too. This conservative propagation prevents ambiguous FP attribution. A failed isoform can therefore mask an overlapping qualified isoform/gene; exclusion reasons and regions are provided. Opposite strands are independent.

| Prediction location after filtering | Treatment |
| --- | --- |
| Within complete territory, or elsewhere in scoreable background | Complete comparison; unmatched model is FP. |
| Wholly within truncated territory | Truncated-region evaluation; report scores for Carbon-A and Helixer. |
| Overlaps an unscored region | Unscored; report count, model CDS bp and proportion. |
| Overlaps truncated territory but extends outside it | Cross-scope; report separately without clipping it into a new model. |

Background FP is a subset of complete-scope predictions whose span overlaps no complete reference territory. Predictions crossing a complete-gene/background boundary remain in that common complete scope. Predictions crossing excluded/truncated boundaries do not. A reference not recovered by a scoreable prediction still contributes FN. No prediction is silently discarded; mutually exclusive bucket counts reconcile to the raw model count. Diagnostics `cds_bp` sums model coding lengths (not unique genomic coverage); overlapping models can contribute the same base more than once. `percent_of_raw` uses all parsed native CDS models as denominator. Background and unresolved subtype rows are additional subsets, not additive buckets.

## Metrics

Gene/CDS TP requires identical record, strand, full visible CDS interval chain and initial phase, with a consistent predicted exon-phase chain. Adjacent/overlapping CDS/stop features are normalized to disjoint intervals. On a truncated transcript, the expected phase is `(original_phase - removed_5prime_coding_bases) mod 3`; residual boundary bases are retained in the coordinate comparison. Identical sequence at a different locus is not a TP.

A maximum one-to-one matching permits any eligible isoform but counts a reference gene only once per fragment. Duplicate predictions cannot create extra TPs. A separate SciPy maximum matching checks each count. `TP + FP` equals the number of scoreable predicted models; `TP + FN` equals the number of scoreable reference genes.

Nucleotide and exon metrics keep the frozen official evaluator definitions: choose the longest eligible visible isoform per reference gene (transcript ID breaks length ties), merge intervals on each record/strand, calculate nucleotide overlap, and match exact merged exon intervals. Thus nucleotide/exon counts do not have the same denominators as gene counts, and exon unions are not counts of transcript-specific exon occurrences. No prediction is clipped into a new model. All three levels use the same confidence/length-filtered prediction collection within their respective scope.

Carbon-A and Helixer truncated precision/F1 are conditional within truncated-reference regions, not whole-genome precision/F1. For subtype results, matched predictions follow their globally matched reference gene's subtype. Unmatched predictions wholly attributable to one subtype territory contribute its FP; ambiguous subtype assignments remain in aggregate truncated counts and `subtype_unresolved`. The three named subtype counts therefore need not sum to the aggregate. Background and cross-scope predictions are always reported alongside these scores.

Each cell includes TP, FP, FN, predicted/reference denominator, precision, recall and F1. No eligible reference gives `NA` rates (even if FP are present). A reference with no prediction gives recall/F1 zero. Counts remain available when rates are NA. Under the revised aggregation policy (v0.2.1), each tool, input condition, stratum and metric is summarized across the unique species in the declared result set by arithmetic mean and sample SD (n−1), giving each species equal weight. Undefined values are excluded, never zero-filled; `n_species` states the number of defined values. With zero defined values, mean and SD are NA; with one, mean is retained and SD is NA. Original assemblies are a separate condition. All inference/validation must finish before final summary tables are produced. The exact species × input × tool × stratum × level grid is asserted.

## Original-assembly length groups

`compare_original_assemblies.R` joins `results/input_lengths.tsv` (`group=genome`, `scope=scored`) to complete-reference F1 scores in `results/metrics.tsv` (`group=genome`, `stratum=complete`). Species are assigned by the unrounded median length of their ≥1 kb nuclear records: <10,000 bp, 10,000–100,000 bp inclusive, or >100,000 bp. The groups contain 5, 4 and 8 species. Each species is paired between Tiberius and Carbon-A at each of the three scoring levels; means give every species equal weight.

Differences are Carbon-A − Tiberius on the F1 percentage scale. Two-sided paired Wilcoxon signed-rank tests use exact P values without continuity or multiple-comparison correction. The current nine comparisons have no zero differences or tied absolute differences; the script checks this and independently verifies every P value by enumerating all signed-rank assignments. With 5 or 4 pairs, the smallest possible two-sided exact P values are 0.0625 and 0.125, respectively, so P ≥0.05 does not establish equivalence. These comparisons are exploratory: species, gene structure, reference quality and training overlap differ among length groups, and a median by record count does not describe the length distribution experienced by most genes.

Outputs are `results/original_assembly_pairs.tsv` (51 species-level pairs), `results/original_assembly_tests.tsv` (9 comparisons) and `figures/original_assembly_comparison.{png,pdf}`. The figure uses the existing palette, `theme_minimal()`, solid circles (`shape=16`) and an 11.5 × 4 inch canvas.

## Validation and outputs

`tests/test_stratified.py` covers positive/negative strand clipping, intron cuts, all three missing-end types, phase, 199/200/201-bp CDS and 999/1,000-bp inputs, confidence boundaries and missing attributes, native incomplete predictions, duplicate and alternative-isoform matching, wrong locus/strand/frame, overlap conflicts, mixed fragments, ignored short references and empty-reference cases.

`validate_stratified.py` independently parses raw prediction and reference GFFs without importing the evaluator. It verifies every assigned gene TP's coordinates and phase, checks duplicate reference assignments, independently recounts complete/aggregate-truncated gene TP/FP/FN and rates, and rehashes all sources. Subtype and interval calculations are covered by boundary tests and the pinned interval evaluator; the independent validator does not independently reimplement those metrics. `validation_examples.tsv` retains real complete/truncated and cross-scope/unscored examples.

| File | Contents |
| --- | --- |
| `metrics.tsv` | Per-species, input, tool, stratum and level counts/rates. |
| `mean_sd.tsv` | Final mean and sample SD across the unique species; defined-value count per row. |
| `diagnostics.tsv` | Filtering, background, cross-scope, ignored and subtype-unresolved predictions. |
| `inputs.tsv` | Included/excluded record counts and bp. |
| `reference_regions.tsv` | Reference genes/regions, strata and exclusion reasons. |
| `reference/<species>/qualification.tsv` | Original reference transcript qualification. |
| `cases/<species>/<group>/<tool>/prediction_audit.tsv` | Every native model, confidence/status/phase, disposition and matched reference. |
| `validation.tsv`, `independent_validation.json` | Per-case checks and independent validation receipt. |
| `source_sha256.json` | Exact inputs and native predictions read. |
| `STATUS.json` | Evaluator version, completion and explicit missing cases. |

Current figures use R/ggplot2, `theme_minimal()`, the accepted muted palette, species-equal means and sample-SD error bars. The final 17-species/four-tool grid contains 340 inference cases and 3,060 metric rows. All 3,060 metric rows and 540 summary rows are complete; independent validation checked 1,603,921 gene TPs and 544 unchanged source files. No partial run is labeled complete. Native source files are reused; derived outputs have one maintained location.

## Implementation

| Script | Purpose | Status |
| --- | --- | --- |
| `install_helixer.sh` | Install fixed Helixer, HelixerPost, GeenuFF and fungi model in isolation. | Installed; 33 official tests and model MD5 passed. |
| `prepare_new_inputs.py` | Fetch exact RefSeq DNA/GFF, verify MD5, exclude organelles and validate all current partitions. | Five added species frozen; official MD5 and all partition reconstruction checks passed. |
| `prepare_benchmark.py` | Link reusable inputs and native predictions into one current view; write the 160-task list. | 17-species view ready; 180 existing cases reused without copying. |
| `run_benchmark.py` | Run pilot or pending cases on selected cooperative GPU locks; score and validate on completion. | All 340 predictions complete; CPU scoring and independent validation passed on 2026-11-24. |
| `run_hybrid_efficiency.py`, `hybrid_efficiency_ais.sh` | Measure CTC on a frozen original assembly. | Uses the [standalone CTC repository](https://github.com/WeStrainGroup/combined_tiberius_carbon_a) via `--ctc-script`; version and source checksum recorded; frozen deployed jobs unchanged. |
| `run_predictor.sh` | One native predictor invocation, four tool choices. | Existing parameters preserved; Helixer CLI checked against v0.3.8. |
| `prepare_efficiency.py` | Keep all original records >=1 kb; sample 50 frozen fragments per species/range in five equal-width length bins. | All 85 inputs prepared; 3,400 sampled fragments match frozen sequence hashes. |
| `check_efficiency_resources.py` | Verify A100 40GB, BF16, 32 GiB memory, and either four physical cores (Slurm) or four-CPU quota (AIS). | Both platform modes record allocation evidence; AIS does not claim exclusive cores. |
| `run_efficiency.py` | Check frozen inputs/software; run 20 serial cases per species and record runtime, CPU time and maximum individual-process RSS. Optional path configuration supports external installations. | All three rounds passed: 1,020 measurements; portable path configuration checked with a mocked 20-case run. |
| `efficiency_ais.sh` | Run one species on AIS with the frozen software prefixes; 20 serial cases per task, separate round outputs. | All 51 species-round tasks completed; output and hardware checks passed. |
| `summarize_efficiency.py` | Validate three full rounds; calculate within-input replicate medians, then equal-species means and sample SDs for throughput/time; export original-assembly relative elapsed times. | Final three-round summary verified; 340 medians and 20 throughput mean/SD summaries independently checked. |
| `collect_efficiency.py` | Collect one complete 17-species runtime round and verify native GFFs, hardware receipts and timing files. | All 1,020 GFFs and 51 hardware receipts verified; software/model receipts identical across tasks. |
| `plot_efficiency.R` | Plot throughput mean ± species SD by input condition, using the accepted palette and theme_minimal; assembly time is a table. | Final three-run-median throughput figure rendered and inspected. |
| `common.py` | FASTA helpers and frozen interval evaluator. | Current dependency. |
| `evaluate_stratified.py` | Reference strata, shared filtering, exact structural matching and metrics. | Four-tool support; inherited partial flag retained. |
| `validate_stratified.py` | Independent native-GFF coordinate/phase checks, count reconciliation and hashes. | Supports both truncated-prediction tools. |
| `summarize_all_species.py` | Deduplicate, verify selected tool/species grid, calculate species-equal mean/sample SD. | Four-tool grid regression passed. |
| `summarize_input_lengths.py` | FASTA length counts, quartiles and N50, per species and pooled; all/scored subsets. | All 85 input totals reconciled; 180 quantile summaries independently checked against NumPy and fragment maps. |
| `plot_all_species.R` | Complete comparison, shared truncated comparison and subtype plots. | Four tool positions in the main plots; ANNEVO/Tiberius truncated scores shown as NA (not evaluated). |
| `compare_original_assemblies.R` | Original-assembly length groups, paired F1 comparisons, exact raw P values and figure. | 51 pairs and 9 tests reproduced; exact P values independently enumerated; layout inspected. |
| `validate_summary.R` | Independently recompute means, sample SD and denominators in R. | All 540 expanded summary rows independently verified in R. |
| `tests/test_stratified.py`, `tests/test_species_summary.py` | Scientific boundary and summary regression tests. | 34 tests passed; 30 real cases retain all 210 metric rows exactly, with 130,930 gene TPs independently checked. |

Run status is recorded in the active compute root's `STATUS.json`; pilot status is `PILOT_STATUS.json`. The input, software and prediction manifests identify the exact sources. README and this registry must be updated whenever script status changes.

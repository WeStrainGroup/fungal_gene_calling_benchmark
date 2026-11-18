# Benchmark methods

This protocol reuses frozen inputs and native predictions. It does not retrain models, change their native postprocessing, or rerun successful inference.

## Inputs and cohorts

The source datasets contain 2 ANNEVO-selected, 6 Tiberius-selected and 5 Carbon-A-selected species. The current summary combines their 12 distinct genomes, counting the shared Saccharomyces cerevisiae inputs and predictions once. Source-specific tables remain provenance records, not the final reporting groups. Carbon-A's five species are all officially `seen`. See [accessions and training exposure](datasets.md); this is not a three-model unseen-species benchmark.

Original assemblies and full-genome random partitions of 1–2 kb, 5–10 kb and 50–100 kb use the existing frozen maps, seeds, annotations and genetic codes. Cut points are not moved to preserve genes. All tools receive identical DNA. Main scoring uses records >=1,000 bp. Smaller original records are retained and reported, not scored: Agaricus bisporus has 771 such records totaling 560,320 bp. Whole-genome reconstruction checks include these records.

The evaluator checks RefSeq MD5 files, frozen input SHA256 where available, source DNA equality, exact consecutive fragment coverage, each fragment's sequence checksum, and reconstruction of every original CDS from its visible pieces. Native prediction hashes must remain unchanged during evaluation and independent validation.

## Reference qualification and strata

A CDS means the spliced coding region of one transcript. Original reference qualification is unchanged: the record must be in the frozen nuclear input, initial phase must be zero, no original partial annotation, valid introns, ATG start, length divisible by three, a terminal stop under the sample's genetic code, and no internal stop. Ambiguous bases follow the historical rules; no new ambiguity filter is introduced. Explicit stop-codon features are unioned with CDS coordinates to harmonize GFF conventions without double counting them.

Completeness is derived from the original qualified annotation and fragment map, never from predicted codons or native model status:

- **Complete:** all original coding bases remain in the fragment.
- **Truncated:** some original coding bases fall outside it; distinguish missing 5' end, missing 3' end and missing both ends in transcriptional orientation.
- Cuts in introns count as truncation when coding exons lie outside the fragment, even if no visible CDS touches its physical edge.
- A fragment may contribute to both strata. No whole-fragment exclusion based on the presence of a truncated gene is applied.

Retain all eligible isoforms for gene matching, with at most one TP per original gene in each fragment. A gene with qualified isoforms in both complete and truncated classes is unscored. If all its qualified visible isoforms are <200 bp, it is unscored. Short isoform extensions outside an eligible gene span are masked; shared regions already represented by an eligible isoform are not converted to background. Truncated genes whose isoforms have different missing-end types remain in aggregate truncated scores, with subtype `mixed_ends`.

## Filters and strand-aware scoring domains

All tools retain their native filters. Apply the following external filters once, before every metric:

1. Input record >=1,000 bp.
2. Visible spliced CDS >=200 bp, counting retained terminal-stop bases and partial boundary bases.
3. Carbon-A confidence >=0.1 or explicit `NA`. A missing confidence attribute is an error, not `NA`.

No model is removed because its predicted start/stop, complete/incomplete status, or phase is wrong. An inconsistent predicted exon-phase chain cannot earn a gene TP, but remains in counts and interval metrics. Raw confidence, native status, gene Parent and phase are read directly from the original GFF; prediction IDs remain unchanged.

Reference territory is the visible CDS span on each record and strand, including internal introns. Same-strand overlapping spans are resolved as connected components. A component containing complete and truncated references, or any unscored reference, is conservatively unscored in full; affected references are removed from scoring too. This conservative propagation prevents ambiguous FP attribution. A failed isoform can therefore mask an overlapping qualified isoform/gene; exclusion reasons and regions are provided. Opposite strands are independent.

| Prediction location after filtering | Treatment |
| --- | --- |
| Within complete territory, or elsewhere in scoreable background | Complete comparison; unmatched model is FP. |
| Wholly within truncated territory | Truncated-region evaluation; report scores only for Carbon-A. |
| Overlaps an unscored region | Unscored; report count, model CDS bp and proportion. |
| Overlaps truncated territory but extends outside it | Cross-scope; report separately without clipping it into a new model. |

Background FP is a subset of complete-scope predictions whose span overlaps no complete reference territory. Predictions crossing a complete-gene/background boundary remain in that common complete scope. Predictions crossing excluded/truncated boundaries do not. A reference not recovered by a scoreable prediction still contributes FN. No prediction is silently discarded; mutually exclusive bucket counts reconcile to the raw model count. Diagnostics `cds_bp` sums model coding lengths (not unique genomic coverage); overlapping models can contribute the same base more than once. `percent_of_raw` uses all parsed native CDS models as denominator. Background and unresolved subtype rows are additional subsets, not additive buckets.

## Metrics

Gene/CDS TP requires **identical record, strand, full visible CDS interval chain and initial phase**, with a consistent predicted exon-phase chain. Adjacent/overlapping CDS/stop features are normalized to disjoint intervals. On a truncated transcript, the expected phase is `(original_phase - removed_5prime_coding_bases) mod 3`; residual boundary bases are retained in the coordinate comparison. Identical sequence at a different locus is not a TP.

A maximum one-to-one matching permits any eligible isoform but counts a reference gene only once per fragment. Duplicate predictions cannot create extra TPs. A separate SciPy maximum matching checks each count. `TP + FP` equals the number of scoreable predicted models; `TP + FN` equals the number of scoreable reference genes.

Nucleotide and exon metrics keep the frozen official evaluator definitions: choose the longest eligible visible isoform per reference gene (transcript ID breaks length ties), merge intervals on each record/strand, calculate nucleotide overlap, and match exact merged exon intervals. Thus nucleotide/exon counts do not have the same denominators as gene counts, and exon unions are not counts of transcript-specific exon occurrences. No prediction is clipped into a new model. All three levels use the same confidence/length-filtered prediction collection within their respective scope.

Carbon-A truncated precision/F1 are **conditional within truncated-reference regions**, not whole-genome precision/F1. For subtype results, matched predictions follow their globally matched reference gene's subtype. Unmatched predictions wholly attributable to one subtype territory contribute its FP; ambiguous subtype assignments remain in aggregate truncated counts and `subtype_unresolved`. The three named subtype counts therefore need not sum to the aggregate. Background and cross-scope predictions are always reported alongside these scores.

Each cell includes TP, FP, FN, predicted/reference denominator, precision, recall and F1. No eligible reference gives `NA` rates (even if FP are present). A reference with no prediction gives recall/F1 zero. Counts remain available when rates are NA. Under the revised aggregation policy (v0.2.1), each tool, input condition, stratum and metric is summarized across the 12 unique species by arithmetic mean and sample SD (n−1), giving each species equal weight. Undefined values are excluded, never zero-filled; `n_species` states the number of defined values. With zero defined values, mean and SD are NA; with one, mean is retained and SD is NA. Original assemblies are a separate condition. All inference/validation must finish before final summary tables are produced. An explicitly labeled interim report may instead retain only species-condition pairs with all three tools and all scoring rows validated. The same species are then used by all tools within a condition, while species counts and composition can differ between conditions. Interim mean/SD is calculated separately within each condition; it is not a paired estimate of fragmentation loss. Unpaired completed tool cases remain outside interim scores, without being treated as zero. The exact species × input × tool × stratum × level grid is asserted.

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

Current figures use R/ggplot2, `theme_minimal()`, the accepted muted tool palette, species-equal means and sample-SD error bars. `summarize_all_species.py` checks duplicate species for exact agreement, preserves zero values and asserts the final 1,008-row metric grid (144 unique inference cases). Current public tables live in `results/`, and final figures are generated in `figures/`. Detailed reference and prediction audits remain with the validated computation outputs. Obsolete derived results and temporary experiments are removed; primary inputs and native predictions are retained. Incomplete results are explicitly labeled and are not summarized as complete.

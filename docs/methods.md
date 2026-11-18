# Methods

## Matched inputs

The initial set uses S. cerevisiae GCF_000146045.2_R64 (16 nuclear records; 12,071,326 bp; code 1) and C. albicans GCF_000182965.3_ASM18296v3 (8 records; 14,282,666 bp; code 12). All three callers receive identical uppercase DNA. Source soft masking is removed; no evidence-based masking is added.

Every nuclear base is retained in nonoverlapping partitions targeting 1–2 kb, 5–10 kb and 50–100 kb. Original assemblies form the baseline. The frozen maps record source record, coordinates and sequence SHA256. Seeds are 20261117 and 20261118 for the initial two species; the six-species extension uses 20261119–20261124. When an original record is shorter than a target interval, it is retained unchanged and explicitly counted as below target; records are never joined across assembly boundaries. No 1 Mb subsampling or sequence loss is allowed.

## Reference qualification

Before fragmentation, retain all RefSeq CDS isoforms satisfying: a retained nuclear record and a single record/strand per transcript; initial phase 0 (or unspecified); no partial annotation; introns at least 2 bp; ATG start; length divisible by three; terminal stop under the sample genetic code; no internal stop. Exclusions are recorded by transcript. This is our reference policy, not the official Carbon-A evaluator's default filtering.

The initial set retains 5,955 and 5,974 genes, each with one qualified CDS isoform. An otherwise biological prediction corresponding to an excluded reference can count as a false positive under this policy.

Project **every** qualified isoform to each fragment, keeping every positive-length CDS intersection. No post-clipping completeness or minimum-length filter is applied. A reference gene unit is `(fragment record, original gene)`. Keep separate gene identities even when clipped sequences coincide. The phase accounts for removed coding bases, and all clipped pieces must reconstruct the original CDS sequence.

## Prediction normalization

Use each caller's native GFF. Extract the union of CDS and explicitly annotated `stop_codon` intervals from the exact shared input DNA. Use full IUPAC reverse complementation. Never extend coordinates using reference truth. Keep terminal stop codons in nucleotide CDS sequences; do not trim initial phase for exact CDS-sequence equality.

This scores the annotated canonical CDS rather than trusting native FASTA exports. Carbon-A's original complement helper handles only ATCGN; Tiberius's native export can mask other ambiguous bases and translate using the standard code. Those native export differences do not automatically imply different annotated coordinates.

All native predicted models are retained: confidence threshold **0**, external minimum CDS length **0**, no external codon filter. Native caller postprocessing still applies: Carbon-A minimum CDS 96 bp; Tiberius minimum CDS 200 bp/internal-stop filtering and input minimum 1,000 bp; ANNEVO native score and intron rules. Thus the comparison includes each tool's native pipeline policy.

## Three metrics

| Level | Unit and matching |
| --- | --- |
| Nucleotide | Choose the longest visible qualified CDS isoform per reference gene. Union intervals by record and strand; TP is the overlapping base count. |
| Exon | Using the same reference selection, merge overlapping or adjacent CDS intervals on each record/strand, including across genes. TP requires exact merged interval boundaries. This is not a per-transcript exon score. |
| Gene | A predicted entire spliced CDS nucleotide sequence must equal any qualified reference isoform on the same record. Use official maximum one-to-one matching; no extra locus or strand restriction. |

P = TP / (TP + FP); R = TP / (TP + FN); F1 = 2TP / (2TP + FP + FN). Percentages are reported on a 0–100 scale. The complete table includes every TP/FP/FN and denominator.

**Truncation example:** a 1,500 bp original CDS leaves 900 bp on an input fragment. Predicting exactly those 900 bp can be a gene-level TP. Predicting only 600 bp does not exactly match that reference. Restoring the missing 600 bp outside the input is neither possible nor required. A partial-CDS TP is not a full-length protein suitable for structure prediction.

The unchanged Carbon-A evaluator is pinned to revision `ee6359990109eef6adccfd19aba8ca64b48ae514`. Each wrapper result must agree with the official `evaluate_prediction` entry point. SciPy's maximum bipartite matcher independently validates gene TP counts. Fixtures cover alternative isoforms, duplicate sequences, record restriction, strand separation and adjacent-interval merging.

## Aggregation

Keep the ANNEVO-selected, Tiberius-selected and any Carbon-A-selected datasets separate. For each set × tool × input length × evaluation level, compute P, R and F1 per species first, then take their arithmetic means. Each species has equal weight. Mean F1 is the mean of species F1 values, not the harmonic mean of mean P and mean R. Do not pool base counts, pool the benchmark sets, compute SD/error bars, or silently average an incomplete species grid.

## Limits

This is an in-silico fragmentation experiment on reference assemblies. It does not test assembly indels, chimeras, contamination, strain mixtures, missing MAG sequence or taxonomically representative fungal sampling. There is one fixed partition per size condition. Training exposure is documented separately: the initial two yeasts trained Tiberius, and Carbon-A's official fungal evaluation set is entirely `seen`. The analysis is informative about fragmentation sensitivity but does not establish universally held-out species generalization.

## Native CDS output policies

The tools do not share a universal completeness filter. In the installed versions, Tiberius postprocessing removes models below 200 coding bp and models with internal stops. ANNEVO applies native score thresholds, with a stricter threshold for short or single-exon predictions. Carbon-A uses a 96 bp native minimum CDS length, explicitly exports edge-truncated models as `status=incomplete` with `confidence=NA`, and checks for internal stops in the specified reading frame.

An independent audit of all 24 pilot outputs is in `results/native_codon_audit.tsv`. Reproduce it after input reconstruction with `python scripts/check_output_codons.py`. CDS is extracted from GFF plus explicit stop-codon features; translation uses initial phase and the sample code. Literal stop codons and stops resolved from ambiguous IUPAC codons are reported separately. An incomplete trailing codon in a boundary model is not automatically a prediction error.

For S. cerevisiae 1–2 kb inputs, ANNEVO has 15,833/15,833 and Tiberius 1,777/1,777 models passing start/stop/internal-stop checks; Carbon-A has 786/5,256, retaining 4,470 models missing start and/or stop. These checks do not validate reference agreement. The reported partial-CDS-inclusive F1 measures visible sequence recovery; it must not be described as full-length protein accuracy. This audit does not change or filter the benchmark results.

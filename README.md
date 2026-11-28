# Fungal gene calling benchmark

[Release v1.3](https://github.com/WeStrainGroup/fungal_gene_calling_benchmark/releases/tag/v1.3)

## Introduction

We compare four deep learning–based *ab initio* gene predictors—[Helixer](https://github.com/usadellab/Helixer), [ANNEVO](https://github.com/xjtu-omics/ANNEVO), [Tiberius](https://github.com/Gaius-Augustus/Tiberius) and [Carbon-A](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B)—on 17 fungal genomes, from original assemblies to short random fragments. The benchmark evaluates complete CDS for all four tools and truncated CDS separately for Helixer and Carbon-A. Frozen inputs, native predictions and evaluation scripts are provided for reproduction.

## Evaluation

### Tools and models

Versions are frozen across all five input conditions.

| Tool | Tested software | Model |
| --- | --- | --- |
| Helixer | v0.3.8; HelixerPost 0.4.0; GeenuFF 0.4.0 | `fungi_v0.3_a_0100.h5` |
| ANNEVO | v2.3.3, commit [`42c920f`](https://github.com/xjtu-omics/ANNEVO/commit/42c920f70f52f96c7c06c665a639ec47d9304fff) | `ANNEVO_Fungi.pt` |
| Tiberius | 2.0.7 | Official fungi weights, downloaded 2026-09-27 |
| Carbon-A | Carbon-A-1.2B, pipeline/evaluator revision [`ee63599`](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/tree/ee6359990109eef6adccfd19aba8ca64b48ae514) | Published 1.2B checkpoint |

### Genomes and fragmentation

The 17 fungal RefSeq assemblies, accession versions and source datasets are listed in the [species table](metadata/species.tsv). Assemblies shared between source datasets are counted once. Training overlap is documented in [Datasets](docs/datasets.md); these are not uniformly unseen species.

Each genome is evaluated as its original assembly and as full-genome random partitions of 10–100 kb, 5–10 kb, 2–5 kb and 1–2 kb. Cut points are independent of gene positions. All tools receive identical DNA without RNA-seq or protein evidence. Original assemblies, annotations, genetic codes, seeds and fragment maps are included in the [data release](https://github.com/WeStrainGroup/fungal_gene_calling_benchmark/releases/tag/v1.3).

Only records ≥1 kb enter scoring. Shorter original records remain in the release; *A. bisporus* has 771 such records totaling 560,320 bp. Within-genome sequence-length distributions after this cutoff are reported in [Assembly length statistics](docs/input_lengths.md), together with reproduction commands.

### Complete and truncated reference CDS

A CDS is the spliced coding region of one transcript. Original reference CDS are qualified for complete coding structure, valid introns, start/stop codons and reading frame under each genome’s genetic code, then mapped onto the fragments.

| Class | Definition | Comparison |
| --- | --- | --- |
| Complete (C) | All original coding bases remain in the fragment. | Helixer, ANNEVO, Tiberius and Carbon-A. |
| Truncated (T) | Some original coding bases fall outside the fragment. | Helixer and Carbon-A: Overall, missing 5′ end, missing 3′ end and missing both ends. |

Completeness follows reference coordinates, not predicted status. A fragment may contain both classes, and a cut within an intron can truncate a CDS when coding exons lie outside the fragment.

The table summarizes reference gene–fragment instances across all 17 species after qualification, the ≥200 bp visible-CDS cutoff and scoring-region exclusions. Isoforms are not counted separately. One original gene may contribute truncated instances in several fragments, so C + T need not equal the original gene count.

| Input condition | Complete (C) | Truncated (T) |
| --- | --- | --- |
| Original assembly | 192,423 | 0 |
| 10–100 kb | 187,311 | 8,497 |
| 5–10 kb | 153,090 | 65,822 |
| 2–5 kb | 113,078 | 136,415 |
| 1–2 kb | 51,192 | 275,133 |

Per-species counts are available in the [metrics table](results/metrics.tsv), using `level=gene` and `reference`. Complete counts agree across all four tools; truncated counts agree between Helixer and Carbon-A (`stratum=partial`).

### Prediction filters

Released models and native postprocessing are retained. The benchmark applies a visible, spliced CDS cutoff of ≥200 bp to references and predictions, including retained terminal-stop bases and boundary residues. Carbon-A also requires confidence ≥0.1 or explicit NA; missing confidence is an error.

The confidence cutoff follows Carbon-A’s [reference-evaluation rule for genes](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/blob/ee6359990109eef6adccfd19aba8ca64b48ae514/annotation_pipeline/README.md#reference-evaluation). Its native pipeline exports models before this filter and assigns NA to record-edge-truncated genes. We apply the same filtered prediction set at all three scoring levels; upstream applies confidence filtering only to gene/protein metrics.

Helixer’s native partial models are retained; the confidence cutoff applies only to Carbon-A. No additional filter removes predictions for incorrect codons, phase or completeness status. [Run settings](docs/run_predictors.md) document the departures from defaults: Carbon-A’s input-length filter is disabled, its native 96 bp minimum CDS precedes the common cutoff, and Tiberius uses a 1,000 bp minimum input.

### Scoring

Nucleotide and exon scoring reuse the frozen interval definitions from the [Carbon-A technical report](https://www.biorxiv.org/content/10.64898/2026.10.20.760448v1). Gene scoring requires exact coding structure and reading frame.

| Level | Requirement for a true positive |
| --- | --- |
| Nucleotide | A coding base overlaps reference coding sequence on the same record and strand. |
| Exon | A merged coding-exon interval matches both reference boundaries on the same record and strand. These are CDS intervals, not UTR-containing transcript exons. |
| Gene / CDS | The full visible CDS exon-coordinate chain and reading frame match on the same record and strand. For a truncated reference, the entire visible portion and its initial phase must match. |

Nucleotide and exon scoring use the longest eligible visible isoform per gene and merge intervals by record and strand. Gene scoring matches any eligible isoform one-to-one, with at most one TP per reference gene per fragment. Duplicate predictions cannot add TPs, and identical sequence at another locus is not a match. Explicit stop-codon features are normalized with CDS coordinates.

Scoring regions follow reference CDS spans, including internal introns, on each record and strand:

- Complete-reference regions and scoreable background form the main comparison; unmatched predictions count as FP.
- Predictions wholly within truncated-reference regions enter the Helixer and Carbon-A truncated comparison; unmatched predictions count as FP within that scope.
- Failed, too-short or ambiguously classified references and conflicting same-strand overlap regions are excluded. Predictions overlapping excluded regions or crossing truncated-scope boundaries are reported separately, without clipping. Unrecovered scoreable references still count as FN.

Truncated precision and F1 are conditional on truncated-reference regions, not whole-genome estimates. Full isoform, phase and exclusion rules are in [Methods](docs/methods.md); [diagnostics](results/diagnostics.tsv) report background FP and excluded predictions.

## Results

We compare nucleotide-, exon- and gene-level F1 for four tools across 17 fungal genomes and five input conditions. These are descriptive comparisons for the selected species, not evidence of statistical equivalence or performance on unseen species.

### Complete reference CDS

Carbon-A achieves the highest mean F1 at all three levels for original assemblies and 10–100 kb inputs, whereas Tiberius leads at all three levels for 2–5 kb and 1–2 kb inputs. For 5–10 kb inputs, their mean F1 differs by less than 1.6 percentage points at each level. Helixer and ANNEVO have lower mean F1 than both tools in all 15 comparisons.

![Complete reference CDS](figures/complete_cds.png)

Bars show species-equal mean F1 and error bars show sample SD across 17 species. Each genome’s score is calculated from its combined TP/FP/FN across fragments; counts are not pooled across species. SD describes between-species variation. [Summary statistics](results/mean_sd.tsv) and [per-species results](results/metrics.tsv) include precision, recall and F1.

### Truncated reference CDS

Carbon-A achieves higher mean F1 than Helixer at all three levels across all four fragment conditions.

![Helixer and Carbon-A on truncated reference CDS](figures/truncated_cds.png)

Bars show species-equal mean F1 ± sample SD within truncated-reference regions. ANNEVO and Tiberius are marked NA because they are not evaluated in this stratum; NA is not a zero F1 score. Overall uses aggregate truncated-region counts within each genome, including unresolved subtypes, rather than averaging subtype F1. Results by truncation type are shown separately for [Helixer](figures/helixer_partial_cds.png) and [Carbon-A](figures/carbon_partial_cds.png). In those figures, NA means no eligible reference exists; no 10–100 kb case contains a scoreable CDS missing both ends.

### Computational efficiency

The runtime benchmark completed all 1,020 measurements: 17 species × five input conditions × four tools × three independent runs. Each run used one A100 40GB GPU, a four-CPU quota and 32 GiB host memory; physical CPU cores were shared. Original inputs contained all nuclear assembly records ≥1 kb. For each fragment range, 50 frozen records were sampled per species, with ten records from each of five equal-width length bins. The four tools received identical inputs and retained the frozen models and prediction settings. Timing includes model loading, native preprocessing, prediction, postprocessing and output writing, but excludes installation, transfer and validation.

ANNEVO had the highest mean throughput in every input condition, whereas Carbon-A had the lowest. For original assemblies, Carbon-A’s mean elapsed time was about nine times that of ANNEVO and about three times that of Tiberius or Helixer. All four tools achieved higher throughput on longer inputs. Because each fragment set contained 50 records rather than a fixed number of bases, this trend also reflects larger total inputs and the amortization of fixed costs such as model loading.

![Runtime throughput: medians of three runs](figures/efficiency/throughput.png)

Bars show species-equal mean throughput (Mb/min) and sample SD across 17 species. For each species, tool and input condition, throughput is calculated from the median elapsed time of three runs; the mean and SD are then calculated across species. Panel y-axis scales differ. SD describes between-species variation, not variation among the three repeats.

| Tool | Time per original assembly (min, mean ± SD) | Relative elapsed time |
| --- | --- | --- |
| ANNEVO | 1.2 ± 0.7 | 1.0× |
| Tiberius | 3.5 ± 2.0 | 2.9× |
| Helixer | 3.3 ± 2.1 | 2.8× |
| Carbon-A | 11.0 ± 12.5 | 9.1× |

Relative elapsed time is each tool’s mean time divided by the fastest tool’s mean time; 1.0× is the fastest. Values are calculated before rounding and displayed to one decimal place. Three-run medians reduce sensitivity to run-to-run fluctuations, but shared CPU resources still limit interpretation of small timing differences.

## Suggestion

**For complete-CDS prediction in fungal assemblies or individual contigs, we recommend Tiberius for sequences of 1–10 kb and Carbon-A for sequences >10 kb.** Across the tested 1–10 kb ranges, Tiberius achieved similar or higher mean complete-reference F1 than Carbon-A and took approximately one fifth of its mean runtime. Carbon-A achieved the highest mean F1 on longer inputs. Route each sequence by its own length. [CTC: Combined Tiberius and Carbon-A](https://github.com/WeStrainGroup/combined_tiberius_carbon_a) implements this strategy. Both tools apply pretrained models to individual sequence records, allowing CTC to route contigs by length and combine selected complete models. Installation, parameters and output details are available in the [CTC README](https://github.com/WeStrainGroup/combined_tiberius_carbon_a#readme).

**If recovering truncated CDS is a priority and computational cost is less limiting, Carbon-A is the preferred single-tool choice**: it achieved higher mean truncated-reference F1 than Helixer across all tested fragment ranges.

## Reproduce

The published scores can be reproduced on CPU from released native predictions; GPU inference is optional. Run commands from the repository root.

### Get the repository

```bash
git clone --branch v1.3 --depth 1 https://github.com/WeStrainGroup/fungal_gene_calling_benchmark.git
cd fungal_gene_calling_benchmark
Rscript -e 'install.packages("ggplot2", repos="https://cloud.r-project.org")'
```

### Redraw published figures

```bash
Rscript scripts/plot_all_species.R results figures
Rscript scripts/compare_original_assemblies.R results figures
Rscript scripts/plot_efficiency.R results/efficiency/summary.tsv figures/efficiency
```

These commands use the published tables in `results/` and write PNG/PDF files to `figures/`. No data download or rescoring is required.

### Reproduce runtime summaries

The repository includes all 1,020 timing observations. Recalculate each input’s three-run median, then the species means, SDs and relative elapsed times:

```bash
python scripts/summarize_efficiency.py \
  --metrics results/efficiency/round1_metrics.tsv \
            results/efficiency/round2_metrics.tsv \
            results/efficiency/round3_metrics.tsv \
  --expected-rounds 3 --output reproduced_runtime
Rscript scripts/plot_efficiency.R reproduced_runtime/summary.tsv reproduced_runtime/figures
python -m unittest discover -s tests -p 'test_efficiency*.py' -v
```

Python 3.10+ and R/ggplot2 are sufficient; no GPU or model download is needed to reproduce these summaries. [Runtime methods](docs/efficiency.md) also describe reconstructing the 85 exact inputs from the v1.3 data release, configuring locally installed predictors, rerunning timing, and collecting independent rounds. Exact elapsed times depend on hardware and host load.

### Recalculate scores

Create a Python 3.12 environment, install dependencies and run the tests:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
```

Download `fungal_benchmark_data_v1.3.tar.gz` and `SHA256SUMS` from the [data release](https://github.com/WeStrainGroup/fungal_gene_calling_benchmark/releases/tag/v1.3). The compact archive (917 MB) contains all 17 original assemblies, reference annotations, fragment maps, genetic codes and 340 native prediction GFFs. It omits 68 fragment FASTAs that can be restored exactly from the original assemblies and frozen coordinates. The reconstruction step verifies source files, fragment sequences and final FASTA SHA256 values; it does not resample cut points or rerun prediction. Verify the download, reconstruct the fragments and check all extracted files before scoring:

```bash
sha256sum -c SHA256SUMS
tar -xzf fungal_benchmark_data_v1.3.tar.gz
python scripts/reconstruct_fragments.py --data data
(cd data && sha256sum -c SHA256SUMS)
python scripts/evaluate_stratified.py --data data \
  --samples data/samples.tsv --cohort fungal_benchmark --output data/evaluation
python scripts/validate_stratified.py --data data --output data/evaluation
```

Summarize the validated results and plot:

```bash
python scripts/summarize_all_species.py --inputs data/evaluation \
  --species data/samples.tsv --output combined_results
Rscript scripts/plot_all_species.R combined_results combined_results/figures
Rscript scripts/validate_summary.R combined_results
```

Recalculated tables, validation reports and figures are written under `combined_results/`; published tables remain in `results/`. The main outputs are `metrics.tsv` (per-species counts and rates), `mean_sd.tsv` (species-equal summaries) and `diagnostics.tsv` (filtering and exclusions).

To run models on new inputs, see [Running the predictors](docs/run_predictors.md). The [script registry](docs/methods.md#implementation) and [validation evidence](results/validation.json) document the maintained implementation and completed checks.

## Contact

Xinyu Wang: [wangxinyu30@westlake.edu.cn](mailto:wangxinyu30@westlake.edu.cn)

## License

Copyright © 2026 WeStrainGroup.

[CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/) for our benchmark scripts and documentation. The vendored Carbon-A evaluator retains its Apache 2.0 license; external predictors, model weights and reference datasets retain their own terms.

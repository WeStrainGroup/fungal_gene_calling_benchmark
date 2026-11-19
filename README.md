# Fungal gene calling benchmark

This benchmark compares ANNEVO, Tiberius and Carbon-A on 12 fungal species drawn from the three tools’ test datasets. For each tool and fragmentation condition, we report the mean ± SD of each metric across species, with equal weight per species.

Conditions comprise the **original assembly**, **1–2 kb**, **5–10 kb** and **50–100 kb** full-genome random partitions. All tools receive identical DNA; no RNA or protein evidence is supplied. Saccharomyces cerevisiae is counted once.

| Species | Assembly | Source |
| --- | --- | --- |
| *Saccharomyces cerevisiae* | GCF_000146045.2 | ANNEVO, Carbon-A |
| *Candida albicans* | GCF_000182965.3 | ANNEVO |
| *Agaricus bisporus* | GCF_000300555.1 | Tiberius |
| *Cryphonectria parasitica* | GCF_011745365.1 | Tiberius |
| *Parastagonospora nodorum* | GCF_000146915.1 | Tiberius |
| *Puccinia striiformis* | GCF_021901695.1 | Tiberius |
| *Punctularia strigosozonata* | GCF_000264995.1 | Tiberius |
| *Tilletiopsis washingtonensis* | GCF_003144115.1 | Tiberius |
| *Schizosaccharomyces pombe* | GCF_000002945.2 | Carbon-A |
| *Neurospora crassa* | GCF_000182925.2 | Carbon-A |
| *Coprinopsis cinerea* | GCF_000182895.1 | Carbon-A |
| *Rhizophagus irregularis* | GCF_026210795.1 | Carbon-A |

## Evaluation

The v0.2 scoring protocol is unchanged: compare all three tools on complete reference CDS, and assess Carbon-A separately on truncated reference CDS. Gene/CDS matches require identical visible CDS coordinates, strand and reading frame. Main scoring uses input records >=1 kb and visible CDS >=200 bp; Carbon-A retains confidence >=0.1 or explicit NA.

Metrics are first computed per species. The final summary combines all 12 species using the arithmetic mean and **sample SD (n−1)** of defined values. Missing values are never replaced with zero; each row reports `n_species`. SD is NA for fewer than two defined species. Original assemblies are a separate condition. This replaces the earlier source-specific means without SD.

## Interim benchmark

Updated **2026-11-19**. This release includes **33 species × condition pairs (99 prediction cases)** for which all three tools have completed inference and independent validation. It covers 11 distinct species; only completed pairs enter each condition. The final update will include all 12 species.

Complete-reference gene/CDS F1 (%), species-equal mean ± sample SD:

| Condition | Species | ANNEVO | Tiberius | Carbon-A |
| --- | ---: | ---: | ---: | ---: |
| Original assembly | 7 | 57.58 ± 20.93 | 59.67 ± 24.33 | 73.58 ± 18.60 |
| 50–100 kb | 8 | 55.41 ± 20.10 | 57.60 ± 22.80 | 71.35 ± 18.24 |
| 5–10 kb | 10 | 56.48 ± 20.67 | 61.32 ± 22.70 | 63.27 ± 21.17 |
| 1–2 kb | 8 | 24.88 ± 15.78 | 63.10 ± 19.66 | 43.83 ± 17.66 |

![Complete CDS: interim benchmark](figures/complete_cds.png)

**Species composition differs between conditions**, so these cross-condition means are not a paired estimate of fragmentation loss. Within each condition, all three tools use the same species. The species included in each condition are listed in [STATUS.json](results/STATUS.json). Undefined metric values are excluded, with their defined-species counts reported in the summary table.

[Mean ± SD: all metrics](results/mean_sd.tsv) · [Per-species TP/FP/FN and P/R/F1](results/metrics.tsv) · [Background and exclusions](results/diagnostics.tsv) · [Truncated CDS figure](figures/carbon_partial_cds.pdf) · [Independent validation](results/validation.json)

[Methods](docs/methods.md) · [Dataset details](docs/datasets.md) · [Script registry](docs/scripts.md)

The final CPU follow-up remains active. After all native predictions and source validations finish, it will generate the complete 12-species summary for replacement of this interim report.

## Reproduce

Install the Python dependencies in `requirements.txt`; plotting requires R packages `ggplot2` and `jsonlite`. The [current data release](https://github.com/WeStrainGroup/fungal_gene_calling_benchmark/releases/tag/v0.3.2) contains frozen pilot inputs and native predictions, with checksums. Expanded inputs use the accession manifests in `metadata/`.

```bash
python scripts/fetch_data.py
python scripts/prepare_inputs.py
python -m unittest discover -s tests -v
python scripts/evaluate_stratified.py --data data --output data/evaluation_v0.2 --cohort annevo_set
python scripts/validate_stratified.py --data data --output data/evaluation_v0.2
```

For an interim report using completed species-condition pairs:

```bash
python scripts/summarize_all_species.py \
  --inputs annevo_results tiberius_results carbon_results \
  --output combined_results --paired-conditions
Rscript scripts/plot_all_species.R combined_results combined_results/figures
```

For the final report, replace `--paired-conditions` with `--wait --rscript Rscript`. This waits for all source validations and produces the complete 12-species summary and plots. The repository contains only the current analysis. Reusable inputs, native predictions and validation records are retained; obsolete scores, figures and temporary experiments are removed.

## Sources and license

[ANNEVO](https://github.com/xjtu-omics/ANNEVO) · [Tiberius](https://github.com/Gaius-Augustus/Tiberius) · [Carbon-A](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B) · [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/)

Benchmark code: [MIT](LICENSE). Vendored evaluator: [Apache 2.0](vendor/carbon_a/LICENSE). See [third-party notices](THIRD_PARTY.md), versioned metadata and `CITATION.cff`.

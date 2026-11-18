# Fungal gene calling benchmark

We benchmarked ANNEVO, Tiberius and Carbon-A on 12 fungal species drawn from the three tools’ test datasets. For each tool and fragmentation condition, we report the mean ± SD of each metric across species, with equal weight per species.

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

## Status and results

The unified snapshot currently contains **81/144 unique prediction cases**. The remaining predictions and their CPU validation are pending. Final 12-species mean ± SD tables and figures are withheld until all cases have been validated.

[Per-species metrics](results/metrics.tsv) · [Background and exclusions](results/diagnostics.tsv) · [Completion status](results/STATUS.json) · [Methods](docs/methods.md) · [Dataset details](docs/datasets.md) · [Script registry](docs/scripts.md)

The locked CPU follow-up combines the three validated input directories, verifies that duplicate yeast rows agree, and generates `mean_sd.tsv` and R figures after completion. It does not rerun inference.

## Reproduce

Install the Python dependencies in `requirements.txt`; plotting requires R packages `ggplot2` and `jsonlite`. The [current data release](https://github.com/WeStrainGroup/fungal_gene_calling_benchmark/releases/tag/v0.3.0) contains frozen pilot inputs and native predictions, with checksums. Expanded inputs use the accession manifests in `metadata/`.

```bash
python scripts/fetch_data.py
python scripts/prepare_inputs.py
python -m unittest discover -s tests -v
python scripts/evaluate_stratified.py --data data --output data/evaluation_v0.2 --cohort annevo_set
python scripts/validate_stratified.py --data data --output data/evaluation_v0.2
```

After independently validating all three source datasets:

```bash
python scripts/summarize_all_species.py \
  --inputs annevo_results tiberius_results carbon_results \
  --output combined_results
Rscript scripts/plot_all_species.R combined_results combined_results/figures
```

Use `--wait --rscript Rscript` with live result directories to wait for all three validations and then produce the final summary and plots. The repository contains only the current analysis. Reusable inputs, native predictions and validation records are retained; obsolete scores, figures and temporary experiments are removed.

## Sources and license

[ANNEVO](https://github.com/xjtu-omics/ANNEVO) · [Tiberius](https://github.com/Gaius-Augustus/Tiberius) · [Carbon-A](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B) · [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/)

Benchmark code: [MIT](LICENSE). Vendored evaluator: [Apache 2.0](vendor/carbon_a/LICENSE). See [third-party notices](THIRD_PARTY.md), versioned metadata and `CITATION.cff`.

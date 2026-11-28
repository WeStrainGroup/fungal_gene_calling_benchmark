# Runtime benchmark

The current release contains all three rounds: 1,020 successful measurements covering 17 species × five input conditions × four tools. Each input is represented by its median elapsed time across three independent runs, followed by species-equal mean and sample SD across 17 species.

## Reproduce the current tables and figure

Run from the repository root with Python 3.10+ and R with ggplot2. No GPU, model download or external Python package is needed for these summaries.

```bash
python scripts/summarize_efficiency.py \
  --metrics results/efficiency/round1_metrics.tsv \
            results/efficiency/round2_metrics.tsv \
            results/efficiency/round3_metrics.tsv \
  --expected-rounds 3 --output reproduced_runtime
Rscript scripts/plot_efficiency.R \
  reproduced_runtime/summary.tsv reproduced_runtime/figures
python -m unittest discover -s tests -p 'test_efficiency*.py' -v
```

Each source table records one successful observation per species/tool/condition in that round: input size, elapsed seconds, CPU seconds, maximum individual-process RSS, execution order, task ID, container hostname, host load and CPU throttling. The script checks all 340 combinations in each of the three rounds against the frozen input manifest and rejects failed, missing, duplicated or non-finite observations, including reused task IDs across rounds. `species_medians.tsv` contains the 340 within-input medians; `summary.tsv` contains the resulting mean and sample SD across 17 species. `assembly_time.md` is the one-decimal presentation table; TSVs retain full precision. Relative elapsed time is each tool's mean assembly time divided by the fastest tool's mean assembly time; the fastest is 1.0×.

The throughput figure uses Original → 10–100 kb → 5–10 kb → 2–5 kb → 1–2 kb, with independent y-axis scales and the established tool colors. SD describes variation between species, not uncertainty from repeated runs.

## Reconstruct the exact inputs

Use the [v1.3 data release](https://github.com/WeStrainGroup/fungal_gene_calling_benchmark/releases/tag/v1.3), shared with the accuracy benchmark. The compact archive contains frozen assemblies, coordinate maps and genetic codes. Restore the omitted fragment FASTAs from these maps and verify their original checksums before preparing runtime inputs; no new random cutting or reference download is needed.

```bash
mkdir -p data_download
curl -L --fail -o data_download/fungal_benchmark_data_v1.3.tar.gz \
  https://github.com/WeStrainGroup/fungal_gene_calling_benchmark/releases/download/v1.3/fungal_benchmark_data_v1.3.tar.gz
curl -L --fail -o data_download/SHA256SUMS \
  https://github.com/WeStrainGroup/fungal_gene_calling_benchmark/releases/download/v1.3/SHA256SUMS
(cd data_download && sha256sum -c SHA256SUMS)
tar -xzf data_download/fungal_benchmark_data_v1.3.tar.gz
python scripts/reconstruct_fragments.py --data data
(cd data && sha256sum -c SHA256SUMS)
python scripts/prepare_efficiency.py --data data --output runtime_work
```

- Original: all nuclear assembly records ≥1 kb, separately for each species.
- Fragments: 50 records per species per range, sampled without replacement from existing frozen fragments. Each range has five equal-width length bins, with ten records selected per bin; seed `20261124`.
- The sampling script writes `runtime_work/cases.tsv`, `selected_sequences.tsv` and `inputs/<species>/`. Compare the cases to `metadata/efficiency_inputs.tsv`; `metadata/efficiency_sampled_records.tsv` records all 3,400 selected fragment IDs and sequence hashes. Absolute source paths in the generated selection log are machine-specific and not part of sequence identity.

| Input condition | Species | Records | bp |
| --- | ---: | ---: | ---: |
| Original | 17 | 2,664 | 690,167,379 |
| 10–100 kb | 17 | 850 | 46,505,838 |
| 5–10 kb | 17 | 850 | 6,376,270 |
| 2–5 kb | 17 | 850 | 2,972,880 |
| 1–2 kb | 17 | 850 | 1,277,565 |

## Rerun timing on a GPU

The measured allocation was one A100-SXM4-40GB, 4 CPU quota, 32 GiB host RAM, with an eight-hour limit per species task. AIS CPU cores, caches and memory bandwidth were shared. Each task runs five conditions × four tools serially. Allocate these resources through your scheduler; the runner verifies the actual device, BF16 support, memory limit and CPU allocation before timing.

Install the frozen versions and obtain the authors' models following [Running the predictors](run_predictors.md). Package versions, source revisions and model SHA256 values are in `metadata/software.json`; the hardware block there describes the original accuracy runs, not this A100 runtime experiment. No models or proprietary platform image are bundled. The original AIS base image was `192.168.108.1:5000/pytorch/pytorch2.3-admin:cuda_12.4`; this private registry is not required for CPU summary reproduction or the configurable runner.

Copy and edit the path configuration. Paths may be absolute or relative to the JSON file. Each environment must have its executables in `bin/`, including Helixer.py and HelixerPost for Helixer. Use GNU Time, not the BSD/macOS `time` command.

```bash
cp metadata/efficiency_config.example.json runtime_work/config.json
# Edit runtime_work/config.json to match your installed environments and models.
# Run inside a scheduler allocation; use the Carbon-A Python environment for the runner.
export BENCHMARK_RUN_ID=runtime_round1_species0
/path/to/envs/carbon/bin/python scripts/run_efficiency.py \
  --platform ais --config runtime_work/config.json \
  --inputs runtime_work/inputs --cases runtime_work/cases.tsv \
  --metadata metadata/software.json \
  --output runtime_work/results/round1 --species-index 0
```

Use `--platform ais` for a cgroup CPU quota of four, or `--platform slurm` for four bound physical CPU cores. Export a unique `BENCHMARK_RUN_ID` per species/round outside Slurm; Slurm supplies `SLURM_JOB_ID`. Species indices are the alphabetical order of the 17 species in `cases.tsv`, from 0 to 16. Set a new output round for each repetition. Existing species output directories are refused. The original AIS wrapper, `efficiency_ais.sh WORK SPECIES_INDEX ROUND`, retains our frozen server layout; external installations should use `--config` above. The configurable path option does not change predictor commands or timing boundaries.

The five conditions are shuffled with seed `20261124 + species_index`, and tool order rotates by species and condition. All three rounds use the same within-species schedule. Timing includes model loading, native preprocessing, prediction, postprocessing and output writing through `run_predictor.sh`; input/model checks, installation, transfer, resource checks and scratch cleanup are excluded. Host page caches are not flushed, and native compilation caches remain enabled.

For each species, outputs include `metrics.tsv`, `order.tsv`, `hardware.json`, `software.json`, and native prediction/log/timing files under `<tool>/<condition>/`. GNU-time RSS is the maximum individual-process RSS, not aggregate process-tree peak memory or GPU memory. Host load and CPU throttling are diagnostics; throttling can reflect a task reaching its own quota and does not by itself establish interference. Exact elapsed times are not expected to reproduce across machines or host load.

## Collect and summarize three rounds

```bash
for round in 1 2 3; do
  python scripts/collect_efficiency.py \
    --results runtime_work/results/round${round} \
    --output runtime_work/round${round}_metrics.tsv
done
python scripts/summarize_efficiency.py \
  --metrics runtime_work/round1_metrics.tsv runtime_work/round2_metrics.tsv runtime_work/round3_metrics.tsv \
  --output reproduced_runtime
Rscript scripts/plot_efficiency.R reproduced_runtime/summary.tsv reproduced_runtime/figures
```

All three rounds passed collection, with 51 passing hardware receipts, consistent software/model checksums and 1,020 nonempty native GFFs. The retrieved first-round table matches the previously published table byte for byte.

Collection requires 20 successful cases for each of 17 species, passing hardware receipts and nonempty native GFFs; it saves GFF hashes alongside the combined table. The final summary requires 1,020 successful observations with independent task IDs. For each species/tool/condition, take the median of its three elapsed times, calculate Mb/min from that median, then report the mean and sample SD across 17 species (ddof=1). Each species has equal weight. Do not pool 51 runs as independent species, substitute missing runs with zero, or combine HPC preflight timings with AIS results. Peak process RSS uses the same median-then-mean convention. `species_medians.tsv` retains the 340 within-input values.

## CTC timing

CTC code, tests and usage are maintained in the [standalone repository](https://github.com/WeStrainGroup/combined_tiberius_carbon_a). Clone that repository and pass its `ctc.py` path to `run_hybrid_efficiency.py --ctc-script /path/to/combined_tiberius_carbon_a/ctc.py`; the AIS launcher reads this path from `CTC_SCRIPT`. Timing records identify the CTC version and confidence threshold; existing version-isolated jobs are unchanged. Current CTC timing includes native protein export, confidence filtering, ID mapping and unified translation of selected complete CDS. Existing measurements from earlier CTC versions must not be relabelled as current-version timings.

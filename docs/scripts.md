# Scripts and status

| Script | Purpose | Status |
| --- | --- | --- |
| `summarize_all_species.py` | Current final aggregation: 12 unique species, equal-weight mean and sample SD; optional locked wait for final validation. | Deduplication, NA/zero handling and sample SD tested; final summary waiting for inference. |
| `plot_all_species.R` | Current combined complete/truncated F1 figures with sample-SD error bars. | Uses the final combined summary; preserves the accepted colors and theme. |
| `tests/test_species_summary.py` | Duplicate species, sample SD, undefined values and missing species checks. | Six tests passed. |
| `evaluate_stratified.py` | v0.2 truth strata, common filters and coordinate/frame matching; all metrics and exclusions. | Pilot 24/24 cases verified; extension snapshots/follow-ups kept separate. |
| `validate_stratified.py` | Independent raw-GFF gene TP and count reconciliation plus source rehashing. | Pilot PASS; 106,556 assigned gene matches checked. |
| `plot_stratified.R` | Historical source-specific v0.2.0 complete and truncated F1 plots. | Retained for reproducibility; current final figures use `plot_all_species.R`. |
| `finish_stratified.py` | Locked CPU snapshot/wait/evaluate/validate/plot follow-up for live native outputs. | Separate extension workers; no inference or automatic publication. |
| `tests/test_stratified.py` | Scientific boundaries, matching and output-grid regression tests. | 24 tests passed. |
| `fetch_data.py` | Download or read the release archive, verify SHA256 and extract regular files. | Validated with the release archive. |
| `prepare_inputs.py` | Reconstruct the exact completed inputs from frozen maps and public RefSeq DNA. | All inference-time input hashes reproduced. |
| `evaluate.py` | Historical v0.1 reference and sequence-matching evaluator; shared frozen helpers for v0.2. | All 72 historical metric rows reproduced; retained unchanged. |
| `check_output_codons.py` | Audit pilot predicted CDS start/stop, initial phase, internal stops and ambiguity without filtering. | All 24 pilot cases checked. |
| `compare_inputs.py` | Compare original-genome and single-record annotations, ignoring generated IDs. | All four initial tool/species comparisons reproduced. |
| `run_predictor.sh` | Run one input in an already installed environment with the tested scientific settings. | Adapted from completed inference; used by the six-species extension. |
| `prepare_new_inputs.py` | Download fixed RefSeq versions and generate complete nuclear partitions from a sample manifest. | Six extension inputs downloaded and verified. |
| `summarize.py` | Historical single-source mean P/R/F1 summary. | Retained unchanged; superseded by `summarize_all_species.py` for final reporting. |
| `plot.R` | Historical mean-F1 bars for one benchmark set. | Retained unchanged; current combined figures include sample SD. |
| `plot_pilot.R` | Reproduce the initial per-species figure. | Completed. |

`vendor/carbon_a/evaluation.py`, `utils.py` and `constants.py` are frozen upstream dependencies, not independent analysis entry points. Their original Apache 2.0 license is included. `metadata/` identifies versions, sources and models; `results/` contains the completed outputs. Status applies to this release, not to an assertion that queued extension results already exist.

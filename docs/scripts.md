# Scripts and status

| Script | Purpose | Status |
| --- | --- | --- |
| `fetch_data.py` | Fetch and verify the current frozen pilot input/prediction bundle. | Current archive validated. |
| `prepare_inputs.py` | Reconstruct pilot inputs from frozen source records and fragment maps. | Input hashes checked. |
| `prepare_new_inputs.py` | Prepare fixed RefSeq accessions and whole-genome partitions. | Twelve-species inputs prepared. |
| `run_predictor.sh` | Run one native predictor on one input. | Existing outputs retained; expanded inference running. |
| `common.py` | Shared FASTA helpers and the pinned interval evaluator. | Used by the current evaluator. |
| `evaluate_stratified.py` | Reference strata, common prediction filters, structural matching and metrics. | Pilot and available extension cases verified. |
| `validate_stratified.py` | Independent coordinate/phase and count checks; rehash inputs. | Latest refresh: 66 Tiberius-set and 31 Carbon-A-set cases verified; pilot 24 cases retained. |
| `finish_stratified.py` | Wait for native outputs, then evaluate and validate a source dataset. | CPU follow-ups active; no source-specific summary plots. |
| `summarize_all_species.py` | Deduplicate species and calculate mean/sample SD; wait for final validation. | Interim report: 33 paired species-condition units, 99 cases; final summary waiting. |
| `plot_all_species.R` | Final complete/truncated CDS figures with sample-SD error bars. | Interim figures validated; actual species counts labeled per condition. |
| `tests/test_stratified.py` | Scientific scoring boundary tests. | 24 tests passed. |
| `tests/test_species_summary.py` | Species weighting, deduplication, NA and sample SD. | Seven tests passed. |

The unchanged files in `vendor/carbon_a/` provide the pinned interval definitions and retain their license. Runtime environments, primary inputs and native predictions are reusable assets, not disposable intermediates. Server scheduling is operational: two existing jobs continue, and the extra queue waits for shared locks on GPUs 2 and 3; after the six-species inference finishes it may use GPUs 0–3, up to four jobs total. GPU/CPU scientific parameters remain unchanged.

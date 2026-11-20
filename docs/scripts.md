# Scripts and status

| Script | Purpose | Status |
| --- | --- | --- |
| `fetch_data.py` | Fetch and verify the current frozen pilot input/prediction bundle. | Current archive validated. |
| `prepare_inputs.py` | Reconstruct pilot inputs from frozen source records and fragment maps. | Input hashes checked. |
| `prepare_new_inputs.py` | Prepare fixed RefSeq accessions and whole-genome partitions. | Twelve-species inputs prepared. |
| `run_queue.py` | Schedule both expansion datasets on shared GPUs 0–7; adopt existing processes and skip completed cases. | Current: one controller, one remaining Carbon-A case on GPU 2; no queued or failed cases. |
| `run_predictor.sh` | Run one native predictor on one input. | Existing outputs retained; expanded inference running. |
| `common.py` | Shared FASTA helpers and the pinned interval evaluator. | Used by the current evaluator. |
| `evaluate_stratified.py` | Reference strata, common prediction filters, structural matching and metrics. | Pilot and available extension cases verified. |
| `validate_stratified.py` | Independent coordinate/phase and count checks; rehash inputs. | Latest snapshot: all 72 Tiberius-set and 59 Carbon-A-set cases verified; pilot 24 cases retained. |
| `finish_stratified.py` | Wait for native outputs, then evaluate and validate a source dataset. | Tiberius-set validation complete; Carbon-A follow-up awaits the last native case. |
| `summarize_all_species.py` | Deduplicate species and calculate mean/sample SD; wait for final validation. | Interim report: 47 paired species-condition units, 141 cases; final summary waiting. |
| `plot_all_species.R` | Final complete/truncated CDS figures with sample-SD error bars. | Interim figures validated; actual species counts labeled per condition. |
| `tests/test_stratified.py` | Scientific scoring boundary tests. | 24 tests passed. |
| `tests/test_species_summary.py` | Species weighting, deduplication, NA and sample SD. | Seven tests passed. |

The unchanged files in `vendor/carbon_a/` provide the pinned interval definitions and retain their license. Runtime environments, primary inputs and native predictions are reusable assets, not disposable intermediates. Server scheduling uses one controller for both expansion datasets, up to eight simultaneous cases on GPUs 0–7. Shared GPU locks and a fresh occupancy check are required before launch. Existing inference processes retain their locks during adoption; their exit status and nonempty native GFF are checked before completion. Per-case scientific parameters remain unchanged.

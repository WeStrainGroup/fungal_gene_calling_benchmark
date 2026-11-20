# Scripts and status

| Script | Purpose | Status |
| --- | --- | --- |
| `fetch_data.py` | Fetch and verify the current frozen pilot input/prediction bundle. | Current archive validated. |
| `prepare_inputs.py` | Reconstruct pilot inputs from frozen source records and fragment maps. | Input hashes checked. |
| `prepare_new_inputs.py` | Prepare fixed RefSeq accessions and whole-genome partitions. | Twelve-species inputs prepared. |
| `run_carbon_multigpu.py` | Reserve eight shared GPUs, replace the final single-GPU inference and run the unchanged native pipeline with GPU_COUNT=8. | Active for Rhizophagus 1–2 kb; prior in-memory progress could not be resumed, old log retained. |
| `run_queue.py` | Schedule both expansion datasets on shared GPUs 0–7; adopt existing processes and skip completed cases. | Completed all other cases; retired for the final native eight-GPU case. |
| `run_predictor.sh` | Run one native predictor on one input. | Existing outputs retained; expanded inference running. |
| `common.py` | Shared FASTA helpers and the pinned interval evaluator. | Used by the current evaluator. |
| `evaluate_stratified.py` | Reference strata, common prediction filters, structural matching and metrics. | Pilot and available extension cases verified. |
| `validate_stratified.py` | Independent coordinate/phase and count checks; rehash inputs. | Latest snapshot: all 72 Tiberius-set and 59 Carbon-A-set cases verified; pilot 24 cases retained. |
| `finish_stratified.py` | Wait for native outputs, then evaluate and validate a source dataset. | Tiberius-set validation complete; Carbon-A follow-up awaits the last native case. |
| `summarize_all_species.py` | Deduplicate species and calculate mean/sample SD; wait for final validation. | Interim report: 47 paired species-condition units, 141 cases; final summary waiting. |
| `plot_all_species.R` | Final complete/truncated CDS figures with sample-SD error bars. | Interim figures validated; actual species counts labeled per condition. |
| `tests/test_stratified.py` | Scientific scoring boundary tests. | 24 tests passed. |
| `tests/test_species_summary.py` | Species weighting, deduplication, NA and sample SD. | Seven tests passed. |

The unchanged files in `vendor/carbon_a/` provide the pinned interval definitions and retain their license. Runtime environments, primary inputs and native predictions are reusable assets, not disposable intermediates. All other cases have completed. The final Carbon-A Rhizophagus 1–2 kb case now uses native eight-GPU dynamic scheduling (dispatch chunk size 4, batch size 4, SDPA). The official implementation writes probabilities only after all input records finish, so this case restarted from the same full input; the earlier single-GPU log and handoff record are retained. The supervisor holds all eight shared locks and marks COMPLETE only after the native pipeline exits successfully and all three output files exist. Model, input, scoring thresholds and post-processing parameters are unchanged.

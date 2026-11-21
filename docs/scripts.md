# Scripts and status

| Script | Purpose | Status |
| --- | --- | --- |
| `fetch_data.py` | Download and verify all 12 frozen inputs and 144 native GFFs. | All archived files verified. |
| `prepare_inputs.py` | Reconstruct the two pilot inputs from source records and maps. | Pilot input hashes checked; release already includes inputs. |
| `prepare_new_inputs.py` | Prepare additional RefSeq accessions and whole-genome partitions. | Expansion inputs frozen; not needed to score released data. |
| `run_predictor.sh` | Run one native predictor on one input; Carbon-A supports GPU_COUNT. | 144 unique inference cases complete. |
| `common.py` | FASTA helpers and pinned interval evaluator. | Current evaluator dependency. |
| `evaluate_stratified.py` | Reference strata, common filters, structural matching and metrics. | All three source datasets complete. |
| `validate_stratified.py` | Independent raw coordinate/phase and count checks; source hashes. | 24 + 72 + 60 source cases passed; 144 unique cases. |
| `summarize_all_species.py` | Deduplicate shared species and calculate mean/sample SD. | Final 1,008 metric rows; 252 summary rows independently checked in R. |
| `plot_all_species.R` | Complete and truncated CDS figures. | Final twelve-species figures inspected. |
| `tests/test_stratified.py` | Scoring boundary tests. | 24 tests passed. |
| `tests/test_species_summary.py` | Weighting, deduplication, NA and sample SD tests. | Seven tests passed. |

One-time GPU scheduling and completion-waiting helpers have been removed. The final Carbon-A Rhizophagus 1–2 kb case used the official native pipeline on eight GPUs, with the same model, inputs, thresholds and postprocessing. Other inference results were reused. Vendored interval definitions retain their original license.

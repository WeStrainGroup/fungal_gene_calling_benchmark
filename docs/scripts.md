# Scripts and status

| Script | Purpose | Status |
| --- | --- | --- |
| `fetch_data.py` | Download or read the release archive, verify SHA256 and extract regular files. | Validated with the release archive. |
| `prepare_inputs.py` | Reconstruct the exact completed inputs from frozen maps and public RefSeq DNA. | All inference-time input hashes reproduced. |
| `evaluate.py` | Rebuild references and calculate nucleotide/exon/gene metrics using the frozen official evaluator. | All 72 initial metric rows reproduced; independent matcher and official entry point agree. Supports additional sample manifests. |
| `compare_inputs.py` | Compare original-genome and single-record annotations, ignoring generated IDs. | All four initial tool/species comparisons reproduced. |
| `run_predictor.sh` | Run one input in an already installed environment with the tested scientific settings. | Adapted from completed inference; used by the six-species extension. |
| `prepare_new_inputs.py` | Download fixed RefSeq versions and generate complete nuclear partitions from a sample manifest. | Six extension inputs downloaded and verified. |
| `summarize.py` | Species-equal mean P/R/F1 for one complete set; no SD or error bars. | Requires a complete species × tool × length × level grid. |
| `plot.R` | R/ggplot2 mean-F1 bars for one benchmark set. | Uses the accepted muted colors and theme_minimal. |
| `plot_pilot.R` | Reproduce the initial per-species figure. | Completed. |

`vendor/carbon_a/evaluation.py`, `utils.py` and `constants.py` are frozen upstream dependencies, not independent analysis entry points. Their original Apache 2.0 license is included. `metadata/` identifies versions, sources and models; `results/` contains the completed outputs. Status applies to this release, not to an assertion that queued extension results already exist.

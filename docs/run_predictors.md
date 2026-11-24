# Running the predictors

Evaluation reuses native predictions without GPU inference. These commands describe the tested scientific parameters for running new inputs. Install each tool in its own environment using its upstream instructions; the wrapper does not install software, acquire a GPU, or change a scheduler.

## Versions and models

| Tool | Tested version/source | Model |
| --- | --- | --- |
| ANNEVO | 2.3.3; commit `42c920f70f52f96c7c06c665a639ec47d9304fff` | `ANNEVO_Fungi.pt` |
| Tiberius | 2.0.7; bricks2marble 0.1.1; TensorFlow 2.20.0 | Official fungi weights |
| Helixer | 0.3.8; HelixerPost 0.4.0; GeenuFF 0.4.0; TensorFlow 2.15.1 | `fungi_v0.3_a_0100.h5` |
| Carbon-A | Published pipeline/evaluator revision `ee6359990109eef6adccfd19aba8ca64b48ae514`; PyTorch 2.13.0+cu126; Transformers 4.57.6 | Carbon-A 1.2B |

Model SHA256 values and installed versions are in `metadata/software.json`. The Tiberius model URL and file hashes are in `metadata/tiberius_model_source.json`. Source fingerprints are in `metadata/source_sha256.json`. The Carbon-A weights originated from the authors' earlier GenerTeam preview and match the published model hash. Its model-loading compatibility changes are author-supplied: the local `modeling_generanno.py` SHA256 (`d4ee34258c1a3bdad43ea3c3d687e4f210c758790859e0d671033b181809e0ef`) matches the official `provenance.json`; this is not an additional benchmark patch. The public evaluator is unmodified. Full inference is not claimed bit-identical across software revisions, kernels or hardware.

Version check (2026-11-23): ANNEVO matches the current official main commit. Tiberius 2.0.7 is older than the latest release, 2.0.9 (2026-10-05); no equivalence of predictions is assumed. Carbon-A matches the revision shown by the official repository and its published checkpoint hash. Exact tested revisions remain frozen across all input conditions.

Official installation sources: [ANNEVO](https://github.com/xjtu-omics/ANNEVO), [Tiberius](https://github.com/Gaius-Augustus/Tiberius), [Carbon-A](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B). Obtain models from the authors and verify their hashes. No model weights are bundled here.

## One input per invocation

Select a genuinely free GPU through your scheduler or resource-management system. `CUDA_VISIBLE_DEVICES` alone does not reserve a device. Activate the corresponding environment before each command, or set `PREDICTOR_PYTHON` explicitly. Carbon-A's pipeline invokes Python through `PATH`.

```bash
export CUDA_VISIBLE_DEVICES=0
input=data/inputs/Saccharomyces_cerevisiae/1-2kb.fna
code=data/inputs/Saccharomyces_cerevisiae/codon_table.json

# In the ANNEVO environment:
export ANNEVO_DIR=/path/to/ANNEVO
bash scripts/run_predictor.sh annevo "$input" data/new_predictions/annevo "$code"

# In the Tiberius environment:
export TIBERIUS_MODEL=/path/to/tiberius/fungi
bash scripts/run_predictor.sh tiberius "$input" data/new_predictions/tiberius "$code"

# In the Carbon-A environment:
export CARBON_DIR=/path/to/Carbon-A
export GPU_COUNT=1
bash scripts/run_predictor.sh carbon "$input" data/new_predictions/carbon "$code"
```

Reuse successful native predictions; each pending case has one output directory. For evaluation, the layout is `out/<tool>/<species>/<group>/prediction.gff` for Carbon-A and `prediction.gff3` for the other tools. Current groups are `genome`, `10-100kb`, `5-10kb`, `2-5kb`, `1-2kb`. Only the current five conditions enter the maintained workflow.

For Helixer, activate its isolated environment and run:

```bash
export HELIXER_MODEL=/path/to/fungi_v0.3_a_0100.h5
bash scripts/run_predictor.sh helixer "$input" data/new_predictions/helixer "$code"
```

## Fixed settings

- Carbon-A: batch size 4, SDPA attention, CPU Viterbi, minimum record length 0, native minimum CDS 96 bp, explicit sample codon table. Allocated GPU counts are recorded per case; the current queue schedules one case per GPU, with long Carbon-A cases started early. The pipeline distributes records over devices. SDPA does not produce the authors' FlashAttention2-identical probabilities.
- ANNEVO: Fungi model, batch size 32, 4 threads, 2 workers, native decoding thresholds, no `--overlap_pred`.
- Tiberius: fungi model, batch size 8, minimum input 1,000 bp, native postprocessing. Native protein exports may use the standard code; common evaluation independently extracts CDS from the input and GFF.
- Helixer: fungi model, subsequence length 21,384 bp, batch size 32, overlap and all native postprocessing defaults (HelixerPost 0.4.0). Explicit model path prevents silent model updates. Preserve partial predictions; no completeness filter.
- Native predictions remain unchanged. The common evaluator applies visible CDS >=200 bp and Carbon-A confidence >=0.1 or explicit NA, without filtering by predicted completeness.

Native run resources and logs record allocated devices and wall-clock execution. These are not hardware-normalized speed comparisons or aggregate peak process-tree memory measurements.

## Compute orchestration on the original server

```bash
base=/ssd/wangxinyu/carbon_a_benchmark_20261117
data=$base/helixer_extension
bash scripts/install_helixer.sh "$data"
python scripts/prepare_new_inputs.py --data "$data" --samples metadata/helixer_samples.tsv
python scripts/prepare_benchmark.py --base "$base" --data "$data" --metadata metadata
python scripts/run_benchmark.py --base "$base" --data "$data" \
  --locks /ssd/wangxinyu/spire_ecm087/logs --gpus 4 5 6 7 --pilot
# After checking pilot GFFs, GPU operation and independent scores:
python scripts/run_benchmark.py --base "$base" --data "$data" \
  --locks /ssd/wangxinyu/spire_ecm087/logs --gpus 2 3 4 5 6 7
```

The pilot consists of S. cerevisiae and N. crassa at Original and 1–2 kb; its four outputs are reused. Full inference skips successful new cases and all 180 reused cases. A failed output requires inspection before retry. Shared per-GPU locks and live process checks prevent collisions with cooperating jobs. CPU evaluation and independent validation follow successful inference.

The installer supports pre-fetched original tagged source trees when the compute server cannot clone GitHub. It installs the official pinned Python requirements and CUDA libraries, unmodified GeenuFF and Helixer, and compiles the tagged HelixerPost. Dependency versions, source commits and model hash are recorded in `provenance/`.

Main comparison uses reference-defined complete CDS. Carbon-A and Helixer also receive identical truncated-reference scoring; the Carbon-A confidence filter applies only to Carbon-A, because Helixer does not provide that score. Labels such as `partial=true` are audit metadata and never determine the reference stratum or remove a prediction.

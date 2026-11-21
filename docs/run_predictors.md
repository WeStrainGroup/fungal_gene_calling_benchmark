# Running the predictors

The release reproduces evaluation without GPU inference. These commands describe the tested scientific parameters for running new inputs. Install each tool in its own environment using its upstream instructions; the wrapper does not install software, acquire a GPU, or change a scheduler.

## Versions and models

| Tool | Tested version/source | Model |
| --- | --- | --- |
| ANNEVO | 2.3.3; commit `42c920f70f52f96c7c06c665a639ec47d9304fff` | `ANNEVO_Fungi.pt` |
| Tiberius | 2.0.7; bricks2marble 0.1.1; TensorFlow 2.20.0 | Official fungi weights |
| Carbon-A | Published pipeline/evaluator revision `ee6359990109eef6adccfd19aba8ca64b48ae514`; PyTorch 2.13.0+cu126; Transformers 4.57.6 | Carbon-A 1.2B |

Model SHA256 values and installed versions are in `metadata/software.json`. The Tiberius model URL and file hashes are in `metadata/tiberius_model_source.json`. Source fingerprints are in `metadata/source_sha256.json`. The recorded Carbon-A weights originated from its earlier GenerTeam preview and match the published model hash; its local model-loading code had Transformers-compatibility changes. The public evaluator is unmodified. Full inference is not claimed bit-identical across software revisions, kernels or hardware.

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

Use separate output directories for new runs so the released predictions remain available. For evaluation, the layout is `out/<tool>/<species>/<group>/prediction.gff` for Carbon-A and `prediction.gff3` for the other tools. Groups are `genome`, `50-100kb`, `5-10kb`, `1-2kb`.

## Fixed settings

- Carbon-A: batch size 4, SDPA attention, CPU Viterbi, minimum record length 0, native minimum CDS 96 bp, explicit sample codon table. Initial runs used 2 GPUs per case, or 4 for 1–2 kb; the six-species extension used one GPU per case. The final Rhizophagus 1–2 kb case used eight GPUs. The pipeline distributes records over devices. SDPA does not produce the authors' FlashAttention2-identical probabilities.
- ANNEVO: Fungi model, batch size 32, 4 threads, 2 workers, native decoding thresholds, no `--overlap_pred`.
- Tiberius: fungi model, batch size 8, minimum input 1,000 bp, native postprocessing. Native protein exports may use the standard code; common evaluation independently extracts CDS from the input and GFF.
- Native predictions remain unchanged. The common evaluator applies visible CDS >=200 bp and Carbon-A confidence >=0.1 or explicit NA, without filtering by predicted completeness.

Runtime tables record measured allocated GPU seconds and GNU time process RSS; this is not a hardware-normalized speed competition or aggregate peak process-tree memory measurement.

## Six-species extension

```bash
python scripts/prepare_new_inputs.py \
  --data data_tiberius --samples metadata/tiberius_samples.tsv
```

This downloads only the fixed RefSeq genomic DNA/GFF, assembly reports and checksums, excludes organellar records, and partitions the entire retained assembly. Run each tool over all four inputs for every listed species, using the layout above, then:

```bash
python scripts/evaluate_stratified.py --data data_tiberius --samples metadata/tiberius_samples.tsv \
  --cohort tiberius_set --output data_tiberius/evaluation
python scripts/validate_stratified.py --data data_tiberius --output data_tiberius/evaluation
```

All 12 unique species must finish before computing the combined mean and sample SD. Original records shorter than a target fragmentation range remain unchanged; `inputs.tsv` records their counts and base contributions.

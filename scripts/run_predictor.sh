#!/usr/bin/env bash
# Run one input with an already installed predictor and local model.
set -euo pipefail
tool=$1
input=$2
output=$3
codon_table=${4:-}
python=${PREDICTOR_PYTHON:-python}
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
export TF_NUM_INTRAOP_THREADS=4 TF_NUM_INTEROP_THREADS=2 TF_FORCE_GPU_ALLOW_GROWTH=true
mkdir -p "$output"
case "$tool" in
  carbon)
    bash "$CARBON_DIR/annotation_pipeline/run_annotation.sh" \
      --input "$input" --accession benchmark --model "$CARBON_DIR" \
      --output-dir "$output/work" --codon-table "$codon_table" \
      --gpu-count "${GPU_COUNT:-1}" --batch-size 4 --attn-implementation sdpa \
      --viterbi-backend cpu --min-record-length 0 --min-cds-length 96
    work=$output/work
    [[ -f $work/prediction.gff ]] || work=${work}_WARNING
    cp "$work"/prediction.{gff,fna,faa} "$output/"
    ;;
  annevo)
    mkdir -p "$output/work"
    "$python" "$ANNEVO_DIR/annotation.py" -g "$input" \
      -m "$ANNEVO_DIR/saved_model/ANNEVO_Fungi.pt" -l Fungi \
      -o "$output/prediction.gff3" -t 4 --batch_size 32 --num_workers 2 \
      --tmp_path "$output/work" --show_log
    ;;
  tiberius)
    "$python" -m tiberius.main --genome "$input" --model "$TIBERIUS_MODEL" \
      --out "$output/prediction.gff3" --codingseq "$output/native.cds.fa" \
      --protseq "$output/native.proteins.fa" --batch_size 8 --min_genome_seqlen 1000
    ;;
  *) echo 'Choose carbon, annevo or tiberius' >&2; exit 2 ;;
esac

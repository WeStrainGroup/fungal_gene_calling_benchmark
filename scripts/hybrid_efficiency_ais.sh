#!/usr/bin/env bash
# Keep original environment prefixes inside the AIS container.
set -euo pipefail
work=$1
index=$2
round=${3:-1}
[[ "$round" =~ ^[123]$ ]]
mkdir -p /ssd
if [[ ! -e /ssd/wangxinyu && ! -L /ssd/wangxinyu ]]; then
  ln -s "$work/software" /ssd/wangxinyu
fi
[[ $(readlink -f /ssd/wangxinyu) == "$work/software" ]]
export GNU_TIME=/ssd/wangxinyu/gnu_time
export BENCHMARK_IMAGE=192.168.108.1:5000/pytorch/pytorch2.3-admin:cuda_12.4
exec /ssd/wangxinyu/carbon_a_benchmark_20261117/runtime/bin/python \
  "$work/scripts/run_hybrid_efficiency.py" \
  --inputs "$work/inputs" --metadata "$work/software.json" \
  --cases "$work/cases.tsv" --output "$work/hybrid/round$round" --species-index "$index"

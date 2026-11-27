#!/usr/bin/env python3
"""Time one original assembly: split at 10 kb, predict serially, merge GFFs."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['inputs', 'cases', 'metadata', 'output']:
        p.add_argument('--'+key, required=True, type=Path)
    p.add_argument('--species-index', required=True, type=int)
    p.add_argument('--config', type=Path)
    a = p.parse_args()
    with a.cases.open() as handle:
        cases = [r for r in csv.DictReader(handle, delimiter='\t') if r['group'] == 'genome']
    cases.sort(key=lambda r:r['species'])
    if not 0 <= a.species_index < len(cases):
        p.error('species-index out of range')
    row = cases[a.species_index]; species = row['species']
    source = a.inputs/species/'genome.fna'; code = a.inputs/species/'codon_table.json'
    assert hashlib.sha256(source.read_bytes()).hexdigest() == row['fasta_sha256']
    from ctc import codon_table
    genetic_code = int(row['genetic_code'])
    assert codon_table(genetic_code) == json.loads(code.read_text()), 'Frozen genetic code differs'
    output = a.output/species; output.mkdir(parents=True, exist_ok=False)
    scripts = Path(__file__).resolve().parent
    subprocess.run([sys.executable, str(scripts/'check_efficiency_resources.py'),
                    '--output', str(output), '--platform', 'ais'], check=True)
    base = Path('/ssd/wangxinyu/carbon_a_benchmark_20261117')
    old = base/'three_predictor_full_genome/software'
    config = dict(environments={'carbon':str(base/'runtime'), 'tiberius':str(old/'tiberius_env')},
                  carbon_dir=str(base/'Carbon-A'), tiberius_model=str(old/'tiberius_model/fungi'),
                  gnu_time=os.environ.get('GNU_TIME', '/ssd/wangxinyu/gnu_time'))
    if a.config:
        raw = json.loads(a.config.read_text())
        def path(value):
            return str((a.config.resolve().parent/Path(value).expanduser()).resolve())
        config = {k:path(raw[k]) for k in ['carbon_dir', 'tiberius_model', 'gnu_time']}
        config['environments'] = {t:path(raw['environments'][t]) for t in ['carbon','tiberius']}
    expected = json.loads(a.metadata.read_text())
    models = {'carbon_a':Path(config['carbon_dir'])/'model.safetensors',
              'tiberius':Path(config['tiberius_model'])/'weights.h5'}
    observed = {}
    for tool, model in models.items():
        sha = hashlib.sha256()
        with model.open('rb') as handle:
            for chunk in iter(lambda:handle.read(8*1024**2), b''):
                sha.update(chunk)
        assert sha.hexdigest() == expected['models'][tool]['sha256'], tool
        observed[tool+'_model_sha256'] = sha.hexdigest()
    for tool in ['tiberius', 'carbon']:
        env = Path(config['environments'][tool])
        query = 'import importlib.metadata as m,json,sys;print(json.dumps({k:m.version(k) for k in sys.argv[1:]}))'
        packages = json.loads(subprocess.check_output([str(env/'bin/python'), '-c', query, *expected[tool]['packages']],text=True))
        assert packages == expected[tool]['packages'], tool
        observed[tool] = packages
    observed['script_sha256'] = {name:hashlib.sha256((scripts/name).read_bytes()).hexdigest()
        for name in ['run_hybrid_efficiency.py','ctc.py','check_efficiency_resources.py']}
    observed['gnu_time_sha256'] = hashlib.sha256(Path(config['gnu_time']).read_bytes()).hexdigest()
    (output/'software.json').write_text(json.dumps(observed,indent=2)+'\n')
    env = dict(os.environ, LC_ALL='C', HF_HUB_OFFLINE='1', TRANSFORMERS_OFFLINE='1',
        OPENBLAS_NUM_THREADS='4', OMP_NUM_THREADS='4', MKL_NUM_THREADS='4',
        TF_NUM_INTRAOP_THREADS='4', TF_NUM_INTEROP_THREADS='2', TF_FORCE_GPU_ALLOW_GROWTH='true')
    command = [sys.executable, str(scripts/'ctc.py'), '--input', str(source),
        '--output', str(output/'ctc'), '--min-length', '1000', '--cutoff-length', '10000', '--confidence-threshold', '0.1',
        '--tiberius-path', str(Path(config['environments']['tiberius'])/'bin/python'),
        '--tiberius-model', config['tiberius_model'], '--tiberius-args=--batch_size 8',
        '--carbon-a-path', str(Path(config['environments']['carbon'])/'bin/python'),
        '--carbon-a-model', config['carbon_dir'], '--codon-table', str(genetic_code),
        '--carbon-a-args=--accession benchmark --gpu-count 1 --batch-size 4 --attn-implementation sdpa '
        '--viterbi-backend cpu --min-record-length 0 --min-cds-length 96']
    load_start = os.getloadavg()[0]
    start = time.monotonic()
    with (output/'run.log').open('w') as log:
        proc = subprocess.run([config['gnu_time'], '-f', '%U\t%S\t%M', '-o', str(output/'time.tsv'),
            *command], env=env, stdout=log, stderr=subprocess.STDOUT)
    elapsed = time.monotonic()-start
    exit_code = proc.returncode
    details_path = output/'ctc/combined/run.json'
    details = json.loads(details_path.read_text()) if details_path.exists() else {}
    counts = details.get('routing', {})
    valid = (details.get('status') == 'complete' and
        counts.get('skipped', {}).get('records') == 0 and
        sum(x['bp'] for x in counts.values()) == int(row['bp']) and
        sum(x['records'] for x in counts.values()) == int(row['records']) and
        (output/'ctc/combined/prediction.gff3').is_file())
    if not valid:
        exit_code = exit_code or 1
    split_seconds = details.get('split_seconds', 0)
    merge_seconds = details.get('merge_seconds', 0)
    cds = details.get('cds_features', 0)
    result = dict(species=species,group='genome',tool='ctc',ctc_version=details.get('version'),confidence_threshold=details.get('confidence_threshold'),records=row['records'],bp=row['bp'],
        elapsed_seconds=elapsed,bp_per_second=int(row['bp'])/elapsed,split_seconds=split_seconds,
        merge_seconds=merge_seconds,exit_code=exit_code,passed=exit_code==0,cds_features=cds,
        job_id=os.environ.get('BENCHMARK_RUN_ID'),node=os.uname().nodename,
        platform='ais',host_load1_start=load_start,host_load1_end=os.getloadavg()[0])
    with (output/'metrics.tsv').open('w') as handle:
        w=csv.DictWriter(handle,fieldnames=list(result),delimiter='\t',lineterminator='\n');w.writeheader();w.writerow(result)
    print(json.dumps(result),flush=True)
    if exit_code:
        raise SystemExit(exit_code)

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Run one species runtime benchmark, or the eight compatibility cases."""
import argparse
import csv
import hashlib
import json
import os
import random
import shutil
from pathlib import Path
import subprocess
import sys
import time

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--inputs', type=Path, required=True)
p.add_argument('--metadata', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--cases', type=Path, required=True)
p.add_argument('--platform', choices=['slurm', 'ais'], default='slurm')
p.add_argument('--config', type=Path, help='JSON paths to installed tool environments and models')
mode = p.add_mutually_exclusive_group(required=True)
mode.add_argument('--species-index', type=int)
mode.add_argument('--preflight', action='store_true')
a = p.parse_args()
rows = list(csv.DictReader(a.cases.open(), delimiter='\t'))
species_names = sorted({r['species'] for r in rows})
tools = ['helixer', 'carbon', 'annevo', 'tiberius']
if a.preflight:
    selected = [r for r in rows if (r['species'], r['group']) in {
        ('Saccharomyces_cerevisiae', '1-2kb'), ('Rhizophagus_irregularis', 'genome')}]
else:
    if not 0 <= a.species_index < len(species_names):
        p.error('--species-index is outside the manifest')
    species = species_names[a.species_index]
    selected = [r for r in rows if r['species'] == species]
    if len(selected) != 5 or len({r['group'] for r in selected}) != 5:
        raise ValueError(f'{species}: expected five distinct input conditions')
    random.Random(20261124 + a.species_index).shuffle(selected)
    a.output = a.output/species
schedule = []
for index, row in enumerate(selected):
    shift = 0 if a.preflight else (a.species_index + index) % len(tools)
    for tool in tools[shift:] + tools[:shift]:
        schedule.append(dict(row, tool=tool, order=len(schedule)+1))
a.output.mkdir(parents=True, exist_ok=False)
with (a.output/'order.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(schedule[0]), delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(schedule)
# Input verification is outside predictor timing.
for row in selected:
    path = a.inputs/row['species']/f"{row['group']}.fna"
    if hashlib.sha256(path.read_bytes()).hexdigest() != row['fasta_sha256']:
        raise ValueError(f'{path}: input hash mismatch')
if a.platform == 'slurm' and not a.config:
    shutil.copyfile(a.metadata.parent/'software/ubuntu22.sif.sha256', a.output/'image.sha256')
scripts = Path(__file__).resolve().parent
subprocess.run([sys.executable, str(scripts/'check_efficiency_resources.py'),
                '--output', str(a.output), '--platform', a.platform], check=True)
hardware = json.loads((a.output/'hardware.json').read_text())
# The Slurm wrapper binds the staged tree at its original environment prefix.
base = Path('/ssd/wangxinyu/carbon_a_benchmark_20261117')
old = base/'three_predictor_full_genome/software'
helixer = base/'helixer_extension/software'
inputs = a.inputs
environments = {'helixer': helixer/'env', 'carbon': base/'runtime',
                'annevo': old/'annevo_env', 'tiberius': old/'tiberius_env'}
expected = json.loads(a.metadata.read_text())
models = {'carbon_a': base/'Carbon-A/model.safetensors',
          'annevo': old/'ANNEVO/saved_model/ANNEVO_Fungi.pt',
          'tiberius': old/'tiberius_model/fungi/weights.h5',
          'helixer': helixer/'models/fungi/fungi_v0.3_a_0100.h5'}
carbon_dir, annevo_dir = base/'Carbon-A', old/'ANNEVO'
tiberius_model, helixer_model = old/'tiberius_model/fungi', helixer/'models/fungi/fungi_v0.3_a_0100.h5'
if a.config:
    config = json.loads(a.config.read_text())
    def configured_path(value):
        return (a.config.resolve().parent/Path(value).expanduser()).resolve()
    environments = {tool: configured_path(config['environments'][tool]) for tool in tools}
    carbon_dir = configured_path(config['carbon_dir'])
    annevo_dir = configured_path(config['annevo_dir'])
    tiberius_model = configured_path(config['tiberius_model'])
    helixer_model = configured_path(config['helixer_model'])
    models = {'carbon_a': carbon_dir/'model.safetensors',
              'annevo': annevo_dir/'saved_model/ANNEVO_Fungi.pt',
              'tiberius': tiberius_model/'weights.h5', 'helixer': helixer_model}
    os.environ['GNU_TIME'] = str(configured_path(config['gnu_time']))
for tool, path in models.items():
    sha = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(8*1024**2), b''):
            sha.update(chunk)
    if sha.hexdigest() != expected['models'][tool]['sha256']:
        raise ValueError(f'{tool}: model hash differs from the frozen accuracy model')
observed = {}
for tool, env in environments.items():
    query = ('import importlib.metadata as m,json,sys;'
             'print(json.dumps(dict(python=sys.version,packages={k:m.version(k) for k in sys.argv[1:]})))')
    observed[tool] = json.loads(subprocess.check_output(
        [str(env/'bin/python'), '-c', query, *expected[tool]['packages']], text=True))
    if observed[tool]['packages'] != expected[tool]['packages']:
        raise ValueError(f'{tool}: packages differ from the frozen accuracy environment')
observed['model_sha256'] = {key: value['sha256'] for key, value in expected['models'].items()}
observed['script_sha256'] = {name: hashlib.sha256((scripts/name).read_bytes()).hexdigest()
                             for name in ['run_efficiency.py', 'run_predictor.sh', 'check_efficiency_resources.py']}
observed['gnu_time_sha256'] = hashlib.sha256(Path(os.environ.get(
    'GNU_TIME', '/ssd/wangxinyu/gnu_time')).read_bytes()).hexdigest()
(a.output/'software.json').write_text(json.dumps(observed, indent=2)+'\n')

fields = ['order', 'species', 'group', 'tool', 'records', 'bp', 'elapsed_seconds',
          'bp_per_second', 'user_seconds', 'system_seconds', 'max_process_rss_kib',
          'exit_code', 'cds_features', 'passed', 'job_id', 'node', 'platform',
          'host_load1_start', 'host_load1_end', 'cpu_throttled_seconds']
def throttled_seconds():
    values = []
    for filename in hardware.get('cpu_stat_paths', []):
        path = Path(filename)
        if path.is_file():
            stats = dict(line.split() for line in path.read_text().splitlines())
            if 'throttled_usec' in stats:
                values.append(int(stats['throttled_usec'])/1e6)
            elif 'throttled_time' in stats:
                values.append(int(stats['throttled_time'])/1e9)
    return max(values, default=0)


results = []
with (a.output/'metrics.tsv').open('w', buffering=1) as metrics:
    writer = csv.DictWriter(metrics, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for row in schedule:
        species, group, tool = row['species'], row['group'], row['tool']
        environment = environments[tool]
        dest = a.output/tool/group
        if a.preflight:
            dest = a.output/tool/species/group
        dest.mkdir(parents=True, exist_ok=False)
        env = dict(os.environ, GPU_COUNT='1', LC_ALL='C',
                   CARBON_DIR=str(carbon_dir), ANNEVO_DIR=str(annevo_dir),
                   TIBERIUS_MODEL=str(tiberius_model),
                   HELIXER_MODEL=str(helixer_model),
                   HF_HUB_OFFLINE='1', TRANSFORMERS_OFFLINE='1', OPENBLAS_NUM_THREADS='4',
                   PREDICTOR_PYTHON=str(environment/'bin/python'))
        env['PATH'] = str(environment/'bin')+':'+env.get('PATH', '')
        env['LD_LIBRARY_PATH'] = str(environment/'lib')+':'+env.get('LD_LIBRARY_PATH', '')
        command = [os.environ.get('GNU_TIME', '/ssd/wangxinyu/gnu_time'), '-f', '%U\t%S\t%M', '-o', str(dest/'time.tsv'),
                   'bash', str(scripts/'run_predictor.sh'), tool,
                   str(inputs/species/f'{group}.fna'), str(dest), str(inputs/species/'codon_table.json')]
        print(f'START {row["order"]}/{len(schedule)} {tool}/{species}/{group}', flush=True)
        load_start, throttle_start = os.getloadavg()[0], throttled_seconds()
        with (dest/'run.log').open('w') as log:
            start = time.monotonic()
            proc = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
            elapsed = time.monotonic()-start
        load_end, throttle_end = os.getloadavg()[0], throttled_seconds()
        # Output validation and scratch cleanup are outside the timed interval.
        gff = dest/('prediction.gff' if tool == 'carbon' else 'prediction.gff3')
        cds = 0
        if gff.is_file():
            for line in gff.open():
                parts = line.rstrip().split('\t')
                if len(parts) == 9 and parts[2] == 'CDS':
                    cds += 1
        passed = proc.returncode == 0 and gff.is_file() and gff.stat().st_size > 0
        if a.preflight:
            passed = passed and cds > 0
        usage = ['NA', 'NA', 'NA']
        if (dest/'time.tsv').is_file():
            for line in (dest/'time.tsv').read_text().splitlines():
                parts = line.split('\t')
                if len(parts) == 3:
                    usage = parts
        result = dict(order=row['order'], tool=tool, species=species, group=group,
                      records=row['records'], bp=row['bp'], elapsed_seconds=round(elapsed, 6),
                      bp_per_second=round(int(row['bp'])/elapsed, 3),
                      user_seconds=usage[0], system_seconds=usage[1], max_process_rss_kib=usage[2],
                      exit_code=proc.returncode, cds_features=cds, passed=passed,
                      job_id=os.environ.get('SLURM_JOB_ID', os.environ.get('BENCHMARK_RUN_ID')), node=os.environ.get('SLURMD_NODENAME', os.uname().nodename),
                      platform=a.platform, host_load1_start=load_start, host_load1_end=load_end,
                      cpu_throttled_seconds=round(throttle_end-throttle_start, 6))
        writer.writerow(result)
        results.append(result)
        if a.preflight:
            (a.output/'cases.json').write_text(json.dumps(results, indent=2)+'\n')
        print(json.dumps(result), flush=True)
        if passed:
            for name in ['work', 'work_WARNING']:
                if (dest/name).is_dir():
                    shutil.rmtree(dest/name)
if not all(r['passed'] for r in results):
    raise SystemExit('Some cases failed; inspect metrics.tsv and native logs')
print(f'Completed {len(results)} cases.', flush=True)

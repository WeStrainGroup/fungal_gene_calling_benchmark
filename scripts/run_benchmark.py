#!/usr/bin/env python3
"""Run the pending cases on shared GPUs, then evaluate and validate all cases."""
import argparse
import concurrent.futures
import csv
import fcntl
import json
import os
from pathlib import Path
import queue
import subprocess
import threading
import time

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--base', type=Path, required=True)
p.add_argument('--data', type=Path, required=True)
p.add_argument('--locks', type=Path, required=True)
p.add_argument('--gpus', type=int, nargs='+', default=[2, 3, 4, 5, 6, 7])
p.add_argument('--pilot', action='store_true', help='Run the four reusable Helixer compatibility cases only')
a = p.parse_args()
base, root = a.base.resolve(), a.data.resolve()
code = Path(__file__).resolve().parent
old = base/'three_predictor_full_genome/software'
assert (a.pilot or (root/'BENCHMARK_READY').exists()) and (root/'SOFTWARE_READY').exists()
if not a.pilot:
    assert json.loads((root/'provenance/helixer_pilot_validation.json').read_text())['status'] == 'PASS'
assert len(a.gpus) == len(set(a.gpus)) and a.gpus
control = (root/'logs/controller.lock').open('a')
fcntl.flock(control, fcntl.LOCK_EX | fcntl.LOCK_NB)
tasks = list(csv.DictReader((root/'tasks.tsv').open(), delimiter='\t'))
if a.pilot:
    tasks = [t for t in tasks if t['tool'] == 'helixer'
             and t['species'] in ('Saccharomyces_cerevisiae', 'Neurospora_crassa')
             and t['group'] in ('genome', '1-2kb')]
assert len(tasks) == (4 if a.pilot else 160)
mutex = threading.Lock()
active, failed = {}, []

def destination(task):
    return root/'out'/task['tool']/task['species']/task['group']

def state(status):
    with mutex:
        completed = sum((destination(t)/'COMPLETE').exists() for t in tasks)
        value = dict(status=status, pilot=a.pilot, controller_pid=os.getpid(), allowed_gpus=a.gpus,
                     completed_new=completed, expected_new=len(tasks), reused_cases=180,
                     completed_predictions=180+completed, target_cases=340,
                     active=list(active.values()), failed=failed,
                     updated=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
        target = root/('PILOT_STATUS.json' if a.pilot else 'STATUS.json')
        temporary = target.with_suffix('.tmp')
        temporary.write_text(json.dumps(value, indent=2)+'\n'); temporary.replace(target)

def run_case(task, gpu, lock):
    tool, species, group = task['tool'], task['species'], task['group']
    dest = destination(task)
    if (dest/'COMPLETE').exists():
        return
    assert not dest.exists() or not list(dest.iterdir()), f'Inspect incomplete output: {dest}'
    dest.mkdir(parents=True, exist_ok=True)
    env_path = {'carbon': base/'runtime', 'annevo': old/'annevo_env',
                'tiberius': old/'tiberius_env', 'helixer': root/'software/env'}[tool]
    env = dict(os.environ, CUDA_VISIBLE_DEVICES=str(gpu), GPU_COUNT='1',
               CARBON_DIR=str(base/'Carbon-A'), ANNEVO_DIR=str(old/'ANNEVO'),
               TIBERIUS_MODEL=str(old/'tiberius_model/fungi'),
               HELIXER_MODEL=str(root/'software/models/fungi/fungi_v0.3_a_0100.h5'),
               HF_HUB_OFFLINE='1', TRANSFORMERS_OFFLINE='1', OPENBLAS_NUM_THREADS='4',
               PREDICTOR_PYTHON=str(env_path/'bin/python'))
    env['PATH'] = str(env_path/'bin') + ':' + env.get('PATH', '')
    env['LD_LIBRARY_PATH'] = str(env_path/'lib') + ':' + env.get('LD_LIBRARY_PATH', '')
    command = ['/usr/bin/time', '-v', 'bash', str(code/'run_predictor.sh'), tool,
               str(root/'inputs'/species/f'{group}.fna'), str(dest),
               str(root/'inputs'/species/'codon_table.json')]
    key = f'{tool}/{species}/{group}'
    with (root/'logs'/f'{tool}_{species}_{group}.log').open('w') as log:
        proc = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT,
                                pass_fds=(lock.fileno(),))
        info = dict(task, gpu=gpu, pid=proc.pid, started=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
        with mutex:
            active[key] = info
        (dest/'resources.json').write_text(json.dumps(dict(info, command=command), indent=2)+'\n')
        state('INFERENCE_RUNNING')
        returncode = proc.wait()
    names = ['prediction.gff', 'prediction.fna', 'prediction.faa'] if tool == 'carbon' else ['prediction.gff3']
    success = returncode == 0 and all((dest/n).is_file() and (dest/n).stat().st_size for n in names)
    with mutex:
        del active[key]
        if success:
            (dest/'COMPLETE').write_text(time.strftime('%Y-%m-%dT%H:%M:%S%z')+'\n')
        else:
            failed.append(dict(info, exit_code=returncode))
    state('INFERENCE_RUNNING')

# Start long Carbon-A cases early while other workers clear the faster tools.
carbon, fast = queue.Queue(), queue.Queue()
for task in sorted(tasks, key=lambda t: (t['group'] != '1-2kb',
                                       -(root/'inputs'/t['species']/'genome.fna').stat().st_size)):
    if not (destination(task)/'COMPLETE').exists():
        (carbon if task['tool'] == 'carbon' else fast).put(task)

def worker(index, gpu):
    with (a.locks/f'gpu{gpu}.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        busy = subprocess.check_output(['nvidia-smi', f'--id={gpu}', '--query-compute-apps=pid',
                                        '--format=csv,noheader,nounits'], text=True).strip()
        assert not busy, f'GPU {gpu} is occupied'
        preferred, fallback = (carbon, fast) if index < len(a.gpus)//2 else (fast, carbon)
        while True:
            try:
                task = preferred.get_nowait()
            except queue.Empty:
                try:
                    task = fallback.get_nowait()
                except queue.Empty:
                    return
            run_case(task, gpu, lock)

try:
    state('STARTING')
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(a.gpus)) as pool:
        futures = [pool.submit(worker, i, gpu) for i, gpu in enumerate(a.gpus)]
        for future in futures:
            future.result()
    assert not failed, 'Failed cases need inspection'
    if not a.pilot:
        state('EVALUATING')
        python = str(base/'runtime/bin/python')
        evaluation = root/'evaluation'
        results = root/'results'
        commands = [
            [python, str(code/'evaluate_stratified.py'), '--data', str(root), '--samples', str(root/'samples.tsv'),
             '--cohort', 'fungal_benchmark', '--output', str(evaluation)],
            [python, str(code/'validate_stratified.py'), '--data', str(root), '--output', str(evaluation)],
            [python, str(code/'summarize_all_species.py'), '--inputs', str(evaluation),
             '--species', str(root/'samples.tsv'), '--output', str(results)]]
        with (root/'logs/evaluation.log').open('w') as log:
            for command in commands:
                subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
            rscript = '/home/wangxinyu/miniconda3/envs/vscode-r/bin/Rscript'
            env = dict(os.environ, R_LIBS_USER=str(old/'R_library'))
            for script, extra in [('plot_all_species.R', [str(results/'figures')]), ('validate_summary.R', [])]:
                subprocess.run([rscript, str(code/script), str(results), *extra], env=env,
                               check=True, stdout=log, stderr=subprocess.STDOUT)
    state('COMPLETE')
except Exception as error:
    with mutex:
        failed.append(dict(controller_error=str(error)))
    state('NEEDS_ATTENTION')
    raise

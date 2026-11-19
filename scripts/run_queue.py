#!/usr/bin/env python3
"""Run pending cases on shared idle GPUs, retaining existing inference processes."""
import argparse
import csv
import fcntl
import json
import os
from pathlib import Path
import re
import subprocess
import time

COHORTS = ('fungi_test_six', 'carbon_fungi_set')
TOOLS = ('tiberius', 'annevo', 'carbon')
GROUPS = ('1-2kb', '5-10kb', '50-100kb', 'genome')


def write_json(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def process_identity(pid):
    path = Path('/proc') / str(pid) / 'stat'
    if not path.exists():
        return None
    fields = path.read_text().rsplit(')', 1)[1].split()
    return None if fields[0] == 'Z' else fields[19]


def adopted_exit(item):
    if process_identity(item['pid']) == item['start_ticks']:
        return None
    with Path(item['log']).open('rb') as handle:
        handle.seek(max(0, os.fstat(handle.fileno()).st_size - 4096))
        matches = re.findall(rb'Exit status:\s*(\d+)', handle.read())
    return int(matches[-1]) if matches else -1


def main(args):
    base = args.base.resolve()
    old = base / 'three_predictor_full_genome'
    operational = base / 'public_release'
    locks = []
    for cohort in COHORTS:
        lock = (base / cohort / 'logs/queue.lock').open('a')
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        locks.append(lock)
        (base / cohort / 'logs/queue.pid').write_text(str(os.getpid()) + '\n')
    all_tasks = []
    for cohort in COHORTS:
        samples = list(csv.DictReader((base/cohort/'samples.tsv').open(), delimiter='\t'))
        all_tasks += [(cohort, tool, row['species'], group)
                      for group in GROUPS for row in samples for tool in TOOLS]

    def output(task):
        cohort, tool, species, group = task
        return base / cohort / 'out' / tool / species / group

    active = {}
    if args.adopt:
        for item in json.loads(args.adopt.read_text())['active']:
            assert item['cohort'] in COHORTS and tuple([item['cohort']] + item['task']) in all_tasks
            active[item['gpu']] = item
    claimed = {tuple([v['cohort']] + v['task']) for v in active.values()}
    tasks = [t for t in all_tasks if not (output(t)/'COMPLETE').exists() and t not in claimed]
    failed = []

    def save_state():
        for cohort in COHORTS:
            running = [{k: v[k] for k in ('task', 'gpu', 'pid')} for v in active.values() if v['cohort'] == cohort]
            errors = [v for v in failed if v['cohort'] == cohort]
            pending = sum(t[0] == cohort for t in tasks)
            state = 'RUNNING' if running or pending else 'NEEDS_ATTENTION' if errors else 'INFERENCE_COMPLETE'
            write_json(base/cohort/'STATUS.json', dict(status=state, updated=time.strftime('%Y-%m-%dT%H:%M:%S%z'),
                pending=pending, active=running, failed=errors,
                completed=sum((output(t)/'COMPLETE').exists() for t in all_tasks if t[0] == cohort),
                controller_pid=os.getpid(), allowed_gpus=args.gpus, max_active=len(args.gpus)))
        write_json(operational/'gpu_queue.json', dict(pid=os.getpid(), allowed_gpus=args.gpus,
            active=[{k: v[k] for k in ('cohort', 'task', 'gpu', 'pid', 'start_ticks', 'log')} for v in active.values()],
            pending=len(tasks), failed=failed, updated=time.strftime('%Y-%m-%dT%H:%M:%S%z')))

    while tasks or active:
        for gpu, item in list(active.items()):
            code = item['process'].poll() if 'process' in item else adopted_exit(item)
            if code is None:
                continue
            task = tuple([item['cohort']] + item['task'])
            gff = output(task) / ('prediction.gff' if task[1] == 'carbon' else 'prediction.gff3')
            if code == 0 and gff.exists() and gff.stat().st_size:
                (output(task)/'COMPLETE').write_text(time.strftime('%Y-%m-%dT%H:%M:%S%z')+'\n')
            else:
                failed.append(dict(cohort=item['cohort'], task=item['task'], exit_code=code))
            for key in ('handle', 'lock'):
                if key in item:
                    item[key].close()
            del active[gpu]
        rows = subprocess.check_output(['nvidia-smi', '--query-gpu=index,memory.used,utilization.gpu',
                                        '--format=csv,noheader,nounits'], text=True)
        for row in rows.splitlines():
            if not tasks:
                break
            gpu, memory, utilization = map(int, row.split(','))
            if gpu not in args.gpus or gpu in active or memory > 50 or utilization:
                continue
            lock = (args.locks/f'gpu{gpu}.lock').open('a')
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                lock.close()
                continue
            occupied = subprocess.check_output(['nvidia-smi', f'--id={gpu}', '--query-compute-apps=pid',
                                                '--format=csv,noheader,nounits'], text=True).strip()
            if occupied:
                lock.close()
                continue
            task = tasks.pop(0)
            if (output(task)/'COMPLETE').exists():
                lock.close()
                continue
            cohort, tool, species, group = task
            root = base/cohort
            out = output(task)
            out.mkdir(parents=True, exist_ok=True)
            runtime = base/'runtime' if tool == 'carbon' else old/'software'/('annevo_env' if tool == 'annevo' else 'tiberius_env')
            env = os.environ.copy()
            env.update(CUDA_VISIBLE_DEVICES=str(gpu), GPU_COUNT='1', CARBON_DIR=str(base/'Carbon-A'),
                ANNEVO_DIR=str(old/'software/ANNEVO'), TIBERIUS_MODEL=str(old/'software/tiberius_model/fungi'),
                HF_HUB_OFFLINE='1', TRANSFORMERS_OFFLINE='1', OPENBLAS_NUM_THREADS='4',
                PREDICTOR_PYTHON=str(runtime/'bin/python'), PATH=str(runtime/'bin')+':'+env.get('PATH',''),
                LD_LIBRARY_PATH=str(runtime/'lib')+':'+env.get('LD_LIBRARY_PATH',''))
            log = root/'logs'/f'{tool}_{species}_{group}.log'
            handle = log.open('w')
            command = ['/usr/bin/time', '-v', 'bash', str(root/'scripts/run_predictor.sh'), tool,
                str(root/'inputs'/species/f'{group}.fna'), str(out), str(root/'inputs'/species/'codon_table.json')]
            process = subprocess.Popen(command, env=env, stdout=handle, stderr=subprocess.STDOUT, pass_fds=(lock.fileno(),))
            active[gpu] = dict(cohort=cohort, task=list(task[1:]), gpu=gpu, pid=process.pid,
                start_ticks=process_identity(process.pid), log=str(log), process=process, handle=handle, lock=lock)
            write_json(out/'resources.json', dict(allocated_gpus=1, gpu=gpu, threads=4, pid=process.pid))
        save_state()
        if tasks or active:
            time.sleep(20)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, required=True)
    parser.add_argument('--locks', type=Path, required=True)
    parser.add_argument('--gpus', type=int, nargs='+', default=list(range(8)))
    parser.add_argument('--adopt', type=Path)
    main(parser.parse_args())

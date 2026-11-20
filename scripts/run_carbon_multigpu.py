#!/usr/bin/env python3
"""Replace the final single-GPU case with native eight-GPU Carbon-A inference."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time


def identity(pid):
    try:
        fields = (Path('/proc')/str(pid)/'stat').read_text().rsplit(')', 1)[1].split()
    except FileNotFoundError:
        return None
    return None if fields[0] == 'Z' else fields[19]


def write_json(path, data):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data, indent=2)+'\n')
    temporary.replace(path)


def main(args):
    base = args.base.resolve()
    request = json.loads(args.handoff.read_text())
    cohort, species, group = request['cohort'], request['species'], request['group']
    root = base/cohort
    out = root/'out/carbon'/species/group
    fasta = root/'inputs'/species/f'{group}.fna'
    assert not (out/'COMPLETE').exists()
    for key, item in request['processes'].items():
        assert identity(item['pid']) == item['start_ticks'], key
    queue = json.loads((base/'public_release/gpu_queue.json').read_text())
    assert queue['pending'] == 0 and len(queue['active']) == 1
    assert queue['active'][0]['pid'] == request['processes']['wrapper']['pid']
    assert not list((out/'work').rglob('*')), 'Existing inference files need inspection'
    locks = []
    for gpu in [0, 1, 3, 4, 5, 6, 7]:
        lock = (args.locks/f'gpu{gpu}.lock').open('a')
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        occupied = subprocess.check_output(['nvidia-smi', f'--id={gpu}', '--query-compute-apps=pid',
                                            '--format=csv,noheader,nounits'], text=True).strip()
        assert not occupied, gpu
        locks.append(lock)
    receipt = dict(request, supervisor_pid=os.getpid(), status='REPLACING_SINGLE_GPU',
                   input_sha256=hashlib.sha256(fasta.read_bytes()).hexdigest(),
                   previous_resources=json.loads((out/'resources.json').read_text()))
    write_json(out/'multigpu_handoff.json', receipt)
    # Only the verified scheduler and this case's inference process are signalled.
    os.kill(request['processes']['controller']['pid'], signal.SIGTERM)
    os.kill(request['processes']['inference']['pid'], signal.SIGINT)
    deadline = time.monotonic()+45
    while time.monotonic() < deadline:
        if all(identity(p['pid']) is None for p in request['processes'].values()):
            break
        time.sleep(1)
    remaining = [p for p in request['processes'].values() if identity(p['pid']) == p['start_ticks']]
    for item in remaining:
        os.kill(item['pid'], signal.SIGTERM)
    if remaining:
        time.sleep(5)
    assert all(identity(p['pid']) is None for p in request['processes'].values())
    for path in [args.locks/'gpu2.lock', base/'fungi_test_six/logs/queue.lock', root/'logs/queue.lock']:
        lock = path.open('a')
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        locks.append(lock)
    assert not subprocess.check_output(['nvidia-smi', '--query-compute-apps=pid',
                                       '--format=csv,noheader,nounits'], text=True).strip()
    runtime = base/'runtime'
    env = os.environ.copy()
    env.update(CUDA_VISIBLE_DEVICES='0,1,2,3,4,5,6,7', GPU_COUNT='8', CARBON_DIR=str(base/'Carbon-A'),
        HF_HUB_OFFLINE='1', TRANSFORMERS_OFFLINE='1', OPENBLAS_NUM_THREADS='4',
        PREDICTOR_PYTHON=str(runtime/'bin/python'), PATH=str(runtime/'bin')+':'+env.get('PATH',''),
        LD_LIBRARY_PATH=str(runtime/'lib')+':'+env.get('LD_LIBRARY_PATH',''))
    log = root/'logs'/f'carbon_{species}_{group}.gpu8.log'
    command = ['/usr/bin/time', '-v', 'bash', str(root/'scripts/run_predictor.sh'), 'carbon',
               str(fasta), str(out), str(root/'inputs'/species/'codon_table.json')]
    with log.open('w') as handle:
        process = subprocess.Popen(command, env=env, stdout=handle, stderr=subprocess.STDOUT,
                                   pass_fds=tuple(lock.fileno() for lock in locks))
        active = dict(cohort=cohort, task=['carbon', species, group], gpus=list(range(8)),
                      pid=process.pid, start_ticks=identity(process.pid), log=str(log))
        receipt.update(status='RUNNING_NATIVE_GPU8', replacement=active)
        write_json(out/'multigpu_handoff.json', receipt)
        write_json(out/'resources.json', dict(allocated_gpus=8, gpus=list(range(8)), threads_per_worker=4,
            gpu_count=8, batch_size=4, dispatch_chunk_size=4, attn_implementation='sdpa', pid=process.pid))
        (base/'public_release/gpu_queue.pid').write_text(str(os.getpid())+'\n')
        while True:
            code = process.poll()
            state = dict(pid=os.getpid(), allowed_gpus=list(range(8)), active=[active] if code is None else [],
                         pending=0, failed=[] if code in (None, 0) else [dict(task=active['task'], exit_code=code)],
                         updated=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
            write_json(base/'public_release/gpu_queue.json', state)
            status = json.loads((root/'STATUS.json').read_text())
            status.update(status='RUNNING' if code is None else 'POST_RUN_CHECK',
                          active=state['active'], pending=0, failed=state['failed'],
                          controller_pid=os.getpid(), updated=state['updated'])
            write_json(root/'STATUS.json', status)
            if code is not None:
                break
            time.sleep(20)
    assert code == 0, f'Native pipeline exit code {code}; inspect {log}'
    for name in ['prediction.gff', 'prediction.fna', 'prediction.faa']:
        assert (out/name).is_file() and (out/name).stat().st_size, name
    assert hashlib.sha256(fasta.read_bytes()).hexdigest() == receipt['input_sha256']
    (out/'COMPLETE').write_text(time.strftime('%Y-%m-%dT%H:%M:%S%z')+'\n')
    status = json.loads((root/'STATUS.json').read_text())
    status.update(status='INFERENCE_COMPLETE', completed=60, pending=0, active=[], failed=[],
                  controller_pid=os.getpid(), updated=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
    write_json(root/'STATUS.json', status)
    receipt['status'] = 'INFERENCE_COMPLETE'
    write_json(out/'multigpu_handoff.json', receipt)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, required=True)
    parser.add_argument('--locks', type=Path, required=True)
    parser.add_argument('--handoff', type=Path, required=True)
    main(parser.parse_args())

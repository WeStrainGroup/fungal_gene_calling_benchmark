#!/usr/bin/env python3
"""Run the 2–5 kb supplement on GPU 4–7, then evaluate it on CPU."""
import argparse, concurrent.futures, csv, fcntl, hashlib, json, os, queue, subprocess, threading, time
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--base',type=Path,required=True)
p.add_argument('--data',type=Path,required=True)
p.add_argument('--locks',type=Path,required=True)
a=p.parse_args()
base=a.base.resolve();root=a.data.resolve();code=Path(__file__).resolve().parent
old=base/'three_predictor_full_genome';gpus=[4,5,6,7]
(root/'logs').mkdir(exist_ok=True)
control=(root/'logs/controller.lock').open('a');fcntl.flock(control,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert (root/'INPUTS_READY').exists()
rows=list(csv.DictReader((root/'inputs.tsv').open(),delimiter='\t'))
species=[r['species'] for r in sorted(rows,key=lambda r:int(r['bp']),reverse=True)]
assert len(species)==12
active={};failed=[];mutex=threading.Lock();gpu_locks=[]
all_tasks=[(t,s) for t in ('annevo','tiberius','carbon') for s in species]
def out(tool,sp):return root/'out'/tool/sp/'2-5kb'
def state(status):
 with mutex:
  complete=sum((out(t,s)/'COMPLETE').exists() for t,s in all_tasks)
  value=dict(status=status,controller_pid=os.getpid(),allowed_gpus=gpus,completed=complete,expected=36,
   pending=36-complete-len(active)-len(failed),active=list(active.values()),failed=failed,updated=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
  tmp=root/'STATUS.tmp';tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(root/'STATUS.json')
def run_model(tool,sp,devices):
 dest=out(tool,sp)
 if (dest/'COMPLETE').exists():return
 if dest.exists() and list(dest.iterdir()):raise RuntimeError(f'Inspect incomplete output before retry: {dest}')
 dest.mkdir(parents=True,exist_ok=True)
 runtime=base/'runtime' if tool=='carbon' else old/'software'/('annevo_env' if tool=='annevo' else 'tiberius_env')
 env=os.environ.copy();env.update(CUDA_VISIBLE_DEVICES=','.join(map(str,devices)),GPU_COUNT=str(len(devices)),
  CARBON_DIR=str(base/'Carbon-A'),ANNEVO_DIR=str(old/'software/ANNEVO'),TIBERIUS_MODEL=str(old/'software/tiberius_model/fungi'),
  HF_HUB_OFFLINE='1',TRANSFORMERS_OFFLINE='1',OPENBLAS_NUM_THREADS='4',PREDICTOR_PYTHON=str(runtime/'bin/python'),
  PATH=str(runtime/'bin')+':'+env.get('PATH',''),LD_LIBRARY_PATH=str(runtime/'lib')+':'+env.get('LD_LIBRARY_PATH',''))
 command=['/usr/bin/time','-v','bash',str(code/'run_predictor.sh'),tool,str(root/'inputs'/sp/'2-5kb.fna'),str(dest),str(root/'inputs'/sp/'codon_table.json')]
 with (root/'logs'/f'{tool}_{sp}.log').open('w') as log:
  proc=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT,pass_fds=tuple(h.fileno() for h in gpu_locks))
  key=f'{tool}/{sp}';info=dict(tool=tool,species=sp,gpus=devices,pid=proc.pid)
  with mutex:active[key]=info
  (dest/'resources.json').write_text(json.dumps(dict(info,command=command,started=time.strftime('%Y-%m-%dT%H:%M:%S%z')),indent=2)+'\n')
  state('INFERENCE_RUNNING');returncode=proc.wait()
 names=['prediction.gff','prediction.fna','prediction.faa'] if tool=='carbon' else ['prediction.gff3']
 success=returncode==0 and all((dest/n).is_file() and (dest/n).stat().st_size for n in names)
 with mutex:
  del active[key]
  if success:(dest/'COMPLETE').write_text(time.strftime('%Y-%m-%dT%H:%M:%S%z')+'\n')
  else:failed.append(dict(info,exit_code=returncode))
 state('INFERENCE_RUNNING')

def single_gpu_worker(gpu):
 while True:
  try:task=pending.get_nowait()
  except queue.Empty:return
  try:run_model(*task,[gpu])
  finally:pending.task_done()

try:
 # Shared locks and process checks prevent collisions with other projects.
 for gpu in gpus:
  h=(a.locks/f'gpu{gpu}.lock').open('a');fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB);gpu_locks.append(h)
  busy=subprocess.check_output(['nvidia-smi',f'--id={gpu}','--query-compute-apps=pid','--format=csv,noheader,nounits'],text=True).strip()
  assert not busy,f'GPU {gpu} is occupied'
 pending=queue.Queue()
 for sp in species:
  for tool in ('tiberius','annevo'):pending.put((tool,sp))
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  futures=[pool.submit(single_gpu_worker,gpu) for gpu in gpus]
  for future in futures:future.result()
 # Native multi-GPU scheduling avoids leaving one large Carbon-A genome on one GPU.
 for sp in species:run_model('carbon',sp,gpus)
 for h in gpu_locks:h.close()
 gpu_locks.clear()
 if failed:raise RuntimeError('Failed inference cases require inspection')
 state('EVALUATING')
 python=str(base/'runtime/bin/python');evaluation=root/'evaluation_v0.2';combined=root/'combined_results'
 commands=[
  [python,str(code/'evaluate_stratified.py'),'--data',str(root),'--samples',str(root/'samples.tsv'),'--cohort','supplement_2_5kb','--groups','2-5kb','--output',str(evaluation)],
  [python,str(code/'validate_stratified.py'),'--data',str(root),'--output',str(evaluation)],
  [python,str(code/'summarize_all_species.py'),'--inputs',*[str(base/c/'evaluation_v0.2') for c in ['three_predictor_full_genome','fungi_test_six','carbon_fungi_set']],str(evaluation),'--groups','genome','50-100kb','5-10kb','2-5kb','1-2kb','--output',str(combined)]]
 with (root/'logs/evaluation.log').open('w') as log:
  for command in commands:subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT)
  env=os.environ.copy();env['R_LIBS_USER']=str(old/'software/R_library')
  subprocess.run(['/home/wangxinyu/miniconda3/envs/vscode-r/bin/Rscript',str(code/'plot_all_species.R'),str(combined),str(combined/'figures')],env=env,check=True,stdout=log,stderr=subprocess.STDOUT)
 state('COMPLETE')
except Exception as exc:
 with mutex:failed.append(dict(controller_error=str(exc)))
 state('NEEDS_ATTENTION')
 raise
finally:
 for h in gpu_locks:h.close()

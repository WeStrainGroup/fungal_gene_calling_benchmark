#!/usr/bin/env python3
"""Verify GPU, CPU allocation and host-memory limits on Slurm or AIS."""
import argparse
import json
import os
from pathlib import Path
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--platform', choices=['slurm', 'ais'], default='slurm')
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=True)
proc_cgroup = Path('/proc/self/cgroup').read_text()
limits = {}
swap_limits = {}
for line in proc_cgroup.splitlines():
    _, controllers, relative = line.split(':', 2)
    if not controllers:
        root, filename = Path('/sys/fs/cgroup'), 'memory.max'
    elif 'memory' in controllers.split(','):
        root, filename = Path('/sys/fs/cgroup/memory'), 'memory.limit_in_bytes'
    else:
        continue
    current = root/relative.lstrip('/')
    for directory in [current, *current.parents]:
        if directory != root and root not in directory.parents:
            break
        path = directory/filename
        if path.is_file():
            limits[str(path)] = path.read_text().strip()
        swap = directory/('memory.swap.max' if filename == 'memory.max' else 'memory.memsw.limit_in_bytes')
        if swap.is_file():
            swap_limits[str(swap)] = swap.read_text().strip()
cpu_limits = {}
cpu_stats = []
for line in proc_cgroup.splitlines():
    _, controllers, relative = line.split(':', 2)
    if not controllers:
        root, unified = Path('/sys/fs/cgroup'), True
    elif 'cpu' in controllers.split(','):
        root, unified = Path('/sys/fs/cgroup/cpu'), False
    else:
        continue
    current = root/relative.lstrip('/')
    for directory in [current, *current.parents]:
        if directory != root and root not in directory.parents:
            break
        if unified and (directory/'cpu.max').is_file():
            quota, period = (directory/'cpu.max').read_text().split()
            if quota != 'max':
                cpu_limits[str(directory)] = int(quota)/int(period)
        elif not unified and (directory/'cpu.cfs_quota_us').is_file():
            quota = int((directory/'cpu.cfs_quota_us').read_text())
            period = int((directory/'cpu.cfs_period_us').read_text())
            if quota > 0:
                cpu_limits[str(directory)] = quota/period
        if (directory/'cpu.stat').is_file():
            cpu_stats.append(str(directory/'cpu.stat'))
cpu_quota = min(cpu_limits.values(), default=None)
numeric_limits = [int(v) for v in limits.values() if v != 'max']
memory_limit = min(numeric_limits, default=0)
gpu = subprocess.check_output(['nvidia-smi', '--query-gpu=name,uuid,memory.total,driver_version',
                               '--format=csv,noheader,nounits'], text=True).strip().splitlines()
cpu = subprocess.check_output(['lscpu'], text=True)
affinity = sorted(os.sched_getaffinity(0))
topology = subprocess.check_output(['lscpu', '-p=CPU,CORE,SOCKET'], text=True)
physical = {tuple(line.split(',')[1:]) for line in topology.splitlines()
            if not line.startswith('#') and int(line.split(',')[0]) in affinity}
import torch
device = torch.cuda.get_device_name(0)
bf16 = torch.cuda.is_bf16_supported()
checks = {'one_visible_cuda_gpu': torch.cuda.device_count() == 1,
          'a100_40gb': 'A100' in device and '40GB' in device,
          'bf16': bf16,
          'memory_32gib': memory_limit == 32*1024**3}
checks['four_physical_cores' if a.platform == 'slurm' else 'four_cpu_quota'] = (
    len(physical) == len(affinity) == 4 if a.platform == 'slurm' else cpu_quota == 4)
receipt = dict(platform=a.platform, cpu_quota=cpu_quota, cpu_limits=cpu_limits,
               cpu_stat_paths=cpu_stats, physical_cores_exclusive=a.platform == 'slurm',
               image=os.environ.get('BENCHMARK_IMAGE'), checks=checks, passed=all(checks.values()), gpu=gpu,
               cpu_affinity=affinity, physical_cores=sorted(physical),
               cgroup_limits=limits, proc_cgroup=proc_cgroup,
               swap_limits=swap_limits,
               host_swap=[line for line in Path('/proc/meminfo').read_text().splitlines()
                          if line.startswith(('SwapTotal:', 'SwapFree:'))],
               slurm_job_id=os.environ.get('SLURM_JOB_ID'),
               node=os.environ.get('SLURMD_NODENAME'), container_hostname=os.uname().nodename,
               cuda_visible_devices=os.environ.get('CUDA_VISIBLE_DEVICES'),
               torch_version=torch.__version__, cuda_version=torch.version.cuda,
               cpu=cpu)
(a.output/'hardware.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps(receipt, indent=2), flush=True)
if not all(checks.values()):
    raise SystemExit('Resource preflight failed; inspect hardware.json before timing')
x = torch.randn((1024, 1024), device='cuda', dtype=torch.bfloat16)
y = x @ x
torch.cuda.synchronize()
assert torch.isfinite(y).all().item()
print('BF16 CUDA matrix multiplication passed', flush=True)

#!/usr/bin/env python3
"""Evaluate completed native outputs, then wait for the rest without running inference."""
import argparse
import csv
import fcntl
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

TOOLS = ('annevo', 'tiberius', 'carbon')
GROUPS = ('genome', '50-100kb', '5-10kb', '1-2kb')


def main(args):
    scripts = Path(__file__).resolve().parent
    data, output = args.data.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    species = [r['species'] for r in csv.DictReader(args.samples.open(), delimiter='\t')]
    cases = [(s, g, t) for s in species for g in GROUPS for t in TOOLS]
    with (output/'followup.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        previous, first = None, True
        while True:
            pending = [dict(species=s, group=g, tool=t) for s, g, t in cases
                       if not (data/'out'/t/s/g/'COMPLETE').is_file()]
            if len(pending) != previous:
                print(f'{len(cases)-len(pending)}/{len(cases)} native cases complete', flush=True)
                previous = len(pending)
            receipt = dict(status='WAITING_NATIVE' if pending else 'EVALUATING',
                           updated_utc=datetime.now(timezone.utc).isoformat(), pending=pending)
            (output/'followup.json').write_text(json.dumps(receipt, indent=2)+'\n')
            if first or not pending:
                subprocess.run([sys.executable, str(scripts/'evaluate_stratified.py'),
                    '--data', str(data), '--output', str(output), '--cohort', args.cohort,
                    '--samples', str(args.samples.resolve()), '--completed-only'], check=True)
                subprocess.run([sys.executable, str(scripts/'validate_stratified.py'),
                    '--data', str(data), '--output', str(output)], check=True)
                first = False
                if json.loads((output/'STATUS.json').read_text())['status'] == 'COMPLETE':
                    subprocess.run([args.rscript, str(scripts/'plot_stratified.R'),
                                    str(output), str(output/'figures')], check=True)
                    receipt.update(status='VALIDATED', pending=[],
                                   updated_utc=datetime.now(timezone.utc).isoformat())
                    (output/'followup.json').write_text(json.dumps(receipt, indent=2)+'\n')
                    return
            time.sleep(60)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--samples', type=Path, required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--rscript', default='Rscript')
    main(p.parse_args())

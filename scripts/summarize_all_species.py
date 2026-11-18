#!/usr/bin/env python3
"""Combine the frozen species once and report species-equal mean and sample SD."""
import argparse
import collections
import csv
import fcntl
import hashlib
import json
import math
import statistics
import subprocess
import time
from pathlib import Path

GROUPS = ('genome', '50-100kb', '5-10kb', '1-2kb')
TOOLS = ('annevo', 'tiberius', 'carbon')
SCOPES = ('complete', 'partial', 'missing_5prime', 'missing_3prime', 'missing_both')
LEVELS = ('nucleotide', 'exon', 'gene')
KEY = ('species', 'group', 'tool', 'stratum', 'level')


def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_table(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else [], delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def unique_rows(rows, fields):
    """Shared species must be scientifically identical, not averaged twice."""
    unique = {}
    for original in rows:
        row = {k: v for k, v in original.items() if k != 'cohort'}
        key = tuple(row[k] for k in fields)
        if key in unique and unique[key] != row:
            raise ValueError(f'Conflicting duplicate: {key}')
        unique[key] = row
    return [unique[k] for k in sorted(unique)]


def mean_sd(rows, species):
    groups = collections.defaultdict(list)
    for row in rows:
        groups[tuple(row[k] for k in KEY[1:])].append(row)
    result = []
    for key, items in sorted(groups.items()):
        assert len(items) == len(species) and {r['species'] for r in items} == set(species)
        for metric in ('precision', 'recall', 'f1'):
            values = [float(r[metric]) for r in items if r[metric] not in ('', 'NA')]
            assert all(math.isfinite(v) for v in values)
            n = len(values)
            result.append(dict(zip(KEY[1:], key), metric=metric, n_species=n,
                total_species=len(species), mean=statistics.mean(values) if n else 'NA',
                sd=statistics.stdev(values) if n > 1 else 'NA'))
    return result


def aggregate(roots, manifest, output):
    species = [r['species'] for r in read_table(manifest)]
    assert len(species) == len(set(species)) == 12
    output.mkdir(parents=True, exist_ok=True)
    collected = {name: [] for name in ('metrics', 'diagnostics', 'inputs')}
    snapshots = []
    for root in roots:
        status = json.loads((root/'STATUS.json').read_text())
        check = json.loads((root/'independent_validation.json').read_text())
        assert status['version'] == '0.2.0'
        assert (status['min_record_bp'], status['min_cds_bp'], status['confidence_threshold']) == (1000, 200, 0.1)
        assert check['status'] == 'PASS' and check['cases'] == status['completed_cases']
        assert status['source_files_unchanged'] and check['source_files_unchanged']
        hashes = {}
        for name in collected:
            path = root/f'{name}.tsv'
            hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
            collected[name].extend(read_table(path))
        snapshots.append(dict(cohort=status['cohort'], status=status['status'],
            completed_cases=status['completed_cases'], expected_cases=status['expected_cases'], sha256=hashes))
    metrics = unique_rows(collected['metrics'], KEY)
    diagnostics = unique_rows(collected['diagnostics'], ('species', 'group', 'tool', 'category'))
    inputs = unique_rows(collected['inputs'], ('species', 'group'))
    expected = {(s, g, t, scope, level) for s in species for g in GROUPS for t in TOOLS
                for scope in (SCOPES if t == 'carbon' else ('complete',)) for level in LEVELS}
    observed = {tuple(r[k] for k in KEY) for r in metrics}
    assert observed <= expected, 'Unexpected species, condition or metric'
    complete = all(s['status'] == 'COMPLETE' for s in snapshots)
    if complete:
        assert observed == expected, 'Missing final metric rows'
    else:
        assert not (output/'mean_sd.tsv').exists(), 'Do not replace a final summary with partial inputs'
    for name, rows in [('metrics', metrics), ('diagnostics', diagnostics), ('inputs', inputs)]:
        write_table(output/f'{name}.tsv', rows)
    if complete:
        write_table(output/'mean_sd.tsv', mean_sd(metrics, species))
    status = dict(aggregation_version='0.2.1', scoring_version='0.2.0',
        status='COMPLETE' if complete else 'PARTIAL', expected_species=12,
        completed_cases=len({tuple(r[k] for k in KEY[:3]) for r in metrics}), expected_cases=144,
        metric_rows=len(metrics), expected_metric_rows=len(expected),
        duplicate_metric_rows_removed=len(collected['metrics'])-len(metrics),
        sample_sd_ddof=1, missing_values='excluded, never replaced with zero',
        sources=snapshots)
    (output/'STATUS.json').write_text(json.dumps(status, indent=2)+'\n')
    return status


def run(args):
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    with (output/'summary.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        status = aggregate(args.inputs, args.species, output)
        print(json.dumps(status, indent=2), flush=True)
        while args.wait and status['status'] != 'COMPLETE':
            time.sleep(60)
            # Wait for all independent validators, not just inference completion.
            ready = True
            for root in args.inputs:
                current = json.loads((root/'STATUS.json').read_text())
                check = json.loads((root/'independent_validation.json').read_text())
                ready &= (current['status'] == 'COMPLETE' and check['status'] == 'PASS'
                          and check['cases'] == current['expected_cases'])
            if ready:
                status = aggregate(args.inputs, args.species, output)
                print(json.dumps(status, indent=2), flush=True)
        if status['status'] == 'COMPLETE' and args.rscript:
            subprocess.run([args.rscript, str(Path(__file__).with_name('plot_all_species.R')),
                            str(output), str(output/'figures')], check=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, nargs=3, required=True,
                        help='Validated ANNEVO, Tiberius and Carbon-A cohort result directories')
    parser.add_argument('--species', type=Path, default=Path(__file__).resolve().parents[1]/'metadata/species.tsv')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--wait', action='store_true', help='Wait for final independently validated source tables')
    parser.add_argument('--rscript', help='Optional Rscript executable to draw the completed summary')
    run(parser.parse_args())

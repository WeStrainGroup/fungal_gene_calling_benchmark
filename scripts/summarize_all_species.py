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


def mean_sd(rows, species, species_by_group=None):
    groups = collections.defaultdict(list)
    for row in rows:
        groups[tuple(row[k] for k in KEY[1:])].append(row)
    result = []
    for key, items in sorted(groups.items()):
        eligible = species_by_group[key[0]] if species_by_group is not None else species
        assert len(items) == len(eligible) and {r['species'] for r in items} == set(eligible)
        for metric in ('precision', 'recall', 'f1'):
            values = [float(r[metric]) for r in items if r[metric] not in ('', 'NA')]
            assert all(math.isfinite(v) for v in values)
            n = len(values)
            result.append(dict(zip(KEY[1:], key), metric=metric, n_species=n,
                total_species=len(eligible), mean=statistics.mean(values) if n else 'NA',
                sd=statistics.stdev(values) if n > 1 else 'NA'))
    return result


def paired_units(observed, expected):
    units = {key[:2] for key in expected}
    return {unit for unit in units if {key for key in expected if key[:2] == unit} <= observed}


def aggregate(roots, manifest, output, paired_conditions=False, groups=GROUPS):
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
    expected = {(s, g, t, scope, level) for s in species for g in groups for t in TOOLS
                for scope in (SCOPES if t == 'carbon' else ('complete',)) for level in LEVELS}
    observed = {tuple(r[k] for k in KEY) for r in metrics}
    assert observed <= expected, 'Unexpected species, condition or metric'
    complete = all(s['status'] == 'COMPLETE' for s in snapshots)
    available_cases = len({tuple(r[k] for k in KEY[:3]) for r in metrics})
    included = species
    species_by_group = {g: species for g in groups}
    if paired_conditions:
        units = paired_units(observed, expected)
        if not units:
            raise ValueError('No species-condition pair has all tools and metrics validated')
        included = [s for s in species if any(unit[0] == s for unit in units)]
        species_by_group = {g: [s for s in species if (s, g) in units] for g in groups}
        metrics = [r for r in metrics if (r['species'], r['group']) in units]
        diagnostics = [r for r in diagnostics if (r['species'], r['group']) in units]
        inputs = [r for r in inputs if (r['species'], r['group']) in units]
    if complete:
        assert observed == expected, 'Missing final metric rows'
    elif not paired_conditions:
        assert not (output/'mean_sd.tsv').exists(), 'Do not replace a final summary with partial inputs'
    for name, rows in [('metrics', metrics), ('diagnostics', diagnostics), ('inputs', inputs)]:
        write_table(output/f'{name}.tsv', rows)
    if complete or paired_conditions:
        write_table(output/'mean_sd.tsv', mean_sd(metrics, included, species_by_group))
    status = dict(aggregation_version='0.2.1', scoring_version='0.2.0',
        status='COMPLETE' if complete else 'INTERIM' if paired_conditions else 'PARTIAL', expected_species=12,
        included_species=len(included) if complete or paired_conditions else 0,
        species=included if complete or paired_conditions else [], available_cases=available_cases,
        species_by_condition=species_by_group if complete or paired_conditions else {},
        interim_rule='All three tools validated within each species-condition pair' if paired_conditions else None,
        completed_cases=len({tuple(r[k] for k in KEY[:3]) for r in metrics}), expected_cases=len(species)*len(groups)*len(TOOLS),
        metric_rows=len(metrics), expected_metric_rows=len(expected),
        duplicate_metric_rows_removed=len(collected['metrics'])-len(observed),
        sample_sd_ddof=1, missing_values='excluded, never replaced with zero',
        sources=snapshots)
    (output/'STATUS.json').write_text(json.dumps(status, indent=2)+'\n')
    return status


def run(args):
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    with (output/'summary.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        status = aggregate(args.inputs, args.species, output, args.paired_conditions, args.groups)
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
                status = aggregate(args.inputs, args.species, output, groups=args.groups)
                print(json.dumps(status, indent=2), flush=True)
        if status['status'] in ('COMPLETE', 'INTERIM') and args.rscript:
            subprocess.run([args.rscript, str(Path(__file__).with_name('plot_all_species.R')),
                            str(output), str(output/'figures')], check=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, nargs='+', required=True,
                        help='Validated ANNEVO, Tiberius and Carbon-A cohort result directories')
    parser.add_argument('--species', type=Path, default=Path(__file__).resolve().parents[1]/'metadata/species.tsv')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--groups', nargs='+', default=GROUPS, choices=(*GROUPS, '2-5kb'))
    parser.add_argument('--wait', action='store_true', help='Wait for final independently validated source tables')
    parser.add_argument('--paired-conditions', action='store_true',
                        help='Interim report restricted to species-condition pairs with all three tools validated')
    parser.add_argument('--rscript', help='Optional Rscript executable to draw the completed summary')
    args = parser.parse_args()
    if args.wait and args.paired_conditions:
        parser.error('Choose either an interim paired-condition report or waiting for the final report')
    run(args)

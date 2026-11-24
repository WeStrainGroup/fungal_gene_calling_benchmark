#!/usr/bin/env python3
"""Summarize actual FASTA record lengths, per species and pooled across species."""
import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

GROUPS = ('genome', '10-100kb', '5-10kb', '2-5kb', '1-2kb')


def fasta_lengths(path):
    sizes, names, current = [], set(), None
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if current is not None:
                    sizes.append(current)
                name = line[1:].split()[0]
                assert name not in names, (path, name)
                names.add(name)
                current = 0
            elif line.strip():
                assert current is not None, path
                current += len(''.join(line.split()))
    if current is not None:
        sizes.append(current)
    assert sizes and min(sizes) > 0, path
    return sizes


def quantile(values, fraction):
    # Linear interpolation, equivalent to R quantile(type = 7).
    index = (len(values) - 1) * fraction
    low, high = math.floor(index), math.ceil(index)
    return values[low] + (values[high] - values[low]) * (index - low)


def summarize(species, group, sizes):
    rows = []
    for scope in ('all', 'scored'):
        values = sorted(x for x in sizes if scope == 'all' or x >= 1000)
        row = dict(species=species, group=group, scope=scope,
                   records=len(values), total_bp=sum(values))
        for name, q in [('min_bp', 0), ('q25_bp', .25), ('median_bp', .5),
                        ('q75_bp', .75), ('max_bp', 1)]:
            row[name] = quantile(values, q) if values else 'NA'
        cumulative = 0
        row['n50_bp'] = 'NA'
        for value in reversed(values):
            cumulative += value
            if 2 * cumulative >= row['total_bp']:
                row['n50_bp'] = value
                break
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--samples', type=Path, required=True)
    parser.add_argument('--output', default='-', help='TSV path; default: stdout')
    args = parser.parse_args()
    with args.samples.open() as handle:
        species = list(dict.fromkeys(r['species'] for r in csv.DictReader(handle, delimiter='\t')))
    rows, pooled = [], defaultdict(list)
    for name in species:
        total = None
        for group in GROUPS:
            sizes = fasta_lengths(args.data/'inputs'/name/f'{group}.fna')
            if total is None:
                total = sum(sizes)
            assert sum(sizes) == total, (name, group)
            rows.extend(summarize(name, group, sizes))
            pooled[group].extend(sizes)
    for group in GROUPS:
        rows.extend(summarize('ALL', group, pooled[group]))
    handle = sys.stdout if args.output == '-' else open(args.output, 'w')
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    if handle is not sys.stdout:
        handle.close()


if __name__ == '__main__':
    main()

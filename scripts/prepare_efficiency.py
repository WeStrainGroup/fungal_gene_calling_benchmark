#!/usr/bin/env python3
"""Prepare original assemblies and length-balanced samples from frozen inputs."""
import argparse
import csv
import hashlib
import random
import shutil
from pathlib import Path

RANGES = {'1-2kb': (1000, 2000), '2-5kb': (2000, 5000),
          '5-10kb': (5000, 10000), '10-100kb': (10000, 100000)}


def fasta(path):
    name, parts = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, ''.join(parts)


def digest(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seed', type=int, default=20261124)
    args = parser.parse_args()
    samples = list(csv.DictReader((args.data/'samples.tsv').open(), delimiter='\t'))
    args.output.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(args.data/'samples.tsv', args.output/'samples.tsv')
    with (args.output/'selected_sequences.tsv').open('w') as manifest, \
            (args.output/'cases.tsv').open('w') as cases:
        writer = csv.writer(manifest, delimiter='\t', lineterminator='\n')
        writer.writerow(['species', 'group', 'record', 'bp', 'length_bin',
                         'sequence_sha256', 'source_fasta', 'seed'])
        summary = csv.writer(cases, delimiter='\t', lineterminator='\n')
        summary.writerow(['species', 'group', 'genetic_code', 'records', 'bp',
                          'min_bp', 'median_bp', 'max_bp', 'fasta_sha256'])
        for sample in samples:
            species = sample['species']
            source = args.data/'inputs'/species
            dest = args.output/'inputs'/species
            dest.mkdir(parents=True)
            shutil.copyfile(source/'codon_table.json', dest/'codon_table.json')
            for group in ['genome', *RANGES]:
                path = source/f'{group}.fna'
                records = list(fasta(path))
                assert len({name for name, seq in records}) == len(records), path
                bins = {}
                if group == 'genome':
                    selected = [(name, seq) for name, seq in records if len(seq) >= 1000]
                else:
                    low, high = RANGES[group]
                    pools = [[] for _ in range(5)]
                    for name, seq in records:
                        if low <= len(seq) <= high:
                            index = min(4, (len(seq)-low)*5//(high-low))
                            pools[index].append((name, seq))
                            bins[name] = index+1
                    seed = hashlib.sha256(f'{args.seed}/{species}/{group}'.encode()).hexdigest()
                    rng = random.Random(seed)
                    selected = []
                    for index, pool in enumerate(pools, 1):
                        if len(pool) < 10:
                            raise ValueError(f'{species}/{group}: bin {index} has only {len(pool)} records')
                        selected.extend(rng.sample(sorted(pool), 10))
                    rng.shuffle(selected)
                    frozen = {r['record']: r for r in csv.DictReader(
                        (source/f'{group}.fragments.tsv').open(), delimiter='\t')}
                    for name, seq in selected:
                        row = frozen[name]
                        assert int(row['bp']) == len(seq) and row['sequence_sha256'] == digest(seq), name
                out = dest/f'{group}.fna'
                with out.open('w') as handle:
                    for name, seq in selected:
                        handle.write(f'>{name}\n{seq}\n')
                        writer.writerow([species, group, name, len(seq), bins.get(name, ''),
                                         digest(seq), str(path.resolve()), args.seed])
                sizes = sorted(len(seq) for name, seq in selected)
                n = len(sizes)
                median = (sizes[(n-1)//2]+sizes[n//2])/2
                summary.writerow([species, group, sample['genetic_code'], n, sum(sizes),
                                  sizes[0], median, sizes[-1], hashlib.sha256(out.read_bytes()).hexdigest()])
                print(f'{species}\t{group}\t{n}\t{sum(sizes)}', flush=True)


if __name__ == '__main__':
    main()

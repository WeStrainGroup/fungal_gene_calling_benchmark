#!/usr/bin/env python3
"""Restore frozen fragment FASTAs from original assemblies and coordinate maps."""
import argparse
import csv
import hashlib
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def verify(path, expected):
    if sha256(path) != expected:
        raise ValueError(f'SHA256 mismatch: {path}')


def read_fasta(path):
    records = {}
    name = None
    parts = []
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    records[name] = ''.join(parts)
                name = line[1:].split()[0]
                if name in records:
                    raise ValueError(f'Duplicate record: {name}')
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        records[name] = ''.join(parts)
    return records


def restore(map_path, genome, output, expected):
    """Write in frozen map order; coordinates are 1-based and inclusive."""
    if output.exists():
        verify(output, expected)
        return False
    temporary = output.with_suffix('.fna.tmp')
    try:
        with map_path.open(newline='') as source, temporary.open('wb') as target:
            seen = set()
            for row in csv.DictReader(source, delimiter='\t'):
                record = row['record']
                start, end = int(row['start']), int(row['end'])
                chrom = genome[row['chrom']]
                if record in seen or not 1 <= start <= end <= len(chrom):
                    raise ValueError(f'Invalid record or coordinates: {record}')
                seen.add(record)
                sequence = chrom[start - 1:end].encode('ascii')
                if len(sequence) != int(row['bp']):
                    raise ValueError(f'Length mismatch: {record}')
                if hashlib.sha256(sequence).hexdigest() != row['sequence_sha256']:
                    raise ValueError(f'Sequence SHA256 mismatch: {record}')
                target.write(b'>' + record.encode('ascii') + b'\n' + sequence + b'\n')
        verify(temporary, expected)
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)
    return True


def reconstruct(data):
    checksums = {}
    for line in (data / 'SHA256SUMS').read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        checksums[name] = digest
    maps = sorted((data / 'inputs').glob('*/*.fragments.tsv'))
    if not maps:
        raise ValueError('No fragment maps found')
    previous = None
    restored = 0
    for map_path in maps:
        if map_path.parent != previous:
            genome_path = map_path.parent / 'genome.fna'
            verify(genome_path, checksums[genome_path.relative_to(data).as_posix()])
            genome = read_fasta(genome_path)
            previous = map_path.parent
        verify(map_path, checksums[map_path.relative_to(data).as_posix()])
        output = map_path.with_name(map_path.name.replace('.fragments.tsv', '.fna'))
        restored += restore(map_path, genome, output,
                            checksums[output.relative_to(data).as_posix()])
    print(f'{len(maps)} fragment FASTAs verified; {restored} reconstructed.')
    return len(maps), restored


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=Path('data'))
    args = parser.parse_args()
    reconstruct(args.data)

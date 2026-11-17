#!/usr/bin/env python3
"""Reconstruct the exact tested inputs from RefSeq DNA and frozen fragment maps."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from Bio import SeqIO

REPO = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data', type=Path, default=REPO/'data')
args = parser.parse_args()
data = args.data
records = list(csv.DictReader((REPO/'results/original_records.tsv').open(), delimiter='\t'))
for species in ('Saccharomyces_cerevisiae', 'Candida_albicans'):
    source = data/'provenance'/species
    checks = {Path(line.split(None, 1)[1]).name: line.split()[0] for line in (source/'md5checksums.txt').read_text().splitlines()}
    for path in source.glob('*_genomic.*.gz'):
        assert hashlib.md5(path.read_bytes()).hexdigest() == checks[path.name], path.name
    with gzip.open(next(source.glob('*_genomic.fna.gz')), 'rt') as handle:
        genome = {r.id: str(r.seq).upper() for r in SeqIO.parse(handle, 'fasta')}
    retained = [r for r in records if r['species'] == species]
    folder = data/'inputs'/species
    (folder/'single').mkdir(exist_ok=True)
    with (folder/'genome.fna').open('w') as output:
        for row in retained:
            record = row['record']
            seq = genome[record]
            assert len(seq) == int(row['bp']) and hashlib.sha256(seq.encode()).hexdigest() == row['sha256']
            text = f'>{record}\n{seq}\n'
            output.write(text)
            (folder/'single'/f'{record}.fna').write_text(text)
    for group in ('1-2kb', '5-10kb', '50-100kb'):
        rows = csv.DictReader((folder/f'{group}.fragments.tsv').open(), delimiter='\t')
        with (folder/f'{group}.fna').open('w') as output:
            for row in rows:
                seq = genome[row['chrom']][int(row['start'])-1:int(row['end'])]
                assert len(seq) == int(row['bp']) and hashlib.sha256(seq.encode()).hexdigest() == row['sequence_sha256']
                output.write(f'>{row["record"]}\n{seq}\n')
checks = json.loads((REPO/'metadata/input_sha256.json').read_text())
for name, digest in checks.items():
    assert hashlib.sha256((data/name).read_bytes()).hexdigest() == digest, name
print('All original and fragmented inputs match the inference-time SHA256 checksums.')

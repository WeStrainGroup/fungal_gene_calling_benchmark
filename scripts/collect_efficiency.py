#!/usr/bin/env python3
"""Collect one complete runtime round and verify its native output files."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--results', type=Path, required=True, help='One round directory containing 17 species')
p.add_argument('--output', type=Path, required=True, help='Combined metrics TSV')
a = p.parse_args()
rows, receipts = [], []
for path in sorted(a.results.glob('*/metrics.tsv')):
    hardware = json.loads((path.parent/'hardware.json').read_text())
    assert hardware['passed'], path.parent
    with path.open() as handle:
        cases = list(csv.DictReader(handle, delimiter='\t'))
    assert len(cases) == 20 and all(r['passed'] == 'True' and r['exit_code'] == '0' for r in cases)
    for r in cases:
        dest = path.parent/r['tool']/r['group']
        gff = dest/('prediction.gff' if r['tool'] == 'carbon' else 'prediction.gff3')
        assert gff.is_file() and gff.stat().st_size and (dest/'time.tsv').is_file()
        receipts.append(dict(species=r['species'], tool=r['tool'], group=r['group'],
            gff_sha256=hashlib.sha256(gff.read_bytes()).hexdigest()))
    rows.extend(cases)
assert len(rows) == len({(r['species'], r['group'], r['tool']) for r in rows}) == 340
assert len({r['species'] for r in rows}) == 17
a.output.parent.mkdir(parents=True, exist_ok=True)
with a.output.open('w') as handle:
    w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
    w.writeheader()
    w.writerows(rows)
a.output.with_suffix('.verified.json').write_text(json.dumps(receipts, indent=2)+'\n')
print(f'Collected {len(rows)} successful cases from 17 species.')

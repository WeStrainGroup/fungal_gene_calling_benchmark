#!/usr/bin/env python3
"""Calculate species-equal means from complete per-species metrics."""
import argparse,collections,csv,statistics
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('input',type=Path)
parser.add_argument('output',type=Path)
args=parser.parse_args()
rows=list(csv.DictReader(args.input.open(),delimiter='\t'))
species={r['species'] for r in rows};groups=collections.defaultdict(list)
for row in rows:groups[(row['tool'],row['group'],row['level'])].append(row)
assert len(groups)==36
results=[]
for (tool,group,level),items in groups.items():
    assert len(items)==len(species) and {r['species'] for r in items}==species, 'Incomplete or duplicate species grid'
    for metric in ('precision','recall','f1'):
        values=[float(r[metric]) for r in items]
        results.append(dict(tool=tool,group=group,level=level,metric=metric,n_species=len(values),
                            mean=statistics.mean(values)))
args.output.parent.mkdir(parents=True,exist_ok=True)
with args.output.open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter='\t');w.writeheader();w.writerows(results)
print(f'{len(species)} equally weighted species; means only; no pooled-base weighting.')

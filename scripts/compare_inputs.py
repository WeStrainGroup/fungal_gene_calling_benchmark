#!/usr/bin/env python3
"""Compare bulk and single-record predictions without relying on generated IDs."""
import argparse
import collections
import csv
import json
from pathlib import Path
from Bio.Seq import Seq
from evaluate import fasta, read_models, sequence

REPO = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data', type=Path, default=REPO/'data')
args = parser.parse_args()
rows = []
for tool in ('annevo', 'tiberius'):
    for species in ('Saccharomyces_cerevisiae', 'Candida_albicans'):
        genome = fasta(args.data/'inputs'/species/'genome.fna')
        folder = args.data/'out'/tool/species
        bulk = read_models(folder/'genome/prediction.gff3', include_stops=True)
        single = []
        for record in genome:
            single.extend(read_models(folder/'single'/record/'prediction.gff3', include_stops=True))
        def signatures(models, include_sequence=False):
            result = collections.Counter()
            for m in models:
                key = (m['record'], m['strand'], tuple(map(tuple, m['blocks'])))
                if include_sequence:
                    seq = sequence(m, genome)
                    code = 12 if species == 'Candida_albicans' else 1
                    coding = seq[m['phase']:]
                    protein = str(Seq(coding[:len(coding)//3*3]).translate(table=code))
                    key += (m['phase'], seq, protein)
                result[key] += 1
            return result
        left, right = signatures(bulk), signatures(single)
        matched = sum((left & right).values())
        seq_matched = sum((signatures(bulk, True) & signatures(single, True)).values())
        assert matched == seq_matched
        rows.append(dict(tool=tool, species=species, bulk_models=len(bulk), single_models=len(single),
                         matching_structures=matched, bulk_only=len(bulk)-matched, single_only=len(single)-matched))
expected = list(csv.DictReader((REPO/'results/equivalence_summary.tsv').open(), delimiter='\t'))
for row in rows:
    old = next(r for r in expected if r['tool'] == row['tool'] and r['species'] == row['species'])
    assert all(str(v) == old[k] for k, v in row.items())
(args.data/'input_equivalence_reproduced.json').write_text(json.dumps(rows, indent=2)+'\n')
print('All four bulk/single comparisons reproduced, including CDS, phase and protein agreement among matching structures.')

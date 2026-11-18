#!/usr/bin/env python3
"""Audit predicted CDS without removing models or changing benchmark scores."""
import argparse,csv,json,sys
from collections import Counter
from pathlib import Path
from Bio.Seq import Seq
REPO=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data',type=Path,default=REPO/'data')
ROOT=parser.parse_args().data.resolve()
import evaluate as e
rows=[]
for species,code in [('Saccharomyces_cerevisiae',1),('Candida_albicans',12)]:
    stops={c for c,a in json.loads((ROOT/'inputs'/species/'codon_table.json').read_text()).items() if a=='*'}
    for group in e.GROUPS:
        genome=e.fasta(ROOT/'inputs'/species/(group+'.fna'))
        for tool in e.TOOLS:
            path=ROOT/'out'/tool/species/group/('prediction.gff' if tool=='carbon' else 'prediction.gff3')
            models=e.read_models(path,include_stops=True)
            counts=Counter(total=len(models))
            for m in models:
                raw=e.sequence(m,genome)
                coding=raw[m['phase']:]
                start=m['phase']==0 and raw.startswith('ATG')
                triplet=len(coding)%3==0
                stop=triplet and coding[-3:] in stops
                aa=str(Seq(coding[:len(coding)//3*3]).translate(table=code))
                codons=[coding[i:i+3] for i in range(0,len(coding)-2,3)]
                internal=any(c in stops for c in (codons[:-1] if stop else codons))
                ambiguous_stop=('*' in (aa[:-1] if stop else aa)) and not internal
                edge=m['blocks'][0][0]<=3 or m['blocks'][-1][1]>=len(genome[m['record']])-2
                counts.update(missing_start=not start,missing_stop=not stop,nontriplet_after_phase=not triplet,
                              nonzero_initial_phase=m['phase']!=0,internal_stop=internal,
                              ambiguous_codon_resolves_to_internal_stop=ambiguous_stop,
                              contains_ambiguous_base=bool(set(raw)-set('ACGT')),
                              codon_complete=start and stop and not internal,
                              incomplete_at_edge=(not start or not stop) and edge,
                              incomplete_away_from_edge=(not start or not stop) and not edge)
            rows.append(dict(species=species,tool=tool,group=group,**counts))
            print(species,tool,group,dict(counts),flush=True)
out=ROOT/'native_codon_audit.tsv'
with out.open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)

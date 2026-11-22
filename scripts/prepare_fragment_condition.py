#!/usr/bin/env python3
"""Add an independent 2–5 kb partition without changing frozen inputs."""
import argparse, csv, hashlib, json, random
from pathlib import Path
from Bio import SeqIO

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--base',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--metadata',type=Path,required=True)
a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=True)
assert not (a.output/'INPUTS_READY').exists(), 'Inputs already frozen; do not regenerate'
rows=[]
for manifest,cohort in [('pilot_samples.tsv','three_predictor_full_genome'),('tiberius_samples.tsv','fungi_test_six'),('carbon_samples.tsv','carbon_fungi_set')]:
 for r in csv.DictReader((a.metadata/manifest).open(),delimiter='\t'):
  if r['species'] not in {x['species'] for x in rows}:rows.append(dict(r,source_cohort=cohort))
assert len(rows)==12
summary=[]; hashes={}; source_hashes={}
for row in rows:
 sp=row['species'];old=a.base/row['source_cohort'];folder=a.output/'inputs'/sp;folder.mkdir(parents=True,exist_ok=True)
 for name in ['genome.fna','codon_table.json']:
  src=old/'inputs'/sp/name;dst=folder/name
  assert not dst.exists();dst.symlink_to(src)
  source_hashes[str(src)]=hashlib.sha256(src.read_bytes()).hexdigest()
 ref=a.output/'provenance'/sp;ref.parent.mkdir(exist_ok=True);ref.symlink_to(old/'provenance'/sp,target_is_directory=True)
 genome={r.id:str(r.seq).upper() for r in SeqIO.parse(folder/'genome.fna','fasta')}
 rng=random.Random(int(row['seed']));fragments=[]
 fasta=folder/'2-5kb.fna'
 with fasta.open('w') as f:
  for chrom,seq in genome.items():
   offset=0
   while offset<len(seq):
    remain=len(seq)-offset
    length=remain if remain<=5000 else rng.randint(2000,min(5000,remain-2000))
    name=f'2_5kb_{len(fragments)+1:06}';part=seq[offset:offset+length]
    f.write(f'>{name}\n{part}\n')
    fragments.append(dict(record=name,chrom=chrom,start=offset+1,end=offset+length,bp=length,sequence_sha256=hashlib.sha256(part.encode()).hexdigest()))
    offset+=length
 mapping=folder/'2-5kb.fragments.tsv'
 with mapping.open('w') as f:
  w=csv.DictWriter(f,fieldnames=fragments[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(fragments)
 # Read back every fragment and reconstruct every original record.
 rebuilt={name:[] for name in genome};ends={name:0 for name in genome}
 observed=list(SeqIO.parse(fasta,'fasta'));assert len(observed)==len(fragments)
 for r,m in zip(observed,fragments):
  seq=str(r.seq);assert r.id==m['record'] and m['start']==ends[m['chrom']]+1
  assert len(seq)==m['bp'] and hashlib.sha256(seq.encode()).hexdigest()==m['sequence_sha256']
  assert seq==genome[m['chrom']][m['start']-1:m['end']]
  assert 2000<=len(seq)<=5000 or len(seq)==len(genome[m['chrom']])<2000
  rebuilt[m['chrom']].append(seq);ends[m['chrom']]=m['end']
 assert all(''.join(rebuilt[k])==v for k,v in genome.items())
 for path in [folder/'genome.fna',folder/'codon_table.json',fasta,mapping]:hashes[str(path.relative_to(a.output))]=hashlib.sha256(path.read_bytes()).hexdigest()
 summary.append(dict(species=sp,seed=row['seed'],records=len(fragments),bp=sum(x['bp'] for x in fragments),below_2kb_records=sum(x['bp']<2000 for x in fragments),below_1kb_records=sum(x['bp']<1000 for x in fragments),below_1kb_bp=sum(x['bp'] for x in fragments if x['bp']<1000),reconstruction='PASS'))
 print(sp,len(fragments),summary[-1]['bp'],flush=True)
for name,items in [('samples.tsv',rows),('inputs.tsv',summary)]:
 with (a.output/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=items[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(items)
for path,sha in source_hashes.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==sha
(a.output/'input_sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
(a.output/'preparation.json').write_text(json.dumps(dict(status='PASS',condition='2-5kb',species=12,expected_cases=36,seed_rule='Each species frozen seed, independent RNG for this condition',source_sha256=source_hashes),indent=2)+'\n')
(a.output/'INPUTS_READY').write_text('12 genomes reconstructed from all 2–5 kb fragments; source hashes unchanged.\n')

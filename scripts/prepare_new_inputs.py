#!/usr/bin/env python3
"""Fetch fixed RefSeq assemblies and verify current nuclear sequence partitions."""
import argparse,csv,gzip,hashlib,json,random,re,subprocess,urllib.request
from pathlib import Path
from Bio import SeqIO
from Bio.Data import CodonTable
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data',type=Path,required=True)
parser.add_argument('--samples',type=Path,required=True)
args=parser.parse_args()
ROOT=args.data.resolve();ROOT.mkdir(parents=True,exist_ok=True)
assert not (ROOT/'INPUTS_READY').exists(), 'Inputs already frozen; reuse them'
rows=list(csv.DictReader(args.samples.open(),delimiter='\t'))
records,summary,downloads=[],[],[]
def fetch(url,path):
    if path.is_file() and path.stat().st_size:return
    temporary=path.with_name(path.name+'.part')
    subprocess.run(['curl','-fLsS','--retry','3','--connect-timeout','15','--max-time','600',url,'-o',str(temporary)],check=True)
    temporary.rename(path)
def table(path,rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
for row in rows:
    species,accession=row['species'],row['accession'];dest=ROOT/'provenance'/species;dest.mkdir(parents=True,exist_ok=True)
    digits=accession.split('_')[1].split('.')[0]
    parent='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/'+'/'.join(digits[i:i+3] for i in (0,3,6))+'/'
    listing=subprocess.check_output(['curl','-fLsS','--retry','3','--max-time','60',parent],text=True)
    names=re.findall(r'href="('+re.escape(accession)+r'_[^"/]+/)"',listing)
    assert len(names)==1,(accession,names)
    assembly=names[0].rstrip('/');url=parent+names[0]
    for suffix in ('_genomic.fna.gz','_genomic.gff.gz','_assembly_report.txt'):
        fetch(url+assembly+suffix,dest/(assembly+suffix))
    fetch(url+'md5checksums.txt',dest/'md5checksums.txt')
    checks={Path(line.split(None,1)[1]).name:line.split()[0] for line in (dest/'md5checksums.txt').read_text().splitlines()}
    for path in dest.glob(assembly+'*'):
        assert hashlib.md5(path.read_bytes()).hexdigest()==checks[path.name],path
    excluded=set()
    for line in (dest/(assembly+'_assembly_report.txt')).read_text().splitlines():
        if line.startswith('#'):continue
        c=line.split('\t')
        if any(word in (c[3]+' '+c[7]).lower() for word in ('mitochond','plastid','chloroplast','non-nuclear')):
            excluded.update((c[4],c[6]))
    with gzip.open(dest/(assembly+'_genomic.fna.gz'),'rt') as f:
        genome={r.id:str(r.seq).upper() for r in SeqIO.parse(f,'fasta') if r.id not in excluded}
    # Confirm the nuclear translation table against explicit CDS qualifiers.
    codes=set()
    with gzip.open(dest/(assembly+'_genomic.gff.gz'),'rt') as f:
        for line in f:
            if not line.startswith('#') and line.split('\t')[0] in genome:
                codes.update(int(x) for x in re.findall(r'transl_table=(\d+)',line))
    assert not codes or codes=={int(row['genetic_code'])},(species,codes)
    folder=ROOT/'inputs'/species;folder.mkdir(parents=True,exist_ok=True)
    with (folder/'genome.fna').open('w') as f:
        for name,seq in genome.items():
            f.write(f'>{name}\n{seq}\n');records.append(dict(species=species,record=name,bp=len(seq),sha256=hashlib.sha256(seq.encode()).hexdigest()))
    code=CodonTable.unambiguous_dna_by_id[int(row['genetic_code'])]
    (folder/'codon_table.json').write_text(json.dumps(dict(code.forward_table,**{x:'*' for x in code.stop_codons}),indent=2)+'\n')
    for group,minimum,maximum in (('1-2kb',1000,2000),('2-5kb',2000,5000),('5-10kb',5000,10000),('10-100kb',10000,100000)):
        rng=random.Random(int(row['seed']))
        fragments=[]
        with (folder/f'{group}.fna').open('w') as f:
            for chrom,sequence in genome.items():
                offset=0
                while offset<len(sequence):
                    remaining=len(sequence)-offset
                    length=remaining if remaining<=maximum else rng.randint(minimum,min(maximum,remaining-minimum))
                    name=f'{group.replace("-","_")}_{len(fragments)+1:06}';start,end=offset+1,offset+length;seq=sequence[offset:end]
                    f.write(f'>{name}\n{seq}\n');fragments.append(dict(record=name,chrom=chrom,start=start,end=end,bp=length,sequence_sha256=hashlib.sha256(seq.encode()).hexdigest()));offset=end
        assert sum(r['bp'] for r in fragments)==sum(map(len,genome.values()))
        table(folder/f'{group}.fragments.tsv',fragments)
        rebuilt={name:[] for name in genome};ends={name:0 for name in genome}
        observed=list(SeqIO.parse(folder/f'{group}.fna','fasta'))
        assert len(observed)==len(fragments)
        for rec,m in zip(observed,fragments):
            seq=str(rec.seq);assert rec.id==m['record'] and m['start']==ends[m['chrom']]+1
            assert seq==genome[m['chrom']][m['start']-1:m['end']]
            assert hashlib.sha256(seq.encode()).hexdigest()==m['sequence_sha256']
            assert minimum<=len(seq)<=maximum or len(seq)==len(genome[m['chrom']])<minimum
            rebuilt[m['chrom']].append(seq);ends[m['chrom']]=m['end']
        assert all(''.join(rebuilt[k])==v for k,v in genome.items())
        summary.append(dict(species=species,group=group,records=len(fragments),bp=sum(r['bp'] for r in fragments),min_bp=min(r['bp'] for r in fragments),max_bp=max(r['bp'] for r in fragments),below_target_min_records=sum(r['bp']<minimum for r in fragments),below_target_min_bp=sum(r['bp'] for r in fragments if r['bp']<minimum)))
    summary.append(dict(species=species,group='genome',records=len(genome),bp=sum(map(len,genome.values())),min_bp=min(map(len,genome.values())),max_bp=max(map(len,genome.values())),below_target_min_records=0,below_target_min_bp=0))
    downloads.append(dict(species=species,accession=accession,source_url=url,genetic_code=int(row['genetic_code']),explicit_cds_translation_tables=sorted(codes),excluded_records=sorted(excluded)))
    print('READY',species,len(genome),sum(map(len,genome.values())),flush=True)
table(ROOT/'original_records.tsv',records);table(ROOT/'inputs.tsv',summary)
(ROOT/'source_metadata.json').write_text(json.dumps(downloads,indent=2)+'\n')
hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'inputs').rglob('*') if p.is_file()}
(ROOT/'input_sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
(ROOT/'INPUTS_READY').write_text('Fixed public RefSeq assemblies; complete nuclear partitions.\n')

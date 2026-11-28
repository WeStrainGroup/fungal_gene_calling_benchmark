#!/usr/bin/env python3
"""Independently compare a CTC reconstruction with the selected native models."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
from urllib.parse import unquote
from Bio import SeqIO
from Bio.Seq import Seq


def read_gff(path):
    rows=[]
    for line in path.read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        fields=line.split('\t')
        assert len(fields)==9, path
        attrs=dict(item.split('=',1) for item in fields[8].split(';') if '=' in item)
        rows.append((fields,attrs))
    return rows


def parents(attrs):
    return [unquote(p) for p in attrs.get('Parent','').split(',') if p]


def signature(fields):
    return tuple('mRNA' if i==2 and value=='transcript' else value for i,value in enumerate(fields[:8]))


def sequences(path):
    result={}
    with path.open() as handle:
        for record in SeqIO.parse(handle,'fasta'):
            assert record.id not in result, f'Duplicate FASTA ID: {record.id}'
            result[record.id]=str(record.seq)
    return result


def verify(sources, combined, code=1, threshold=0.1):
    expected={}
    expected_cds_rows=Counter()
    for tool,directory in sources:
        rows=read_gff(directory/'prediction.gff3')
        genes={unquote(a['ID']):(f,a) for f,a in rows if f[2]=='gene'}
        transcripts={unquote(a['ID']):(f,a) for f,a in rows if f[2] in {'mRNA','transcript'}}
        by_gene={}
        for transcript,(_,attrs) in transcripts.items():
            for gene in parents(attrs):
                by_gene.setdefault(gene,[]).append(transcript)
        allowed=set(genes)
        if tool=='carbon_a':
            allowed={g for g,(_,a) in genes.items() if a['status']=='complete' and float(a['confidence'])>=threshold}
        cds_rows={t:[] for t in transcripts}
        for fields,attrs in rows:
            if fields[2]=='CDS':
                for transcript in parents(attrs):
                    cds_rows[transcript].append(fields)
                if all(set(parents(transcripts[t][1])) <= allowed for t in parents(attrs)):
                    expected_cds_rows[(signature(fields),tuple(sorted((tool,t) for t in parents(attrs))))]+=1
        native_proteins=sequences(directory/'proteins.fa')
        for original,dna in sequences(directory/'cds.fa').items():
            assert original in native_proteins
            if tool=='tiberius':
                transcript=original if original in transcripts else original.split('|')[1]
            else:
                candidates=by_gene.get(original,[])
                transcript=original if original in transcripts else candidates[0]
                assert original in transcripts or len(candidates)==1
            fields,attrs=transcripts[transcript]
            ancestor=set(parents(attrs))
            assert ancestor <= set(genes)
            if not ancestor <= allowed:
                continue
            segments=cds_rows[transcript]
            first=sorted(segments,key=lambda f:int(f[3]),reverse=fields[6]=='-')[0]
            phase=int(first[7])
            oriented=dna.upper().replace('U','T')[phase:]
            usable=oriented[:3*(len(oriented)//3)]
            protein=str(Seq(usable).translate(table=code,cds=False,to_stop=False))
            if protein.endswith('*'):
                protein=protein[:-1]
            expected[tool+':'+original]=dict(tool=tool,transcript=transcript,fields=signature(fields),
                genes={g:signature(genes[g][0]) for g in ancestor},segments=Counter(signature(f) for f in segments),
                dna=dna,protein=protein)
    with (combined/'id_map.tsv').open() as handle:
        reader=csv.DictReader(handle,delimiter='\t')
        assert reader.fieldnames==['new_id','original_id']
        mapping=list(reader)
    new_to_old={r['new_id']:r['original_id'] for r in mapping}
    assert len(new_to_old)==len(mapping)==len(expected)
    assert set(new_to_old.values())==set(expected)
    dna=sequences(combined/'cds.fa');protein=sequences(combined/'proteins.fa')
    assert set(dna)==set(protein)==set(new_to_old)
    for name in ['cds.fa','proteins.fa']:
        with (combined/name).open() as handle:
            assert all(r.description==r.id for r in SeqIO.parse(handle,'fasta'))
    rows=read_gff(combined/'prediction.gff3')
    assert all(f[2] in {'gene','mRNA','CDS'} and set(a)<= {'ID','Parent'} for f,a in rows)
    features={}
    for fields,attrs in rows:
        features.setdefault(attrs['ID'],[]).append((fields,attrs))
    for key,parts in features.items():
        assert len({(f[0],f[2],f[6],a.get('Parent')) for f,a in parts})==1
        assert len({(f[3],f[4]) for f,a in parts})==len(parts)
        if parts[0][0][2]!='CDS':
            assert len(parts)==1
        assert all(p in features for _,a in parts for p in parents(a))
    assert set(new_to_old)=={a['ID'] for f,a in rows if f[2]=='mRNA'}
    segments={t:[] for t in new_to_old};observed_cds_rows=Counter()
    for fields,attrs in rows:
        if fields[2]=='CDS':
            for t in parents(attrs):
                segments[t].append(signature(fields))
            original_parents=tuple(sorted((expected[new_to_old[t]]['tool'],expected[new_to_old[t]]['transcript']) for t in parents(attrs)))
            observed_cds_rows[(signature(fields),original_parents)]+=1
    assert expected_cds_rows==observed_cds_rows, 'CDS rows or their Parent associations changed'
    gene_mapping={}
    for new,old in new_to_old.items():
        e=expected[old];fields,attrs=features[new][0]
        assert signature(fields)==e['fields'], new
        assert Counter(segments[new])==e['segments'], new
        assert dna[new]==e['dna'] and protein[new]==e['protein'], new
        native_genes=e['genes']
        combined_genes={g:signature(features[g][0][0]) for g in parents(attrs)}
        assert Counter(native_genes.values())==Counter(combined_genes.values()), new
        # Normal native models have one gene; retain isoform-to-gene associations.
        if len(native_genes)==1:
            origin=(e['tool'],next(iter(native_genes)))
            renamed=next(iter(combined_genes))
            assert gene_mapping.setdefault(origin,renamed)==renamed
    assert len(set(gene_mapping.values()))==len(gene_mapping)
    assert {a['ID'] for f,a in rows if f[2]=='gene'}=={p for f,a in rows if f[2]=='mRNA' for p in parents(a)}
    return dict(models=len(mapping),genes=sum(f[2]=='gene' for f,a in rows),cds_rows=sum(f[2]=='CDS' for f,a in rows),passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tiberius',type=Path)
    p.add_argument('--carbon-a',type=Path)
    p.add_argument('--combined',type=Path,required=True)
    p.add_argument('--codon-table',type=int,default=1)
    p.add_argument('--confidence-threshold',type=float,default=.1)
    a=p.parse_args()
    sources=[(tool,directory) for tool,directory in [('tiberius',a.tiberius),('carbon_a',a.carbon_a)] if directory]
    if not sources:
        p.error('Supply at least one native output directory')
    print(json.dumps(verify(sources,a.combined,a.codon_table,a.confidence_threshold)))

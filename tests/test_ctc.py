"""Synthetic regression tests for CTC routing, translation and GFF3 model merging."""
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('ctc', ROOT/'scripts/ctc.py')
ctc = importlib.util.module_from_spec(spec);spec.loader.exec_module(ctc)


def feature(seq, kind, start, end, attrs, strand='+', phase='.'):
    return f'{seq}\tPredictor\t{kind}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}\n'


def model(seq, gene, transcript, start=1, status='', strand='+'):
    extra=';'+status if status else ''
    return (feature(seq,'gene',start,start+8,f'ID={gene};Name={gene}{extra}',strand)+
            feature(seq,'mRNA',start,start+8,f'ID={transcript};Parent={gene}',strand)+
            feature(seq,'CDS',start,start+8,f'ID={transcript}.part1;Parent={transcript}',strand,'0'))


def native(root, tool, gff, records):
    directory=root/tool;directory.mkdir()
    (directory/'prediction.gff3').write_text('##gff-version 3\n'+gff)
    for filename,col in [('cds.fa',1),('proteins.fa',2)]:
        (directory/filename).write_text(''.join(f'>{r[0]}\n{r[col]}\n' for r in records))
    return directory


def id_map(directory):
    with (directory/'id_map.tsv').open() as handle:
        rows=csv.DictReader(handle,delimiter='\t')
        assert rows.fieldnames==['new_id','original_id']
        return {r['new_id']:r['original_id'] for r in rows}


class CTCTests(unittest.TestCase):
    def test_boundaries_and_duplicate_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for tool in ['tiberius','carbon_a','combined']:(root/tool).mkdir()
            f=root/'input.fa'
            f.write_text(''.join(f'>r{i} description\n'+('A'*n)+'\n' for i,n in enumerate([999,1000,10000,10001])))
            counts=ctc.split_fasta(f,root,1000,10000)
            self.assertEqual(counts,{'tiberius':{'records':2,'bp':11000},'carbon_a':{'records':1,'bp':10001},'skipped':{'records':1,'bp':999}})
            f.write_text('>dup one\nA\n>dup two\nA\n')
            with self.assertRaisesRegex(ValueError,'Duplicate'):ctc.split_fasta(f,root,1,10)

    def test_genetic_codes_stops_and_partial_codons(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'cds.fa'
            source.write_text('>g1|t1|chr:0-9(+)\nATGCTGTAA\n>t2\nATGTAANNN\n>t3\nATGA\n>t4\nATGAA\n>t5\nGCNTAA\n>t6\nATGTAATAA\n')
            for code,expected in [('1','ML'),('12','MS'),('26','MA')]:
                summary=ctc.translate_cds(source,root/'proteins.fa',root/'check.tsv',int(code))
                proteins=dict(ctc.fasta(root/'proteins.fa'))
                self.assertEqual(proteins['g1|t1|chr:0-9(+)'],expected)
                self.assertEqual(proteins['t2'],'M*X')
                self.assertEqual(proteins['t3'],'M')
                self.assertEqual(proteins['t4'],'M')
                self.assertEqual(proteins['t5'],'A')
                self.assertEqual(proteins['t6'],'M*')
                self.assertEqual(summary['internal_stops'],2)
                self.assertEqual(summary['partial_codons'],2)
            source.write_text('>ambiguous\nRAY\n')
            ctc.translate_cds(source,root/'proteins.fa',root/'check.tsv')
            self.assertEqual(dict(ctc.fasta(root/'proteins.fa'))['ambiguous'],'B')

    def test_ncbi_codes_and_native_limits(self):
        supported=[1,3,4,5,6,9,10,11,12,13,14,15,16,21,24,25,26,29,30,32,33]
        for code in supported:
            table=ctc.codon_table(code)
            self.assertEqual(len(table),64)
            self.assertTrue({c for c,a in table.items() if a=='*'} <= {'TAA','TAG','TGA'})
        for code,codon,expected in [(3,'CTG','T'),(4,'TGA','W'),(6,'TAA','Q'),
                (10,'TGA','C'),(25,'TGA','G'),(29,'TAA','Y'),(30,'TAG','E'),(32,'TAG','W')]:
            self.assertEqual(ctc.codon_table(code)[codon],expected)
        for code in [2,22,23,27,28,31]:
            with self.assertRaisesRegex(ValueError,'Carbon-A'):ctc.codon_table(code)
        for code in [0,7,8,17,34,-1]:
            with self.assertRaisesRegex(ValueError,'Unknown or discontinued'):ctc.codon_table(code)

    def test_complete_selection_native_ids_and_incomplete_exclusion(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'combined';out.mkdir()
            t=native(root,'tiberius',model('b','old1','tx1')+model('a','old2','tx2',30),
                [('old1|tx1|b:0-9(+) description','ATGTAATAA','M*'),('old2|tx2|a:29-38(+)','ATGNNNTAA','MX')])
            c=native(root,'carbon_a',model('a','gene1','mrna1',1,'status=complete;confidence=1')+
                model('a','gene2','mrna2',12,'status=incomplete;confidence=1')+
                model('a','gene3','mrna3',50,'status=complete;confidence=0.001'),
                [('gene1 status=complete','ATGCTGTAA','MS'),('gene2 status=incomplete','ATGCTG','MS'),('gene3','ATGCTGTAA','MS')])
            before={p:p.read_bytes() for d in [t,c] for p in d.iterdir()}
            stats=ctc.combine_models([('tiberius',t),('carbon_a',c)],out)
            self.assertEqual(stats['models'],3);self.assertEqual(stats['carbon_a_confidence_excluded'],1);self.assertEqual(stats['carbon_a_incomplete_excluded'],1)
            self.assertEqual(id_map(out),{'g1.t1':'carbon_a:gene1','g2.t1':'tiberius:old2|tx2|a:29-38(+)','g3.t1':'tiberius:old1|tx1|b:0-9(+)'})
            expected=['g1.t1','g2.t1','g3.t1']
            self.assertEqual(list(ctc.fasta_index(out/'cds.fa')),expected)
            self.assertEqual(list(ctc.fasta_index(out/'proteins.fa')),expected)
            self.assertIn('>g3.t1\n',(out/'proteins.fa').read_text())
            self.assertIn('M*',(out/'proteins.fa').read_text())
            self.assertIn('confidence=1',(out/'prediction.gff3').read_text())
            self.assertNotIn('status=incomplete',(out/'prediction.gff3').read_text())
            self.assertTrue(all(p.read_bytes()==data for p,data in before.items()))
            stats=ctc.combine_models([('tiberius',t),('carbon_a',c)],out,confidence_threshold=0)
            self.assertEqual(stats['models'],4);self.assertEqual(stats['carbon_a_incomplete_excluded'],1)
            self.assertNotIn('carbon_a:gene2',id_map(out).values())
            self.assertTrue(all(p.read_bytes()==data for p,data in before.items()))

    def test_confidence_thresholds_na_and_invalid_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'combined';out.mkdir()
            scores=['0.0999','0.1','0.1001','NA','NA']
            gff=''.join(model('x',f'g{i}',f't{i}',1+12*i,
                f"status={'incomplete' if i>=3 else 'complete'};confidence={score}")
                for i,score in enumerate(scores))
            c=native(root,'carbon_a',gff,[(f'g{i}','ATGCTGTAA','ML') for i in range(5)])
            before={p:p.read_bytes() for p in c.iterdir()}
            for threshold,expected in [(.1,{'g1','g2'}),(0,{'g0','g1','g2'}),(1,set())]:
                stats=ctc.combine_models([('carbon_a',c)],out,confidence_threshold=threshold)
                mapping=id_map(out)
                self.assertEqual(set(mapping.values()),{'carbon_a:'+x for x in expected})
                self.assertEqual(set(ctc.fasta_index(out/'cds.fa')),set(mapping))
                self.assertEqual(set(ctc.fasta_index(out/'proteins.fa')),set(mapping))
                rows,_=ctc.read_gff(out/'prediction.gff3')
                self.assertEqual(len({r['attrs']['ID'] for r in rows if r['fields'][2]=='gene'}),len(expected))
                self.assertEqual(stats['models']+stats['carbon_a_incomplete_excluded']+
                    stats['carbon_a_confidence_excluded'],5)
            self.assertTrue(all(p.read_bytes()==data for p,data in before.items()))
            for bad in [-.1,1.1,float('nan'),float('inf')]:
                with self.assertRaisesRegex(ValueError,'confidence-threshold'):
                    ctc.combine_models([('carbon_a',c)],out,confidence_threshold=bad)
            for bad in ['NA','nan','inf','-0.1','1.1','bad',None]:
                attrs='status=complete'+('' if bad is None else ';confidence='+bad)
                (c/'prediction.gff3').write_text(model('x','g0','t0',1,attrs))
                for filename,sequence in [('cds.fa','ATGCTGTAA'),('proteins.fa','ML')]:
                    (c/filename).write_text('>g0\n'+sequence+'\n')
                with self.assertRaises(ValueError):ctc.combine_models([('carbon_a',c)],out)

    def test_native_discontinuous_features_and_shared_references(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'combined';out.mkdir()
            gff=feature('x','gene',1,30,'ID=gene;Name=Biological_name','-')
            for transcript in ['a','b']:
                gff+=feature('x','mRNA',1,30,f'ID={transcript};Parent=gene','-')
            for start,end,phase in [(1,3,'2'),(20,25,'0')]:
                gff+=feature('x','CDS',start,end,'ID=shared;Parent=a,b;confidence=NA','-',phase)
            gff+=feature('x','exon',20,25,'ID=exon;Parent=a,b;Derives_from=shared','-')
            t=native(root,'tiberius',gff,[('gene|a|x:0-30(-)','ATGCTGTAA','ML'),('gene|b|x:0-30(-)','ATGCTGTAA','ML')])
            ctc.combine_models([('tiberius',t)],out)
            rows,features=ctc.read_gff(out/'prediction.gff3')
            self.assertEqual(len(features['shared_cds1']),2)
            self.assertEqual(features['shared_exon1'][0]['attrs']['Parent'],'g1.t1,g1.t2')
            self.assertEqual(features['shared_exon1'][0]['attrs']['Derives_from'],'shared_cds1')
            self.assertEqual(features['g1'][0]['attrs']['Name'],'Biological_name')
            self.assertEqual(id_map(out),{'g1.t1':'tiberius:gene|a|x:0-30(-)','g1.t2':'tiberius:gene|b|x:0-30(-)'})
            self.assertEqual([r['fields'][7] for r in features['shared_cds1']],['2','0'])

    def test_unified_translation_preserves_native_outputs_and_excludes_partial(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'combined';out.mkdir()
            t=native(root,'tiberius',model('x','gT','tT'),[('gT|tT|x:0-9(+)','ATGCTGTAA','ML')])
            gff=model('y','gC','tC',1,'status=complete;confidence=1')
            for tag,strand in [('p','+'),('m','-')]:
                gff+=feature(tag,'gene',1,25,f'ID=g{tag};status=incomplete;confidence=NA',strand)
                gff+=feature(tag,'mRNA',1,25,f'ID=t{tag};Parent=g{tag}',strand)
                # The first CDS in transcript direction has phase 2, regardless of row order.
                for start,end,phase in ([(1,5,'2'),(20,25,'0')] if strand=='+' else [(1,6,'0'),(21,25,'2')]):
                    gff+=feature(tag,'CDS',start,end,f'ID=t{tag}.cds;Parent=t{tag}',strand,phase)
            c=native(root,'carbon_a',gff,[('gC description','ATGGCNTAA','MX'),
                ('gp partial','CCATGCTGTAA','MS'),('gm partial','CCATGCTGTAA','MS')])
            before={p:p.read_bytes() for d in [t,c] for p in d.iterdir()}
            stats=ctc.combine_models([('tiberius',t),('carbon_a',c)],out,code=12)
            proteins=ctc.fasta_index(out/'proteins.fa')
            self.assertEqual({k:v[1] for k,v in proteins.items()},
                {'g1.t1':'MS','g2.t1':'MA'})
            with (out/'translation.tsv').open() as handle:
                checks={r['sequence']:r for r in csv.DictReader(handle,delimiter='\t')}
            self.assertEqual(set(checks),{'g1.t1','g2.t1'})
            self.assertEqual(stats['translation']['records'],2)
            self.assertEqual(stats['carbon_a_incomplete_excluded'],2)
            self.assertTrue(all(p.read_bytes()==data for p,data in before.items()))
            self.assertEqual([seq for _,seq in ctc.fasta(out/'cds.fa')],['ATGCTGTAA','ATGGCNTAA'])

    def test_translation_initial_phase_and_trailing_residues(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'cds.fa'
            source.write_text('>one\nAATGCTGTAAA\n>two\nCCATGCTGTAA\n')
            stats=ctc.translate_cds(source,root/'proteins.fa',root/'check.tsv',12,{'one':1,'two':2})
            self.assertEqual(dict(ctc.fasta(root/'proteins.fa')),{'one':'MS','two':'MS'})
            self.assertEqual(stats['partial_codons'],1)
            with self.assertRaisesRegex(ValueError,'initial CDS phase'):
                ctc.translate_cds(source,root/'proteins.fa',root/'check.tsv',12,{'one':3,'two':2})
            with self.assertRaisesRegex(ValueError,'unmatched'):
                ctc.translate_cds(source,root/'proteins.fa',root/'check.tsv',12,{'one':1})

    def test_native_collisions_are_separated_in_combined_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'combined';out.mkdir()
            t=native(root,'tiberius',model('x','g','tx'),[('g|tx|x:0-9(+)','ATGCTGTAA','ML')])
            c=native(root,'carbon_a',model('y','g','mrna',1,'status=complete;confidence=1'),[('g','ATGCTGTAA','ML')])
            ctc.combine_models([('tiberius',t),('carbon_a',c)],out)
            self.assertEqual(id_map(out),{'g1.t1':'tiberius:g|tx|x:0-9(+)','g2.t1':'carbon_a:g'})
            _,features=ctc.read_gff(out/'prediction.gff3')
            self.assertEqual(set(features),{'g1','g1.t1','g1.t1.cds','g2','g2.t1','g2.t1.cds'})

    def test_record_order_discontinuous_cds_and_exon_numbering(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'combined';out.mkdir()
            gff=model('a','late','late_tx')
            gff+=feature('z','gene',1,30,'ID=early;Name=early','-')
            gff+=feature('z','mRNA',1,30,'ID=early_tx;Parent=early','-')
            for index,start,end,phase in [(1,1,3,'2'),(2,20,25,'0')]:
                gff+=feature('z','CDS',start,end,f'ID=old_cds{index};Parent=early_tx','-',phase)
                gff+=feature('z','exon',start,end,f'Parent=early_tx','-')
            t=native(root,'tiberius',gff,[('late|late_tx|a:0-9(+)','ATGCTGTAA','ML'),
                ('early|early_tx|z:0-30(-)','ATGCTGTAA','ML')])
            ctc.combine_models([('tiberius',t)],out,record_order={'z':0,'a':1})
            _,features=ctc.read_gff(out/'prediction.gff3')
            self.assertEqual(features['g1'][0]['fields'][0],'z')
            self.assertEqual(features['g1'][0]['attrs']['Name'],'g1')
            self.assertEqual(len(features['g1.t1.cds']),2)
            self.assertEqual(features['g1.t1.exon1'][0]['fields'][3:5],['20','25'])
            self.assertEqual(features['g1.t1.exon2'][0]['fields'][3:5],['1','3'])
            self.assertEqual(id_map(out)['g1.t1'],'tiberius:early|early_tx|z:0-30(-)')
            before={p.name:p.read_bytes() for p in out.iterdir()}
            ctc.combine_models([('tiberius',t)],out,record_order={'z':0,'a':1})
            self.assertEqual(before,{p.name:p.read_bytes() for p in out.iterdir()})

    def test_invalid_references_and_duplicate_features_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'test.gff3'
            path.write_text(feature('x','CDS',1,9,'ID=c;Parent=absent','+','0'))
            with self.assertRaisesRegex(ValueError,'Dangling'):ctc.read_gff(path)
            path.write_text(feature('x','gene',1,9,'ID=g')+feature('x','mRNA',1,9,'ID=g'))
            with self.assertRaisesRegex(ValueError,'Duplicate'):ctc.read_gff(path)

    def test_pipeline_empty_subsets_numeric_code_and_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);binpath=root/'env/bin';binpath.mkdir(parents=True)
            mock=binpath/'python'
            mock.write_text('#!'+sys.executable+'\n'+'''import json,pathlib,sys
args=sys.argv
if '--fail' in args:sys.exit(7)
tiberius='--genome' in args
f=pathlib.Path(args[args.index('--genome' if tiberius else '--input')+1])
p=pathlib.Path(args[args.index('--out')+1]) if tiberius else pathlib.Path(args[args.index('--output-dir')+1])/'prediction.gff'
p.parent.mkdir(parents=True,exist_ok=True)
gene,tx=('gT','tT') if tiberius else ('gC','tC')
seq=next(l[1:].split()[0] for l in f.read_text().splitlines() if l.startswith('>'))
p.write_text('##gff-version 3\\n'+seq+'\\tmock\\tgene\\t1\\t9\\t.\\t+\\t.\\tID='+gene+';status=complete;confidence=0.1\\n'+seq+'\\tmock\\tmRNA\\t1\\t9\\t.\\t+\\t.\\tID='+tx+';Parent='+gene+'\\n'+seq+'\\tmock\\tCDS\\t1\\t9\\t.\\t+\\t0\\tID='+tx+'.cds;Parent='+tx+'\\n')
if tiberius:
 pathlib.Path(args[args.index('--codingseq')+1]).write_text('>gT|tT|chr:0-9(+) description\\nATGCTGTAA\\n')
 pathlib.Path(args[args.index('--protseq')+1]).write_text('>gT|tT|chr:0-9(+) description\\nML\\n')
else:
 table=json.loads(pathlib.Path(args[args.index('--codon-table')+1]).read_text())
 assert args[args.index('--min-record-length')+1]=='0'
 (p.parent/'prediction.fna').write_text('>gC status=complete\\nATGCTGTAA\\n')
 (p.parent/'prediction.faa').write_text('>gC status=complete\\nM'+table['CTG']+'\\n')
''');mock.chmod(0o755)
            pipeline=root/'carbon/annotation_pipeline';pipeline.mkdir(parents=True)
            (pipeline/'run_annotation.sh').write_text('#!/bin/sh\npython "$@"\n')
            for label,lengths in [('mixed',[1000,10001]),('short',[1000]),('long',[10001]),('skipped',[999]),('failure',[1000])]:
                f=root/(label+'.fa');f.write_text(''.join(f'>r{i}\n'+('A'*n)+'\n' for i,n in enumerate(lengths)))
                out=root/label
                command=[sys.executable,str(ROOT/'scripts/ctc.py'),'--input',str(f),'--output',str(out),'--cutoff-length','10000',
                    '--tiberius-path',str(mock),'--tiberius-model',str(root/'model'),'--carbon-a-path',str(mock),'--carbon-a-model',str(root/'carbon'), '--codon-table','12']
                if label=='failure':command+=['--tiberius-args=--fail']
                result=subprocess.run(command,capture_output=True,text=True)
                report=json.loads((out/'combined/run.json').read_text())
                self.assertEqual(report['confidence_threshold'],0.1)
                if label=='failure':
                    self.assertNotEqual(result.returncode,0);self.assertEqual(report['status'],'failed')
                else:
                    self.assertEqual(result.returncode,0,result.stderr)
                    self.assertEqual({p.name for p in out.iterdir()},{'tiberius','carbon_a','combined'})
                    self.assertEqual(report['combined']['models'],sum(n>=1000 for n in lengths))
                    if report['routing']['tiberius']['records']:
                        self.assertEqual((out/'tiberius/proteins.fa').read_text(),'>gT|tT|chr:0-9(+) description\nML\n')
                        self.assertEqual(ctc.fasta_index(out/'combined/proteins.fa')['g1.t1'][1],'MS')
                    self.assertTrue((out/'combined/translation.tsv').exists())
                    self.assertFalse((out/'tiberius/translation.tsv').exists())
                    self.assertNotEqual(subprocess.run(command,capture_output=True).returncode,0)
                    rejected=subprocess.run(command+['--include-incomplete'],capture_output=True,text=True)
                    self.assertIn('unrecognized arguments: --include-incomplete',rejected.stderr)


if __name__=='__main__':unittest.main()

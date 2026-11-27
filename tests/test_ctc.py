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
                table=ctc.codon_table(code)
                summary=ctc.translate_cds(source,root/'proteins.fa',root/'check.tsv',table,{'t1','t2','t3','t4','t5','t6'})
                proteins=dict(ctc.fasta(root/'proteins.fa'))
                self.assertEqual(proteins['g1|t1|chr:0-9(+)'],expected)
                self.assertEqual(proteins['t2'],'M*X')
                self.assertEqual(proteins['t3'],'MX')
                self.assertEqual(proteins['t4'],'MX')
                self.assertEqual(proteins['t5'],'A')
                self.assertEqual(proteins['t6'],'M*')
                self.assertEqual(summary['internal_stops'],2)
                self.assertEqual(summary['partial_codons'],2)
            custom=ctc.codon_table();custom['GCT']='S'
            (root/'custom.json').write_text(json.dumps(custom))
            self.assertEqual(ctc.codon_table(root/'custom.json'),custom)
            source.write_text('>ambiguous\nGCN\n')
            ctc.translate_cds(source,root/'proteins.fa',root/'check.tsv',custom)
            self.assertEqual(dict(ctc.fasta(root/'proteins.fa'))['ambiguous'],'X')
            (root/'bad.json').write_text('{"ATG":"M"}')
            with self.assertRaises(ValueError):ctc.codon_table(root/'bad.json')

    def test_complete_selection_native_ids_and_optional_incomplete(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'combined';out.mkdir()
            t=native(root,'tiberius',model('b','old1','tx1')+model('a','old2','tx2',30),
                [('old1|tx1|b:0-9(+) description','ATGTAATAA','M*'),('old2|tx2|a:29-38(+)','ATGNNNTAA','MX')])
            c=native(root,'carbon_a',model('a','gene1','mrna1',1,'status=complete;confidence=NA')+
                model('a','gene2','mrna2',12,'status=incomplete;confidence=1')+
                model('a','gene3','mrna3',50,'status=complete;confidence=0.001'),
                [('gene1 status=complete','ATGCTGTAA','MS'),('gene2 status=incomplete','ATGCTG','MS'),('gene3','ATGCTGTAA','MS')])
            before={p:p.read_bytes() for d in [t,c] for p in d.iterdir()}
            stats=ctc.combine_models([('tiberius',t),('carbon_a',c)],out)
            self.assertEqual(stats['models'],3);self.assertEqual(stats['carbon_a_confidence_excluded'],1);self.assertEqual(stats['carbon_a_incomplete_excluded'],1)
            self.assertFalse((out/'id_map.tsv').exists())
            expected=['old1|tx1|b:0-9(+)','old2|tx2|a:29-38(+)','gene1']
            self.assertEqual(list(ctc.fasta_index(out/'cds.fa')),expected)
            self.assertEqual(list(ctc.fasta_index(out/'proteins.fa')),expected)
            self.assertIn('>old1|tx1|b:0-9(+) description',(out/'proteins.fa').read_text())
            self.assertIn('M*',(out/'proteins.fa').read_text())
            self.assertIn('confidence=NA',(out/'prediction.gff3').read_text())
            self.assertNotIn('status=incomplete',(out/'prediction.gff3').read_text())
            self.assertTrue(all(p.read_bytes()==data for p,data in before.items()))
            stats=ctc.combine_models([('tiberius',t),('carbon_a',c)],out,include_incomplete=True)
            self.assertEqual(stats['models'],4);self.assertEqual(stats['carbon_a_incomplete_excluded'],0)
            self.assertIn('gene2',ctc.fasta_index(out/'cds.fa'))
            self.assertTrue(all(p.read_bytes()==data for p,data in before.items()))

    def test_confidence_thresholds_na_and_invalid_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'combined';out.mkdir()
            scores=['0.0999','0.1','0.1001','NA','NA']
            gff=''.join(model('x',f'g{i}',f't{i}',1+12*i,
                f"status={'incomplete' if i==4 else 'complete'};confidence={score}")
                for i,score in enumerate(scores))
            c=native(root,'carbon_a',gff,[(f'g{i}','ATGCTGTAA','ML') for i in range(5)])
            before={p:p.read_bytes() for p in c.iterdir()}
            for threshold,partial,expected in [(.1,False,{'g1','g2','g3'}),
                    (.1,True,{'g1','g2','g3','g4'}),(0,False,{'g0','g1','g2','g3'}),
                    (1,True,{'g3','g4'})]:
                stats=ctc.combine_models([('carbon_a',c)],out,partial,threshold)
                self.assertEqual(set(ctc.fasta_index(out/'cds.fa')),expected)
                self.assertEqual(set(ctc.fasta_index(out/'proteins.fa')),expected)
                rows,_=ctc.read_gff(out/'prediction.gff3')
                self.assertEqual({r['attrs']['ID'] for r in rows if r['fields'][2]=='gene'},expected)
                self.assertEqual(stats['models']+stats['carbon_a_incomplete_excluded']+
                    stats['carbon_a_confidence_excluded'],5)
            self.assertTrue(all(p.read_bytes()==data for p,data in before.items()))
            for bad in [-.1,1.1,float('nan'),float('inf')]:
                with self.assertRaisesRegex(ValueError,'confidence-threshold'):
                    ctc.combine_models([('carbon_a',c)],out,confidence_threshold=bad)
            for bad in ['nan','inf','-0.1','1.1','bad',None]:
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
            self.assertEqual((out/'prediction.gff3').read_text(),'##gff-version 3\n'+gff)
            self.assertEqual((out/'cds.fa').read_bytes(),(t/'cds.fa').read_bytes())
            rows,features=ctc.read_gff(out/'prediction.gff3')
            self.assertEqual(len(features['shared']),2)
            self.assertEqual(features['exon'][0]['attrs']['Parent'],'a,b')
            self.assertEqual(features['exon'][0]['attrs']['Derives_from'],'shared')

    def test_native_collisions_are_not_silently_renamed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'combined';out.mkdir()
            t=native(root,'tiberius',model('x','g','tx'),[('g|tx|x:0-9(+)','ATGCTGTAA','ML')])
            c=native(root,'carbon_a',model('y','g','mrna',1,'status=complete;confidence=1'),[('g','ATGCTGTAA','ML')])
            with self.assertRaisesRegex(ValueError,'collide'):
                ctc.combine_models([('tiberius',t),('carbon_a',c)],out)

    def test_invalid_references_and_duplicate_features_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'test.gff3'
            path.write_text(feature('x','CDS',1,9,'ID=c;Parent=absent','+','0'))
            with self.assertRaisesRegex(ValueError,'Dangling'):ctc.read_gff(path)
            path.write_text(feature('x','gene',1,9,'ID=g')+feature('x','mRNA',1,9,'ID=g'))
            with self.assertRaisesRegex(ValueError,'Duplicate'):ctc.read_gff(path)

    def test_pipeline_empty_subsets_custom_table_and_failure(self):
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
                        self.assertEqual((out/'tiberius/proteins.fa').read_text(),'>gT|tT|chr:0-9(+) description\nMS\n')
                    self.assertNotEqual(subprocess.run(command,capture_output=True).returncode,0)


if __name__=='__main__':unittest.main()

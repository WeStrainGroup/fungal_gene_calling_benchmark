"""Test routing, GFF integrity and orchestration without a GPU."""
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

class CTCTests(unittest.TestCase):
    def test_boundaries_and_duplicate_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for tool in ['tiberius','carbon']:(root/tool).mkdir()
            f = root/'input.fa'
            f.write_text(''.join(f'>r{i} description\n'+('A'*n)+'\n' for i,n in enumerate([999,1000,10000,10001])))
            counts = ctc.split_fasta(f,root,1000,10000)
            self.assertEqual(counts,{'tiberius':{'records':2,'bp':11000},'carbon':{'records':1,'bp':10001},'skipped':{'records':1,'bp':999}})
            self.assertIn('description', (root/'tiberius/input.fna').read_text())
            f.write_text('>dup one\nA\n>dup two\nA\n')
            with self.assertRaisesRegex(ValueError,'Duplicate'):ctc.split_fasta(f,root,1,10)

    def test_merge_keeps_phase_and_attributes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);gff=root/'native.gff'
            gff.write_text('##gff-version 3\nchr\ttool\tCDS\t4\t20\t0.1\t-\t2\tID=c;Parent=t1,t2;Derives_from=g;confidence=NA;partial=true\n##FASTA\n>chr\nA\n')
            out=root/'merged.gff3'
            self.assertEqual(ctc.merge_gffs([('tiberius',gff),('carbon',gff)],out),2)
            rows=[l.split('\t') for l in out.read_text().splitlines() if not l.startswith('#')]
            self.assertEqual(rows[0][:8],['chr','tool','CDS','4','20','0.1','-','2'])
            self.assertIn('Parent=tiberius_t1,tiberius_t2',rows[0][8])
            self.assertIn('confidence=NA;partial=true',rows[1][8])

    def test_pipeline_empty_partitions_and_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);binpath=root/'env/bin';binpath.mkdir(parents=True)
            mock=binpath/'python'
            mock.write_text('#!'+sys.executable+'\n'+'''import pathlib,sys
args=sys.argv
if '--input' in args:
 import json
 table=json.loads(pathlib.Path(args[args.index('--codon-table')+1]).read_text())
 assert len(table)==64 and table['TAA']=='*'
 assert args[args.index('--min-record-length')+1]=='0'
 assert args[args.index('--gpu-count')+1]=='1'
if '--fail' in args: sys.exit(7)
f=pathlib.Path(args[args.index('--genome' if '--genome' in args else '--input')+1])
p=pathlib.Path(args[args.index('--out')+1]) if '--out' in args else pathlib.Path(args[args.index('--output-dir')+1])/'prediction.gff'
p.parent.mkdir(parents=True,exist_ok=True)
if '--codingseq' in args:
 pathlib.Path(args[args.index('--codingseq')+1]).write_text('>t description\\nATGCTGTAA\\n')
with p.open('w') as h:
 h.write('##gff-version 3\\n')
 for l in f.read_text().splitlines():
  if l.startswith('>'): h.write(l[1:].split()[0]+'\\tmock\\tCDS\\t1\\t99\\t.\\t+\\t0\\tID=c;Parent=t;confidence=NA\\n')
''');mock.chmod(0o755)
            pipeline=root/'carbon/annotation_pipeline';pipeline.mkdir(parents=True)
            (pipeline/'run_annotation.sh').write_text('#!/bin/sh\npython "$@"\n')
            for label,lengths in [('mixed',[1000,10001]),('short',[1000]),('long',[10001]),('skipped',[999]),('failure',[1000])]:
                f=root/(label+'.fa');f.write_text(''.join(f'>r{i}\n'+('A'*n)+'\n' for i,n in enumerate(lengths)))
                out=root/label
                command=[sys.executable,str(ROOT/'scripts/ctc.py'),'--input',str(f),'--output',str(out),'--cutoff-length','10000',
                    '--tiberius-path',str(mock),'--tiberius-model',str(root/'model'),
                    '--carbon-a-path',str(mock),'--carbon-a-model',str(root/'carbon')]
                if label=='mixed':
                    table=ctc.codon_table();table['CTG']='S'
                    custom=root/'custom.json';custom.write_text(json.dumps(table))
                    command+=['--codon-table',str(custom)]
                if label=='failure':command+=['--tiberius-args=--fail']
                result=subprocess.run(command,capture_output=True,text=True)
                report=json.loads((out/'run.json').read_text())
                if label=='failure':
                    self.assertNotEqual(result.returncode,0)
                    self.assertEqual(report['status'],'failed')
                    self.assertFalse((out/'prediction.gff3').exists())
                else:
                    self.assertEqual(result.returncode,0,result.stderr)
                    if report['routing']['tiberius']['records']:
                        self.assertEqual((out/'tiberius/proteins.fa').read_text(),'>tiberius_t t description\n'+('MS' if label=='mixed' else 'ML')+'\n')
                    actual_table=json.loads((out/'codon_table.json').read_text())
                    self.assertEqual(actual_table['CTG'],'S' if label=='mixed' else 'L')
                    self.assertEqual(report['cds_features'],sum(n>=1000 for n in lengths))
                    self.assertEqual(len(report['components']),len({n<=10000 for n in lengths if n>=1000}))
                    self.assertNotEqual(subprocess.run(command,capture_output=True).returncode,0)

    def test_shared_genetic_codes_and_defects(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'cds.fa'
            source.write_text('>t1 description\nATGCTGTAA\n>t2\nATGTAANNN\n>t3\nATGTAAANN\n>t4\nATGA\n>t5\nGCNTAA\n>t6\nATGTAATAA\n')
            for residue,expected in [('L','ML'),('S','MS'),('A','MA')]:
                table=ctc.codon_table();table['CTG']=residue
                summary=ctc.translate_cds(source,root/'proteins.fa',root/'diagnostics.tsv',table)
                proteins=dict(ctc.fasta(root/'proteins.fa'))
                self.assertEqual(proteins['t1 description'],expected)
                self.assertEqual(proteins['t2'],'M*X')
                self.assertEqual(proteins['t3'],'M*X')
                self.assertEqual(proteins['t4'],'MX')
                self.assertEqual(proteins['t5'],'A')
                self.assertEqual(proteins['t6'],'M*')
                self.assertEqual(summary['internal_stops'],3)
                self.assertEqual(summary['partial_codons'],1)
            custom=ctc.codon_table();custom['GCT']='S'
            source.write_text('>ambiguous\nGCN\n')
            ctc.translate_cds(source,root/'proteins.fa',root/'diagnostics.tsv',custom)
            self.assertEqual(dict(ctc.fasta(root/'proteins.fa'))['ambiguous'],'X')
            (root/'bad.json').write_text('{"ATG":"M"}')
            with self.assertRaises(ValueError):ctc.codon_table(root/'bad.json')

    def test_native_composite_headers_match_merged_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'cds.fa'
            source.write_text('>g1|g1_t1|chr:0-9(+)\nATGCTGTAA\n')
            ctc.translate_cds(source,root/'proteins.fa',root/'diagnostics.tsv',ctc.codon_table(),{'g1_t1'})
            header,seq=next(ctc.fasta(root/'proteins.fa'))
            self.assertEqual(header,'tiberius_g1_t1 g1|g1_t1|chr:0-9(+)')
            self.assertEqual(seq,'ML')
            with self.assertRaisesRegex(ValueError,'Unmatched'):
                ctc.translate_cds(source,root/'proteins.fa',root/'diagnostics.tsv',ctc.codon_table(),{'wrong'})

    def test_prevent_routing_override(self):
        for flag in ['--genome other.fa','--genome=other.fa']:
            with self.assertRaises(ValueError):ctc.extra_args(flag,{'--genome'})

if __name__=='__main__':unittest.main()

"""Exercise portable path configuration and collection without allocating GPUs."""
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

class RunnerTests(unittest.TestCase):
    def test_portable_configuration_and_twenty_cases(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            scripts = root/'scripts';scripts.mkdir()
            shutil.copyfile(ROOT/'scripts/run_efficiency.py', scripts/'run_efficiency.py')
            (scripts/'check_efficiency_resources.py').write_text(
                "import json,sys;from pathlib import Path\n"
                "p=Path(sys.argv[sys.argv.index('--output')+1]);p.mkdir(exist_ok=True)\n"
                "(p/'hardware.json').write_text(json.dumps({'cpu_stat_paths':[]}))\n")
            (scripts/'run_predictor.sh').write_text(
                '#!/bin/sh\nmkdir -p "$3"\n'
                'if [ "$1" = carbon ]; then ext=gff; else ext=gff3; fi\n'
                'printf "##gff-version 3\\n" > "$3/prediction.$ext"\n')
            env = root/'env/bin';env.mkdir(parents=True)
            (env/'python').symlink_to(sys.executable)
            timer = root/'gnu_time'
            timer.write_text('#!/bin/sh\nlog=$4\nshift 4\n"$@"\nrc=$?\nprintf "0.01\\t0.00\\t1\\n" > "$log"\nexit "$rc"\n')
            timer.chmod(0o755)
            config = {'environments':{t:'env' for t in ['carbon','annevo','tiberius','helixer']},
                'carbon_dir':'carbon','annevo_dir':'annevo','tiberius_model':'tiberius',
                'helixer_model':'helixer.h5','gnu_time':'gnu_time'}
            metadata = {t:{'packages':{}} for t in config['environments']}
            metadata['models'] = {}
            for tool,name in [('carbon_a','carbon/model.safetensors'),('annevo','annevo/saved_model/ANNEVO_Fungi.pt'),
                              ('tiberius','tiberius/weights.h5'),('helixer','helixer.h5')]:
                f=root/name;f.parent.mkdir(parents=True,exist_ok=True);f.write_bytes(b'model fixture')
                metadata['models'][tool]={'sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
            (root/'config.json').write_text(json.dumps(config));(root/'software.json').write_text(json.dumps(metadata))
            inputs=root/'inputs/yeast';inputs.mkdir(parents=True);(inputs/'codon_table.json').write_text('{}')
            cases=[]
            for group in ['genome','1-2kb','2-5kb','5-10kb','10-100kb']:
                f=inputs/f'{group}.fna';f.write_text('>contig\nATG\n')
                cases.append(dict(species='yeast',group=group,records=1,bp=3,fasta_sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
            with (root/'cases.tsv').open('w') as h:
                w=csv.DictWriter(h,fieldnames=list(cases[0]),delimiter='\t');w.writeheader();w.writerows(cases)
            command=[sys.executable,str(scripts/'run_efficiency.py'),'--platform','ais','--config',str(root/'config.json'),
                     '--inputs',str(root/'inputs'),'--metadata',str(root/'software.json'),'--cases',str(root/'cases.tsv'),
                     '--output',str(root/'out'),'--species-index','0']
            result=subprocess.run(command,capture_output=True)
            self.assertEqual(result.returncode,0,result.stderr.decode())
            with (root/'out/yeast/metrics.tsv').open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
            self.assertEqual(len(rows),20)
            self.assertTrue(all(r['passed']=='True' and r['cds_features']=='0' for r in rows))
            self.assertNotEqual(subprocess.run(command,capture_output=True).returncode,0)

    def test_collector_checks_native_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            with (ROOT/'results/efficiency/round1_metrics.tsv').open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
            for species in {r['species'] for r in rows}:
                dest=root/species;dest.mkdir();(dest/'hardware.json').write_text('{"passed":true}')
                selected=[r for r in rows if r['species']==species]
                with (dest/'metrics.tsv').open('w') as h:
                    w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(selected)
                for r in selected:
                    case=dest/r['tool']/r['group'];case.mkdir(parents=True)
                    (case/('prediction.gff' if r['tool']=='carbon' else 'prediction.gff3')).write_text('##gff-version 3\n')
                    (case/'time.tsv').write_text('1\t0\t100\n')
            command=[sys.executable,str(ROOT/'scripts/collect_efficiency.py'),'--results',str(root),'--output',str(root/'combined.tsv')]
            self.assertEqual(subprocess.run(command,capture_output=True).returncode,0)
            next(root.glob('*/*/*/prediction.gff')).unlink()
            self.assertNotEqual(subprocess.run(command,capture_output=True).returncode,0)

if __name__ == '__main__':
    unittest.main()

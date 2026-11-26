"""Check replicate medians, species SDs, input identity and failure exclusion."""
import csv
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

class EfficiencyTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        with (ROOT/'results/efficiency/round1_metrics.tsv').open() as handle:
            self.rows = list(csv.DictReader(handle, delimiter='\t'))

    def save(self, index, factor=1, edit=None):
        path = self.root/f'round{index}.tsv'
        rows = [dict(r, elapsed_seconds=str(float(r['elapsed_seconds'])*factor),
                     job_id=r['job_id']+f'_test{index}') for r in self.rows]
        if edit:
            edit(rows)
        with path.open('w') as handle:
            w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
            w.writeheader()
            w.writerows(rows)
        return path

    def run_summary(self, paths):
        return subprocess.run([sys.executable, str(ROOT/'scripts/summarize_efficiency.py'),
            '--metrics', *map(str, paths), '--output', str(self.root/'out')], capture_output=True)

    def test_median_then_species_mean_and_sd(self):
        result = self.run_summary([self.save(1), self.save(2, .5), self.save(3, 2)])
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        with (self.root/'out/summary.tsv').open() as handle:
            summary = list(csv.DictReader(handle, delimiter='\t'))
        for row in summary:
            values = [int(r['bp'])*60/float(r['elapsed_seconds'])/1e6 for r in self.rows
                      if r['tool'] == row['tool'] and r['group'] == row['group']]
            avg = sum(values)/17
            sd = math.sqrt(sum((x-avg)**2 for x in values)/16)
            self.assertAlmostEqual(float(row['mean_mbp_per_minute']), avg)
            self.assertAlmostEqual(float(row['sd_mbp_per_minute']), sd)
        with (self.root/'out/assembly_time.tsv').open() as handle:
            assembly = list(csv.DictReader(handle, delimiter='\t'))
        fastest = min(float(r['mean_minutes']) for r in assembly)
        for r in assembly:
            self.assertAlmostEqual(float(r['relative_time']), float(r['mean_minutes'])/fastest)

    def test_missing_and_reused_rounds_rejected(self):
        first, second = self.save(1), self.save(2)
        self.assertNotEqual(self.run_summary([first, second]).returncode, 0)
        self.assertNotEqual(self.run_summary([first, first, second]).returncode, 0)

    def test_invalid_cases_rejected(self):
        edits = [lambda rows: rows.pop(),
                 lambda rows: rows[0].update(passed='False', exit_code='1'),
                 lambda rows: rows[0].update(bp='1'),
                 lambda rows: rows[0].update(elapsed_seconds='nan'),
                 lambda rows: rows[0].update(elapsed_seconds='0')]
        for edit in edits:
            with self.subTest(edit=edit):
                self.assertNotEqual(self.run_summary([self.save(1), self.save(2), self.save(3, edit=edit)]).returncode, 0)

if __name__ == '__main__':
    unittest.main()

"""Species weighting, sample SD, missing values and duplicate-species safeguards."""
import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
import summarize_all_species as s


def row(species, value, cohort='a'):
    return dict(cohort=cohort, species=species, group='genome', tool='carbon',
                stratum='complete', level='gene', precision=value, recall=value, f1=value)


class SummaryTests(unittest.TestCase):
    def test_shared_species_once(self):
        rows = s.unique_rows([row('yeast', '10'), row('yeast', '10', 'b'), row('other', '30')], s.KEY)
        result = s.mean_sd(rows, ['yeast', 'other'])[0]
        self.assertEqual(result['mean'], 20)
        self.assertAlmostEqual(result['sd'], math.sqrt(200))
        self.assertEqual(result['n_species'], 2)

    def test_conflicting_duplicate_fails(self):
        with self.assertRaises(ValueError):
            s.unique_rows([row('yeast', '10'), row('yeast', '11', 'b')], s.KEY)

    def test_missing_is_not_zero_and_one_has_no_sd(self):
        result = s.mean_sd([row('a', 'NA'), row('b', '20')], ['a', 'b'])[0]
        self.assertEqual((result['mean'], result['sd'], result['n_species']), (20, 'NA', 1))

    def test_zero_is_a_defined_value(self):
        result = s.mean_sd([row('a', '0'), row('b', '20')], ['a', 'b'])[0]
        self.assertEqual(result['mean'], 10)
        self.assertAlmostEqual(result['sd'], math.sqrt(200))

    def test_all_missing(self):
        result = s.mean_sd([row('a', 'NA'), row('b', 'NA')], ['a', 'b'])[0]
        self.assertEqual((result['mean'], result['sd'], result['n_species']), ('NA', 'NA', 0))

    def test_incomplete_species_cell_rejected(self):
        with self.assertRaises(AssertionError):
            s.mean_sd([row('a', '10')], ['a', 'b'])

    def test_supplement_adds_36_cases_without_changing_original_grid(self):
        import json
        import tempfile
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            species = [f'species_{i}' for i in range(12)]
            manifest = root/'species.tsv'
            s.write_table(manifest, [dict(species=sp) for sp in species])
            roots = []
            for label, groups in [('original', s.GROUPS), ('supplement', ('2-5kb',))]:
                folder = root/label
                folder.mkdir()
                roots.append(folder)
                cases = 12*len(groups)*3
                (folder/'STATUS.json').write_text(json.dumps(dict(version='0.2.0', min_record_bp=1000,
                    min_cds_bp=200, confidence_threshold=0.1, source_files_unchanged=True,
                    cohort=label, status='COMPLETE', completed_cases=cases, expected_cases=cases)))
                (folder/'independent_validation.json').write_text(json.dumps(dict(status='PASS',
                    cases=cases, source_files_unchanged=True)))
                metrics = [dict(species=sp, group=g, tool=t, stratum=scope, level=level,
                    precision='50', recall='50', f1='50') for sp in species for g in groups
                    for t in s.TOOLS for scope in (s.SCOPES if t=='carbon' else ('complete',))
                    for level in s.LEVELS]
                s.write_table(folder/'metrics.tsv', metrics)
                s.write_table(folder/'diagnostics.tsv', [dict(species=sp, group=g, tool=t, category='complete')
                    for sp in species for g in groups for t in s.TOOLS])
                s.write_table(folder/'inputs.tsv', [dict(species=sp, group=g) for sp in species for g in groups])
            status = s.aggregate(roots, manifest, root/'combined', groups=(*s.GROUPS, '2-5kb'))
            self.assertEqual((status['completed_cases'], status['expected_cases'], status['metric_rows']), (180, 180, 1260))
            self.assertEqual(len(s.read_table(root/'combined/mean_sd.tsv')), 315)
            self.assertEqual(set(status['species_by_condition']), {*s.GROUPS, '2-5kb'})

    def test_interim_requires_all_tools_within_each_condition(self):
        expected = {(sp,g,t,scope,level) for sp in ('a','b') for g in s.GROUPS for t in s.TOOLS
                    for scope in (s.SCOPES if t=='carbon' else ('complete',)) for level in s.LEVELS}
        observed = expected - {('b','1-2kb','tiberius','complete','gene')}
        units = s.paired_units(observed,expected)
        self.assertNotIn(('b','1-2kb'), units)
        self.assertIn(('b','genome'), units)
        self.assertIn(('a','1-2kb'), units)
        self.assertEqual(len(units), 7)


if __name__ == '__main__':
    unittest.main()

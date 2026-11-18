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

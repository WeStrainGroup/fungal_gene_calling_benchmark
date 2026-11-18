"""Scientific boundary cases for the v0.2 evaluator."""
import copy
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
import evaluate_stratified as e


def model(gene='g', blocks=((101, 400),), strand='+', kind='complete', phase=0, **extra):
    return dict(id=gene+'.t1', gene=gene, record='r', strand=strand, blocks=list(blocks),
                phase=phase, phase_consistent=True, partial=False, qualification='PASS',
                kind=kind, confidence_present=True, confidence=None, status='', incomplete_reason='', **extra)


class StratifiedTests(unittest.TestCase):
    def setUp(self):
        self.genome = {'r': 'A'*2000}

    def test_clip_complete(self):
        m = model(blocks=((101, 160), (201, 260)))
        c = e.clip_model(m, 'f', 80, 300)
        self.assertEqual(c['kind'], 'complete')
        self.assertEqual(c['blocks'], [(22, 81), (122, 181)])

    def test_cut_in_intron_positive_and_negative(self):
        for strand, expected in [('+', 'missing_3prime'), ('-', 'missing_5prime')]:
            m = model(blocks=((101, 160), (201, 260)), strand=strand)
            c = e.clip_model(m, 'f', 1, 190)
            self.assertEqual(c['kind'], expected)
            self.assertEqual(c['phase'], 0)

    def test_both_ends_and_phase(self):
        m = model(blocks=((101, 160), (201, 260)), strand='-')
        c = e.clip_model(m, 'f', 130, 240)
        self.assertEqual((c['kind'], c['removed_5prime'], c['removed_3prime'], c['phase']),
                         ('missing_both', 20, 29, 1))

    def test_no_visible_cds(self):
        self.assertIsNone(e.clip_model(model(blocks=((101, 160), (201, 260))), 'f', 170, 190))

    def test_partition_reconstruction(self):
        original = {'r': 'ACGT'*500}
        genome = {'a': original['r'][:150], 'b': original['r'][150:]}
        mapping = [dict(record=k, chrom='r', start=a, end=b,
                        sequence_sha256=hashlib.sha256(genome[k].encode()).hexdigest())
                   for k, a, b in [('a', 1, 150), ('b', 151, 2000)]]
        for strand in ('+', '-'):
            projected = e.project_reference([model(blocks=((101, 160), (201, 260)), strand=strand)], original, genome, mapping)
            self.assertEqual(len(projected), 2)
        mapping[0]['end'] = 149
        with self.assertRaises(AssertionError):
            e.project_reference([model()], original, genome, mapping)

    def test_length_199_200_201(self):
        for n in (199, 200, 201):
            bucket, _ = e.prediction_bucket(model(blocks=((1, n),)), 'carbon', self.genome, {})
            self.assertEqual(bucket, 'cds_below_200bp' if n < 200 else 'complete')

    def test_input_999_1000(self):
        for n in (999, 1000):
            bucket, _ = e.prediction_bucket(model(), 'tiberius', {'r': 'A'*n}, {})
            self.assertEqual(bucket, 'input_below_1kb' if n < 1000 else 'complete')

    def test_confidence_boundary_and_na(self):
        for value, expected in [(0.099, 'low_confidence'), (0.1, 'complete'), (None, 'complete')]:
            m = model(); m['confidence'] = value
            self.assertEqual(e.prediction_bucket(m, 'carbon', self.genome, {})[0], expected)
        m['confidence_present'] = False
        with self.assertRaises(ValueError):
            e.prediction_bucket(m, 'carbon', self.genome, {})

    def test_native_incomplete_does_not_filter(self):
        m = model(); m['status'] = 'incomplete'
        self.assertEqual(e.prediction_bucket(m, 'carbon', self.genome, {})[0], 'complete')

    def test_same_strand_conflict(self):
        units, domains = e.build_domains([model(), model('p', blocks=((350, 700),), kind='missing_5prime')], self.genome)
        self.assertTrue(all(u['kind'] == 'ignore' for u in units))
        self.assertTrue(domains[('r', '+')]['ignore'].contains((101, 700)))

    def test_opposite_strands_not_conflict(self):
        units, _ = e.build_domains([model(), model('p', strand='-', kind='missing_5prime')], self.genome)
        self.assertEqual({u['kind'] for u in units}, {'complete', 'partial'})

    def test_mixed_isoforms_are_ignored(self):
        a, b = model(), model(kind='missing_3prime')
        b['id'] = 'g.t2'
        units, _ = e.build_domains([a, b], self.genome)
        self.assertEqual(units[0]['reason'], 'mixed_isoform_class')

    def test_short_and_unqualified_reference_are_not_background(self):
        for m in [model(blocks=((101, 299),)), model()]:
            if e.length(m) >= 200:
                m['qualification'] = 'missing_stop'
            _, d = e.build_domains([m], self.genome)
            self.assertEqual(e.prediction_bucket(model(), 'carbon', self.genome, d)[0], 'unscored_reference')

    def test_transitive_ignore_conflict(self):
        a = model(blocks=((100, 400),))
        b = model('p', blocks=((350, 650),), kind='missing_3prime')
        c = model('q', blocks=((600, 900),), kind='missing_3prime')
        units, _ = e.build_domains([a, b, c], self.genome)
        self.assertEqual([u['kind'] for u in units], ['ignore']*3)

    def test_short_isoform_extension_is_not_background(self):
        long = model()
        short = model(blocks=((350, 500),)); short['id'] = 'g.t2'
        units, domains = e.build_domains([long, short], self.genome)
        self.assertEqual(units[0]['kind'], 'complete')
        self.assertTrue(domains[('r', '+')]['ignore'].contains((401, 500)))

    def test_metric_grid_requires_every_level(self):
        case = dict(species='s', group='genome', tool='tiberius')
        rows = [dict(**case, stratum='complete', level=level) for level in e.LEVELS]
        e.validate_metric_grid(rows, [case])
        with self.assertRaises(AssertionError):
            e.validate_metric_grid(rows[:-1], [case])
        with self.assertRaises(AssertionError):
            e.validate_metric_grid(rows+[rows[0]], [case])

    def test_scoring_domains_and_background(self):
        units, d = e.build_domains([model(), model('p', blocks=((501, 800),), kind='missing_5prime')], self.genome)
        self.assertEqual(e.prediction_bucket(model('x', blocks=((450, 750),)), 'carbon', self.genome, d)[0], 'cross_scope')
        self.assertEqual(e.prediction_bucket(model('x', blocks=((550, 749),)), 'carbon', self.genome, d)[0], 'partial')
        self.assertEqual(e.prediction_bucket(model('x', blocks=((901, 1200),)), 'carbon', self.genome, d), ('complete', True))

    def test_maximum_matching_all_isoforms(self):
        a = model('a'); b = model('a', blocks=((501, 800),)); b['id'] = 'a.t2'
        c = model('b')
        refs = {e.gene_key(a): [a, b], e.gene_key(c): [c]}
        self.assertEqual(len(e.matching([a, b], refs)), 2)

    def test_duplicate_predictions(self):
        m = model()
        rows, _ = e.metric_rows([m, copy.deepcopy(m)], {e.gene_key(m): [m]})
        gene = next(r for r in rows if r['level'] == 'gene')
        self.assertEqual((gene['tp'], gene['fp'], gene['fn']), (1, 1, 0))

    def test_wrong_position_same_sequence_strand_and_phase(self):
        truth = model()
        refs = {e.gene_key(truth): [truth]}
        shifted = model(blocks=((501, 800),))
        self.assertEqual(e.legacy.sequence(truth, self.genome), e.legacy.sequence(shifted, self.genome))
        for m in [shifted, model(strand='-'), model(phase=1)]:
            self.assertEqual(len(e.matching([m], refs)), 0)
        m = model(); m['phase_consistent'] = False
        self.assertEqual(len(e.matching([m], refs)), 0)

    def test_partial_coordinates_and_frame(self):
        m = model(kind='missing_5prime', phase=2)
        self.assertEqual(len(e.matching([m], {e.gene_key(m): [m]})), 1)
        wrong = dict(m, phase=0)
        self.assertEqual(len(e.matching([wrong], {e.gene_key(m): [m]})), 0)

    def test_empty_stratum_and_missing_prediction(self):
        rows, _ = e.metric_rows([model()], {})
        self.assertTrue(all(r['f1'] == 'NA' for r in rows))
        rows, _ = e.metric_rows([], {e.gene_key(model()): [model()]})
        self.assertTrue(all(r['recall'] == r['f1'] == 0 for r in rows))

    def test_gff_inheritance_stop_and_phase(self):
        text = ('r\tx\tgene\t101\t400\t.\t+\t.\tID=g;confidence=0.1;status=complete\n'
                'r\tx\tmRNA\t101\t400\t.\t+\t.\tID=t;Parent=g\n'
                'r\tx\tCDS\t101\t397\t.\t+\t0\tParent=t\n'
                'r\tx\tstop_codon\t398\t400\t.\t+\t.\tParent=t\n')
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)/'a.gff'; p.write_text(text)
            m = e.read_models(p)[0]
            self.assertEqual((m['gene'], m['confidence'], m['status'], m['blocks']), ('g', 0.1, 'complete', [(101, 400)]))

    def test_real_mixed_fragment_accounting(self):
        refs = [model(), model('p', blocks=((501, 800),), kind='missing_5prime', phase=1)]
        units, domains = e.build_domains(refs, self.genome)
        predictions = [dict(m) for m in refs] + [model('bg', blocks=((901, 1200),)), model('cross', blocks=((450, 750),))]
        results, stats, audit = e.score_case(predictions, 'carbon', self.genome, units, domains)
        gene = {r['stratum']: r for r in results if r['level'] == 'gene'}
        self.assertEqual((gene['complete']['tp'], gene['complete']['fp']), (1, 1))
        self.assertEqual((gene['partial']['tp'], gene['partial']['fp']), (1, 0))
        self.assertEqual(gene['missing_5prime']['tp'], 1)
        self.assertEqual(sum(r['models'] for r in stats if r['category'] not in ('background_fp', 'subtype_unresolved')), 4)
        self.assertEqual(sum(bool(r['matched_reference']) for r in audit), 2)


if __name__ == '__main__':
    unittest.main()

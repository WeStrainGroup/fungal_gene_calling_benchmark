"""Check reconstruction equivalence and ensure deliberate corruption is detected."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from test_ctc import ctc,model,native

spec=importlib.util.spec_from_file_location('equivalence',Path(__file__).resolve().parents[1]/'scripts/check_ctc_equivalence.py')
equivalence=importlib.util.module_from_spec(spec);spec.loader.exec_module(equivalence)


class EquivalenceTests(unittest.TestCase):
    def fixture(self,root):
        out=root/'combined';out.mkdir()
        t=native(root,'tiberius',model('x','gT','tT'),[('gT|tT|x:0-9(+)','ATGCTGTAA','ML')])
        c=native(root,'carbon_a',model('y','gC','tC',1,'status=complete;confidence=0.1')+
            model('y','gP','tP',20,'status=incomplete;confidence=NA'),
            [('gC','ATGGCTTAA','MA'),('gP','ATGCTG','MS')])
        sources=[('tiberius',t),('carbon_a',c)]
        ctc.combine_models(sources,out,code=12)
        return sources,out

    def test_matches_selected_union(self):
        with tempfile.TemporaryDirectory() as tmp:
            sources,out=self.fixture(Path(tmp))
            result=equivalence.verify(sources,out,12)
            self.assertEqual(result,dict(models=2,genes=2,cds_rows=2,passed=True))
            rows,_=ctc.read_gff(out/'prediction.gff3')
            self.assertTrue(all(set(r['attrs']) <= {'ID','Parent'} for r in rows))

    def test_detects_sequence_or_mapping_corruption(self):
        with tempfile.TemporaryDirectory() as tmp:
            sources,out=self.fixture(Path(tmp))
            for name,old,new in [('cds.fa','ATGGCTTAA','ATGGCCTAA'),('proteins.fa','MA','MF'),
                    ('id_map.tsv','carbon_a:gC','carbon_a:gP')]:
                path=out/name;content=path.read_text();path.write_text(content.replace(old,new))
                with self.assertRaises(AssertionError):equivalence.verify(sources,out,12)
                path.write_text(content)

    def test_detects_coordinate_or_parent_corruption(self):
        with tempfile.TemporaryDirectory() as tmp:
            sources,out=self.fixture(Path(tmp));path=out/'prediction.gff3';content=path.read_text()
            for old,new in [('CDS\t1\t9','CDS\t2\t9'),('ID=g2.t1;Parent=g2','ID=g2.t1;Parent=g1')]:
                path.write_text(content.replace(old,new))
                with self.assertRaises(AssertionError):equivalence.verify(sources,out,12)
                path.write_text(content)


if __name__=='__main__':unittest.main()

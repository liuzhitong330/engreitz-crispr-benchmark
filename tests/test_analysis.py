"""Run: python -m unittest discover -s tests -v"""
import json
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from reproduce import overlap_mask, ranked_ap_setup, verify_official, BINS

class Algorithms(unittest.TestCase):
    def test_overlap_coordinate_boundary(self):
        starts=np.array([0,10,20]); ends=np.array([10,20,30])
        np.testing.assert_array_equal(overlap_mask(10,20,starts,ends),[True,True,True])
        np.testing.assert_array_equal(overlap_mask(10,20,starts,ends,True),[False,True,False])

    def test_ap_ties_and_zero_fill(self):
        y=np.array([1,0,1,0,0]); s=np.array([.5,.5,0,0,0.])
        calc=ranked_ap_setup(y,s)
        self.assertAlmostEqual(calc(np.ones(5)),average_precision_score(y,s),places=14)
        self.assertAlmostEqual(calc(np.ones(5)),.45,places=14)

    def test_weighted_cluster_ap(self):
        y=np.array([1,0,1,0]); s=np.array([.5,.5,.1,0]); weights=np.array([2,2,0,1])
        calc=ranked_ap_setup(y,s)
        self.assertAlmostEqual(calc(weights),average_precision_score(np.repeat(y,weights),np.repeat(s,weights)))
        self.assertTrue(np.isnan(calc(np.zeros(4))))

    def test_bins_do_not_overlap(self):
        for x in [0,19999,20000,99999,100000,2499999]:
            self.assertEqual(sum(lo<=x<hi for _,lo,hi in BINS[1:]),1)
        self.assertEqual(sum(lo<=2500000<hi for _,lo,hi in BINS[1:]),0)

class ReleasedEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.df=pd.read_csv(ROOT/'results/cohort.csv')
        cls.qc=json.loads((ROOT/'results/qc.json').read_text())
        cls.metrics=pd.read_csv(ROOT/'results/metrics.csv')

    def test_same_pairs_and_original_code(self):
        self.assertEqual(verify_official(self.df,ROOT/'results')['status'],'PASS')
        self.assertEqual((len(self.df),self.df.regulated.sum(),self.df.gene.nunique()),(3852,101,54))

    def test_missing_and_labels_are_not_conflated(self):
        missing=self.df[~self.df.matched]
        self.assertEqual((len(missing),missing.regulated.sum()),(521,7))
        self.assertTrue(missing.abc_score.eq(0).all())
        self.assertEqual((~self.df.regulated & self.df.significant & (self.df.effect>0)).sum(),31)
        self.assertTrue((self.df.regulated==(self.df.significant&(self.df.effect<0))).all())

    def test_coordinates_and_coverage(self):
        self.assertEqual(self.qc['endpoint_overlap_differences'],0)
        self.assertEqual(self.qc['name_coordinate_disagreements'],3852)
        self.assertEqual(self.qc['prediction_rows_removed'],10606)
        self.assertEqual(self.qc['missing_by_reason'],{'absent before TSS filter':494,'TSS-filtered only':27})
        self.assertTrue(self.df.cell_type.eq('K562').all())

    def test_metrics_both_methods_same_cohort(self):
        for _,m in self.metrics.iterrows():
            sub=self.df if m.cohort=='all' else self.df[self.df.matched]
            _,lo,hi=next(b for b in BINS if b[0]==m['range'])
            sub=sub[(sub.distance_bp>=lo)&(sub.distance_bp<hi)]
            self.assertEqual(len(sub),m.n)
            self.assertAlmostEqual(average_precision_score(sub.regulated,sub.abc_score),m.abc_ap,places=12)
            self.assertAlmostEqual(average_precision_score(sub.regulated,-sub.distance_bp),m.distance_ap,places=12)
            self.assertLessEqual(m.bootstrap_valid,m.bootstrap_draws)
        main=self.metrics[(self.metrics.cohort=='all')&(self.metrics['range']=='all')].iloc[0]
        self.assertAlmostEqual(main.abc_ap,.6016172397,places=9)
        self.assertAlmostEqual(main.distance_ap,.4140346797,places=9)

    def test_examples_recover_source_rows(self):
        examples=pd.read_csv(ROOT/'results/examples.csv')
        self.assertEqual(len(examples),5)
        self.assertTrue(examples.id.isin(self.df.id).all())

if __name__=='__main__':
    unittest.main()

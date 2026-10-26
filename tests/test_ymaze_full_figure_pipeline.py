"""Run with: python -m unittest discover -s . -p 'test_ymaze_full_figure_pipeline.py'"""
import unittest
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
from ymaze_full_figure_pipeline import (load_and_prepare, calculate_stats,
    exact_permutation, exact_sign_flip, bh_fdr, run)

DATA=Path(__file__).resolve().parents[1]/'Ymaze_AnyMaze_all_measures.csv'

class TestPipeline(unittest.TestCase):
    def test_exact_p(self):
        self.assertEqual(exact_sign_flip([1,2,3]),.25)
        self.assertTrue(np.isnan(exact_sign_flip([1])))
        self.assertAlmostEqual(exact_permutation([1,2],[3,4]), 1/3)
    def test_bh(self):
        q=bh_fdr([.01,.1,np.nan,.02])
        self.assertTrue(np.isnan(q[2]))
        self.assertTrue(np.allclose(q[[0,1,3]],[.03,.1,.03]))
    def test_qc_and_groups(self):
        d,q=load_and_prepare(DATA)
        self.assertEqual(len(d),17)
        self.assertEqual(len(q),23)
        self.assertEqual(int((~q.included).sum()),6)
        self.assertEqual(set(d.loc[d.mouse_id.eq('m78'),'condition']),{'baseline'})
        self.assertEqual(d.groupby(['condition','genotype']).size().to_dict(),{
            ('baseline','WT'):7,('baseline','APP/PS1'):1,
            ('ambtemp','WT'):6,('ambtemp','APP/PS1'):3})
    def test_statistics(self):
        d,_=load_and_prepare(DATA)
        desc,comp,pair,ch=calculate_stats(d)
        c=comp[(comp.metric=='novel_entries_per_mobile_min')&(comp.condition=='ambtemp')].iloc[0]
        self.assertAlmostEqual(c.exact_two_sided_p,1/84,places=7)
        self.assertEqual(int(c.n_APP),3)
        base=comp[(comp.metric=='novel_preference_pct')&(comp.condition=='baseline')].iloc[0]
        self.assertTrue(np.isnan(base.exact_two_sided_p))
        p=pair[(pair.genotype=='WT')&(pair.metric=='novel_preference_pct')].iloc[0]
        self.assertEqual(int(p.n_pairs),3)
        self.assertTrue('m76' in p.paired_mouse_ids)
    def test_export_files(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)
            run(DATA,out)
            for p in ['figures/Figure_Ymaze_memory.png',
                      'figures/Figure_Ymaze_exploratory_all.pdf',
                      'prism/column_tables/Prism_novel_preference_pct.csv',
                      'prism/paired_WT/Prism_WT_paired_novel_preference_pct.csv',
                      'statistics/genotype_comparisons_by_condition.csv',
                      'raw_feature_screen/all_anymaze_raw_measure_screen.csv']:
                self.assertTrue((out/p).exists(),p)
            pair=pd.read_csv(out/'prism/paired_WT/Prism_WT_paired_novel_preference_pct.csv')
            self.assertEqual(set(pair.mouse_id),{'m76','m80','m82'})

if __name__=='__main__':unittest.main()
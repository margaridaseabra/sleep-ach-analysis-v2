"""Unit tests for the exploratory feature module; python -m unittest discover -s tests."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
import pandas as pd

SCRIPT=Path(__file__).resolve().parents[1]/'scripts'/'behavior'/'explore_ymaze_features.py'
spec=importlib.util.spec_from_file_location('explore_ymaze_features',SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class TestYmazeFeatureScreen(unittest.TestCase):
    def test_derived_features(self):
        df=pd.DataFrame([dict(mouse_id='m01',genotype='APP/PS1',condition='ambtemp',novel_entries=8,
                            familiar_entries=4,start_entries=8,novel_time_s=48,familiar_time_s=72,
                            mobile_time_s=120,total_distance_m=4,
                            **{'Average duration of visit to the novel arm zone (s)':6,
                               'novel_latency_s':10})])
        x=mod.features(df).iloc[0]
        self.assertAlmostEqual(x['novel_preference_pct'],40.)
        self.assertAlmostEqual(x['novel_entry_preference_pct'],100*8/12)
        self.assertAlmostEqual(x['novel_entries_per_mobile_min'],4.)
        self.assertAlmostEqual(x['novel_entries_per_meter'],2.)
        self.assertAlmostEqual(x['novel_entries_share_all_arms_pct'],40.)
        self.assertAlmostEqual(x['all_entries_per_mobile_min'],10.)
    def test_permutation_complete_separation(self):
        p=mod.exact_label_p([1,2,3],[9,10,11])
        self.assertAlmostEqual(p,2/20)
    def test_bh(self):
        q=mod.benjamini_hochberg([.01,.04,.2])
        np.testing.assert_allclose(q,[.03,.06,.2])
    def test_excluded_input_rejected(self):
        df=pd.DataFrame([dict(mouse_id='m01',genotype='WT',condition='ambtemp',novel_entries=1,
                              familiar_entries=1,start_entries=1,novel_time_s=50,familiar_time_s=50,
                              mobile_time_s=100,total_distance_m=5,included=False)])
        with self.assertRaisesRegex(ValueError,'QC-filtered'):
            mod.features(df)

if __name__=='__main__':unittest.main()
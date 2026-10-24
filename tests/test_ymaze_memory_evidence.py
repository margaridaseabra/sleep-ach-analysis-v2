"""Run with python -m unittest discover -s tests -p test_ymaze_memory_evidence.py"""
from pathlib import Path
import importlib.util
import unittest
import numpy as np
SCRIPT=Path(__file__).resolve().parents[1]/'scripts'/'behavior'/'analyze_ymaze_memory_evidence.py'
spec=importlib.util.spec_from_file_location('maze_memory',SCRIPT)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class EvidenceTests(unittest.TestCase):
    def test_signflip_symmetric(self):
        self.assertAlmostEqual(m.exact_signflip_mean_p([40,60]),1.0)
    def test_permutation_groups(self):
        self.assertAlmostEqual(m.exact_unpaired_permutation_p([1,2],[3,4]),2/6)
    def test_one_observation_no_inference(self):
        self.assertTrue(np.isnan(m.exact_signflip_mean_p([75])))
        self.assertTrue(np.isnan(m.exact_unpaired_permutation_p([70],[45,50,65])))
    def test_all_mean_bootstrap_constant(self):
        self.assertEqual(m.bootstrap_mean_ci([50,50,50]),(50,50))
if __name__=='__main__':unittest.main()
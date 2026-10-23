"""Run with: python -m unittest discover -s tests -p 'test_ymaze_pipeline.py'."""
import importlib.util
import tempfile
import unittest
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/behavior/run_ymaze.py'
spec = importlib.util.spec_from_file_location('run_ymaze', SCRIPT)
ymaze = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ymaze)


class YMazeTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.folder = Path(self.tempdir.name)
        self.rows = pd.DataFrame({
            'mouse_id': ['m1', 'm2', 'm3', 'm4'],
            'condition': ['baseline', 'baseline', 'ambtemp', 'ambtemp'],
            'genotype': ['WT', 'APP/PS1', 'WT', 'APP/PS1'],
            'tracking_pct': [99., 22., 98., 97.],
            'novel_time_s': [80., 10., 25., 100.],
            'familiar_time_s': [20., 90., 75., 100.],
            'novel_entries': [8, 1, 4, 5],
            'familiar_entries': [2, 7, 6, 5],
            'start_entries': [5, 5, 5, 5],
        })
        self.source = self.folder / 'input.csv'
        self.rows.to_csv(self.source, index=False)

    def test_preference_and_discrimination(self):
        d = ymaze.load_behavior(self.source)
        self.assertAlmostEqual(d.loc[0, 'novel_preference_pct'], 80.)
        self.assertAlmostEqual(d.loc[0, 'discrimination_index'], .6)
        self.assertAlmostEqual(d.loc[0, 'novel_entry_preference_pct'], 80.)

    def test_qc_only_tracking_and_explicit_exclusions(self):
        d = ymaze.load_behavior(self.source)
        ex = self.folder / 'exclude.csv'
        ex.write_text('mouse_id,condition,reason\nm3,*,duplicate_pending_review\n')
        qc = ymaze.apply_qc(d, ex)
        self.assertEqual(qc.included.tolist(), [True, False, False, True])
        self.assertIn('low_tracking', qc.loc[1, 'qc_reason'])
        self.assertIn('manual:duplicate_pending_review', qc.loc[2, 'qc_reason'])
        self.assertTrue(qc.loc[3, 'included'])  # 50% is not a QC failure

    def test_zone_time_discrepancies_are_reviewed_not_excluded(self):
        d = ymaze.load_behavior(self.source)
        d['start_time_s'] = [90., 90., 40., 70.]
        d['center_time_s'] = [110., 110., 20., 30.]
        # Totals: 300, 300, 160, 300. Only m3 needs review.
        qc = ymaze.apply_qc(d, exclusions=None)
        self.assertEqual(qc.zone_time_review_needed.tolist(), [False, False, True, False])
        self.assertTrue(qc.loc[2, 'included'])
        self.assertAlmostEqual(qc.loc[2, 'zone_gap_s'], 140.)

    def test_missing_zone_columns_do_not_trigger_exclusion(self):
        d = ymaze.load_behavior(self.source)
        qc = ymaze.apply_qc(d, exclusions=None)
        self.assertFalse(qc.zone_time_review_needed.any())

    def test_full_export_column_mapping_and_uniqueness(self):
        c = self.rows.rename(columns={'tracking_pct': 'Percentage of frames tracked (%)',
                                       'novel_time_s': 'Time in the novel arm zone (s)',
                                       'familiar_time_s': 'Time in the familiar arm zone (s)'})
        c.to_csv(self.source, index=False)
        self.assertEqual(len(ymaze.load_behavior(self.source)), 4)
        c.loc[3, 'mouse_id'] = 'm1'
        c.loc[3, 'condition'] = 'baseline'
        c.to_csv(self.source, index=False)
        with self.assertRaises(ValueError):
            ymaze.load_behavior(self.source)

    def test_exact_permutation(self):
        from numpy import array
        self.assertAlmostEqual(ymaze.exact_permutation_p(array([1, 2]), array([5, 6])), 2 / 6)


    def test_recognition_summary_respects_individual_pairs(self):
        d = ymaze.load_behavior(self.source)
        # m2's QC failure leaves m1 as the only baseline mouse.
        # m3 and m4 are also unpaired, so output must have no artificial pairs.
        q = ymaze.apply_qc(d, None)
        ymaze.summarize_recognition(q.loc[q.included], self.folder)
        paired = pd.read_csv(self.folder / 'within_mouse_condition_changes.csv')
        self.assertEqual(len(paired), 0)
        summaries = pd.read_csv(self.folder / 'recognition_descriptive_summary.csv')
        self.assertEqual(int(summaries.loc[(summaries.condition == 'baseline') &
                                           (summaries.genotype == 'WT'), 'n_mice'].iloc[0]), 1)
        text = (self.folder / 'behavior_findings.md').read_text()
        self.assertIn('fewer than 2 mice', text)

    def test_spindle_join_rejects_duplicate_animal_sessions(self):
        s = pd.DataFrame({'mouse_id': ['m1', 'm1'], 'condition': ['baseline', 'baseline'],
                          'spindle_density': [4, 5]})
        path = self.folder / 'sp.csv'
        s.to_csv(path, index=False)
        with self.assertRaises(ValueError):
            ymaze.link_spindles(ymaze.load_behavior(self.source), path, self.folder)


if __name__ == '__main__':
    unittest.main()
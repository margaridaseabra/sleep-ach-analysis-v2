"""Tests for the mouse/condition alignment guardrail."""
from __future__ import annotations
import importlib.util
from pathlib import Path
import pandas as pd
import pytest

S = (Path(__file__).resolve().parents[1] / 'scripts/behavior/check_sleep_ymaze_overlap.py')
spec = importlib.util.spec_from_file_location('sleep_bridge', S)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def example_spindles(path: Path):
    pd.DataFrame([
        {'mouse_id': 'mouse69', 'genotype':'APP', 'condition':'warming', 'phase':'WARM_ACTIVE', 'transition_type':'NREM_to_REM', 'occupancy_probability_25_15':.1},
        {'mouse_id': 'mouse69', 'genotype':'APP', 'condition':'warming', 'phase':'WARM_ACTIVE', 'transition_type':'NREM_to_WAKE', 'occupancy_probability_25_15':.05},
        {'mouse_id': 'mouse69', 'genotype':'APP', 'condition':'baseline', 'phase':'BASELINE_ZT3_9', 'transition_type':'NREM_to_REM', 'occupancy_probability_25_15':.2},
        {'mouse_id': 'mouse69', 'genotype':'APP', 'condition':'baseline', 'phase':'BASELINE_ZT3_9', 'transition_type':'NREM_to_WAKE', 'occupancy_probability_25_15':.1},
    ]).to_csv(path,index=False)


def test_mouse_id_aliases():
    assert [m.normalize_mouse_id(s) for s in ['m69','Mouse69','M-069','mouse_0069']] == ['m69']*4


def test_phase_mapping_and_difference(tmp_path):
    p=tmp_path/'sp.csv'; example_spindles(p)
    x=m.load_spindle_probability(p)
    assert len(x)==2
    assert x.loc[x.condition=='ambtemp','spindle_prerem_minus_prewake_25_15'].iloc[0] == pytest.approx(.05)


def test_overlap_avoids_cross_mouse_join(tmp_path):
    p=tmp_path/'sp.csv'; example_spindles(p)
    s=m.load_spindle_probability(p)
    b=pd.DataFrame({'mouse_id':['m83'], 'condition':['ambtemp'], 'genotype':['APP/PS1'], 'novel_preference_pct':[60.]})
    audit,matched=m.build_audit(b,s)
    assert len(matched)==0
    assert audit.spindle_match_status.iloc[0]=='no_matching_sleep_measurement'


def test_overlap_with_genotype_guard(tmp_path):
    p=tmp_path/'sp.csv'; example_spindles(p)
    s=m.load_spindle_probability(p)
    b=pd.DataFrame({'mouse_id':['m69'], 'condition':['ambtemp'], 'genotype':['WT'], 'novel_preference_pct':[60.]})
    with pytest.raises(ValueError,match='Genotype mismatch'):
        m.build_audit(b,s)


def test_overlap_matches_correct_condition(tmp_path):
    p=tmp_path/'sp.csv'; example_spindles(p)
    s=m.load_spindle_probability(p)
    b=pd.DataFrame({'mouse_id':['m69'], 'condition':['ambtemp'], 'genotype':['APP/PS1'], 'novel_preference_pct':[60.]})
    _,matched=m.build_audit(b,s)
    assert len(matched)==1
    assert matched.spindle_pre_rem_occupancy_25_15.iloc[0] == pytest.approx(.1)
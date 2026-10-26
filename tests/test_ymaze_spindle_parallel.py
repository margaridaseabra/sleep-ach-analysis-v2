import importlib.util
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
mod_path = ROOT/'scripts/behavior/summarize_ymaze_parallel.py'
spec = importlib.util.spec_from_file_location('ymaze_parallel',mod_path)
mod=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=mod
spec.loader.exec_module(mod)


def synthetic():
    return pd.DataFrame({
        'mouse_id':['m1','m2','m3'],
        'condition':['ambtemp','ambtemp','ambtemp'],
        'genotype':['WT','APP/PS1','WT'],
        'included':[True,True,False],
        'novel_preference_pct':[40,60,50],
        'novel_entry_preference_pct':[60,66.7,50],
        'novel_entries':[3,8,1],
        'familiar_entries':[2,4,1],
        'start_entries':[1,1,1],
        'total_arm_entries':[6,13,3],
        'mobile_time_s':[100,150,60],
        'novel_time_s':[40,80,20],
        'familiar_time_s':[60,53.333333,20],
        'total_distance_m':[4,6,1],
        'Time mobile in the novel arm zone (s)':[30,65,10],
        'Time immobile in the novel arm zone (s)':[10,15,10],
        'Average duration of visit to the novel arm zone (s)':[13.3,10,20],
    })


def test_filters_and_metrics():
    out=mod.make_features(synthetic())
    assert len(out)==2
    wt=out[out.genotype=='WT'].iloc[0]
    app=out[out.genotype=='APP/PS1'].iloc[0]
    assert wt['novel_entries_per_mobile_min']==pytest.approx(1.8)
    assert app['novel_entries_per_mobile_min']==pytest.approx(3.2)
    assert app['novel_share_all_entries_pct']==pytest.approx(8/13*100)
    assert app['novel_immobility_share_pct']==pytest.approx(15/80*100)
    assert len(mod.summarize(out))==len(mod.FEATURES)*2


def test_missing_columns_and_duplicates():
    with pytest.raises(ValueError,match='Missing required columns'):
        mod.make_features(pd.DataFrame({'mouse_id':['m1']}))
    df=synthetic().iloc[[0,0]].copy()
    with pytest.raises(ValueError,match='Duplicate mouse_id'):
        mod.make_features(df)


def test_zero_mobile_time_produces_nan():
    df=synthetic().iloc[:1].copy()
    df['mobile_time_s']=0
    out=mod.make_features(df)
    assert np.isnan(out['novel_entries_per_mobile_min'].iloc[0])
    assert np.isnan(out['all_entries_per_mobile_min'].iloc[0])
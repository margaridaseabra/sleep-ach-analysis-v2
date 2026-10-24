import importlib.util
from pathlib import Path
import sys
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]
mod_path=ROOT/'scripts/behavior/summarize_ymaze_parallel.py'
spec=importlib.util.spec_from_file_location('ymaze_parallel',mod_path)
mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)

def test_real_dataset_exclusions_and_rate():
    f=Path('/mnt/data/ymaze_feature_results/mouse_level_features.csv')
    if not f.exists():return
    d=mod.make_features(pd.read_csv(f))
    assert len(d)==17
    assert not ((d.mouse_id=='m78')&(d.condition=='ambtemp')).any()
    a=d[(d.condition=='ambtemp')&(d.genotype=='APP/PS1')]
    w=d[(d.condition=='ambtemp')&(d.genotype=='WT')]
    assert len(a)==3 and len(w)==6
    assert a.novel_entries_per_mobile_min.min()>w.novel_entries_per_mobile_min.max()
    assert abs(a.novel_entries_per_mobile_min.mean()-4.591069)<.001
    assert abs(w.novel_entries_per_mobile_min.mean()-2.529957)<.001

def test_missing_required_columns():
    try: mod.make_features(pd.DataFrame({'mouse_id':['m1']}))
    except ValueError:return
    assert False, 'Expected missing-column validation'

def test_no_data_fabrication_when_mobile_time_zero():
    f=Path('/mnt/data/ymaze_feature_results/mouse_level_features.csv')
    if not f.exists():return
    df=pd.read_csv(f).iloc[:1].copy()
    df['mobile_time_s']=0
    out=mod.make_features(df)
    assert np.isnan(out['novel_entries_per_mobile_min'].iloc[0])
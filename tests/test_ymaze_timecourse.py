import importlib.util
from pathlib import Path
import sys
import numpy as np
import pandas as pd

source=Path(__file__).resolve().parents[1]/'scripts/behavior/analyze_ymaze_timecourse.py'
spec=importlib.util.spec_from_file_location('timecourse',source)
mod=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=mod
spec.loader.exec_module(mod)


def demo():
  return pd.DataFrame({'Animal ID':['m76']*5+['m67']*5,'Treatment':['baseline']*10,
    'Segment of test':['0-60 seconds','60-120 seconds','120-180 seconds','180-240 seconds','240-300 seconds']*2,
    'Time in the novel arm zone (s)':[25,26,27,28,29]*2,
    'Time in the familiar arm zone (s)':[15,14,13,12,11]*2,
    'Number of entries to the novel arm zone':[3,4,5,4,3]*2,
    'Number of entries to the familiar arm zone':[3,2,3,4,2]*2,
    'Total distance travelled (m)':[1]*10})


def test_parse_segments():
    assert mod.parse_segment('0-60 seconds')==(0.0,60.0)
    assert mod.parse_segment('60–120 s')==(60.0,120.0)
    assert mod.parse_segment('120 to 180')==(120.0,180.0)


def test_filter_exclusions_and_early_cumulative(tmp_path):
    x=tmp_path/'exclude.csv'
    x.write_text('mouse_id,condition,reason\nm67,baseline,user requested\n')
    df=mod.exclude(mod.parse_export(demo()),x)
    assert df.mouse_id.unique().tolist()==['m76']
    bins,cum,group=mod.build_outputs(df)
    first=cum[cum.window_end_s==60].iloc[0]
    last=cum[cum.window_end_s==300].iloc[0]
    assert round(first.novel_preference_pct,4)==62.5
    assert abs(last.novel_preference_pct - (135/200)*100)<1e-8
    assert bins.expected_60s_bin.all()
    assert len(group)==3


def test_missing_bin_is_not_imputed():
    df=mod.parse_export(demo())
    df=df[(df.mouse_id=='m76') & (df.start_s!=120)]
    _,cum,_=mod.build_outputs(df)
    assert cum.loc[cum.window_end_s==60,'complete'].iloc[0]
    assert cum.loc[cum.window_end_s==120,'complete'].iloc[0]
    assert not cum.loc[cum.window_end_s==300,'complete'].iloc[0]


def test_overlapping_zone_times_flagged():
    df=demo().iloc[:5].copy()
    df.loc[0,'Time in the novel arm zone (s)']=50
    df.loc[0,'Time in the familiar arm zone (s)']=40
    b,c,_=mod.build_outputs(mod.parse_export(df))
    assert b.zone_time_overrun.sum()==1
    assert not c.loc[c.window_end_s==60,'complete'].iloc[0]


def test_repeat_rows_raise():
    df=mod.parse_export(demo().iloc[:5])
    try:
        mod.build_outputs(pd.concat([df,df.iloc[:1]],ignore_index=True))
    except ValueError as e:
        assert 'Duplicate' in str(e)
    else:
        assert False


def test_metadata_mapping(tmp_path):
    df=mod.parse_export(demo())
    fn=tmp_path/'metadata.csv'
    pd.DataFrame({'mouse_id':['m76','m67'],'condition':['baseline','baseline'],'genotype':['WT','APP/PS1']}).to_csv(fn,index=False)
    m=mod.with_metadata(df,fn)
    assert m[m.mouse_id=='m67'].genotype.iloc[0]=='APP/PS1'
    assert m[m.mouse_id=='m76'].genotype.iloc[0]=='WT'
    _,cum,summary=mod.build_outputs(m)
    assert sorted(summary.genotype.unique())==['APP/PS1','WT']
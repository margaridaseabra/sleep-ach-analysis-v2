#!/usr/bin/env python3
"""Hypothesis-led, *descriptive* Y-maze summary for discussion beside independent EEG cohorts.

This never merges mice across cohorts or tests a spindle-memory correlation.
The Y-maze feature differences were discovered exploratorily; previously reported
17-feature FDR results must be retained and may not be replaced by a smaller family.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REQUIRED = ['mouse_id','condition','genotype','included','novel_preference_pct',
            'novel_entry_preference_pct','novel_entries','familiar_entries',
            'start_entries','total_arm_entries','mobile_time_s','novel_time_s',
            'familiar_time_s','total_distance_m',
            'Time mobile in the novel arm zone (s)',
            'Time immobile in the novel arm zone (s)']
FEATURES = {
  'Novel time preference (%)': 'novel_preference_pct',
  'Novel entry preference (%)': 'novel_entry_preference_pct',
  'Novel entries/min mobile': 'novel_entries_per_mobile_min',
  'All entries/min mobile': 'all_entries_per_mobile_min',
  'Novel share of all entries (%)': 'novel_share_all_entries_pct',
  'Mean novel visit duration (s)': 'mean_novel_visit_s',
  'Novel immobility share (%)': 'novel_immobility_share_pct',
}

def make_features(src: pd.DataFrame) -> pd.DataFrame:
    missing = sorted(set(REQUIRED)-set(src.columns))
    if missing: raise ValueError(f'Missing required columns: {missing}')
    d=src.copy()
    if d['included'].dtype == object:
        is_included=d['included'].astype(str).str.lower().isin(['true','1','yes'])
    else: is_included=d['included'].astype(bool)
    d=d.loc[is_included].copy()
    if d.duplicated(['mouse_id','condition']).any():
        raise ValueError('Duplicate mouse_id/condition among included records')
    if not d['condition'].isin(['baseline','ambtemp']).all():
        raise ValueError('Unexpected experimental condition')
    if not d['genotype'].isin(['WT','APP/PS1']).all():
        raise ValueError('Unexpected genotype')
    for c in ['mobile_time_s','total_arm_entries','novel_entries','novel_time_s','familiar_time_s']:
        d[c]=pd.to_numeric(d[c],errors='coerce')
    valid_minutes = d.mobile_time_s/60
    d['novel_entries_per_mobile_min'] = np.where(valid_minutes>0,d.novel_entries/valid_minutes,np.nan)
    d['all_entries_per_mobile_min'] = np.where(valid_minutes>0,d.total_arm_entries/valid_minutes,np.nan)
    d['novel_share_all_entries_pct'] = np.where(d.total_arm_entries>0,d.novel_entries/d.total_arm_entries*100,np.nan)
    d['mean_novel_visit_s'] = pd.to_numeric(d['Average duration of visit to the novel arm zone (s)'],errors='coerce') if 'Average duration of visit to the novel arm zone (s)' in d.columns else np.where(d.novel_entries>0,d.novel_time_s/d.novel_entries,np.nan)
    d['novel_immobility_share_pct'] = np.where(d.novel_time_s>0,pd.to_numeric(d['Time immobile in the novel arm zone (s)'],errors='coerce')/d.novel_time_s*100,np.nan)
    return d

def summarize(d: pd.DataFrame) -> pd.DataFrame:
    out=[]
    for (condition, genotype), group in d.groupby(['condition','genotype'],sort=True):
        for label,col in FEATURES.items():
            v=group[col].dropna().to_numpy(dtype=float)
            out.append(dict(condition=condition,genotype=genotype,feature=label,n=len(v),
                            mean=float(v.mean()) if len(v) else np.nan,
                            median=float(np.median(v)) if len(v) else np.nan,
                            min=float(v.min()) if len(v) else np.nan,
                            max=float(v.max()) if len(v) else np.nan))
    return pd.DataFrame(out)

def leave_one_app_out(d:pd.DataFrame) -> pd.DataFrame:
    a=d[d.condition=='ambtemp']
    wt=a[a.genotype=='WT']
    app=a[a.genotype=='APP/PS1']
    rows=[]
    for c in ['novel_entries_per_mobile_min','all_entries_per_mobile_min','mean_novel_visit_s','novel_entry_preference_pct','novel_preference_pct']:
        for mouse in app.mouse_id:
            rest=app[app.mouse_id!=mouse]
            rows.append(dict(feature=c,omitted_mouse=mouse,n_app_remaining=len(rest),
                             app_minus_wt_mean=float(rest[c].mean()-wt[c].mean()),
                             other_app_all_exceed_all_wt=bool(rest[c].min()>wt[c].max())))
    return pd.DataFrame(rows)

def plot(d:pd.DataFrame, output:Path):
    fig,axs=plt.subplots(2,3,figsize=(14.5,8.4),constrained_layout=True)
    panels=[('Novel-arm time preference (%)','novel_preference_pct', True),
            ('Novel entries/min mobile','novel_entries_per_mobile_min',False),
            ('All-arm entries/min mobile','all_entries_per_mobile_min',False),
            ('Novel entry preference (%)','novel_entry_preference_pct',True),
            ('Novel share of all arm entries (%)','novel_share_all_entries_pct',True),
            ('Mean novel visit duration (s)','mean_novel_visit_s',False)]
    markers={'WT':('o','#6b7280'),'APP/PS1':('D','#377eb8')}
    for n,(name,column,has_chance) in enumerate(panels):
        ax=axs.flat[n]
        groups=[('baseline','WT'),('baseline','APP/PS1'),('ambtemp','WT'),('ambtemp','APP/PS1')]
        if n!=0: groups=groups[2:]
        for j,(c,g) in enumerate(groups):
            subset=d[(d.condition==c)&(d.genotype==g)].sort_values('mouse_id')
            if subset.empty:continue
            symbol,color=markers[g]
            vals=subset[column].to_numpy(float)
            xx=np.full(len(vals),j+1,dtype=float)+np.linspace(-0.09,0.09,len(vals))
            ax.scatter(xx,vals,color=color,marker=symbol,s=52,zorder=3,edgecolors='white',linewidth=.6)
            ax.plot([j+.82,j+1.18],[np.mean(vals),np.mean(vals)],color=color,lw=2)
            if n==0 and c=='baseline' and g=='APP/PS1':
                ax.text(j+1,np.mean(vals)+3,'n=1',ha='center',fontsize=8,color=color)
        if has_chance:
            reference = 100/3 if n == 4 else 50
            ax.axhline(reference,color='#9ca3af',ls='--',lw=1)
        labs=['WT\nBaseline','APP\nBaseline','WT\nWarming','APP\nWarming'] if n==0 else ['WT','APP/PS1']
        ax.set_xticks(range(1,len(groups)+1),labs)
        ax.set_xlim(.55,len(groups)+.45)
        ax.set_title(chr(65+n)+'. '+name,loc='left',fontsize=11,fontweight='semibold')
        ax.grid(axis='y',alpha=.18)
        ax.spines[['top','right']].set_visible(False)
        if n==0: ax.set_ylim(0,100)
        if n==3: ax.set_ylim(0,100)
        if n==4: ax.set_ylim(0,75)
        ax.tick_params(labelsize=9)
    fig.suptitle('Y-maze recognition and exploration | QC-filtered test records',fontsize=15,y=1.03)
    fig.text(.5,-.025,'Independent EEG cohort: these plots cannot establish spindle–memory associations. Exploratory results; no FDR-significant features.',ha='center',fontsize=10)
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=210,bbox_inches='tight')
    plt.close(fig)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,required=True,help='Output from run_ymaze.py after all user-approved exclusions')
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    d=make_features(pd.read_csv(args.input))
    args.output.mkdir(parents=True,exist_ok=True)
    d[['mouse_id','condition','genotype']+list(FEATURES.values())].to_csv(args.output/'hypothesis_led_mouse_features.csv',index=False)
    summarize(d).to_csv(args.output/'hypothesis_led_group_summary.csv',index=False)
    leave_one_app_out(d).to_csv(args.output/'app_warming_leave_one_out.csv',index=False)
    plot(d,args.output/'ymaze_memory_vs_exploration.png')
    a=d[d.condition=='ambtemp']
    nwt=sum(a.genotype=='WT');napp=sum(a.genotype=='APP/PS1')
    print(f'Test records {len(d)}, ambient WT={nwt}, APP/PS1={napp}')
    print('Exploratory only. This does not infer any spindle-memory correlation or memory rescue.')
    print(f'Wrote {args.output}')
if __name__=='__main__':main()
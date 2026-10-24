#!/usr/bin/env python3
"""Exploratory Y-maze phenotype screen. Uses QC-filtered test-trial export.

Analyses are mouse-level. Exact two-sided label permutations are restricted to
independent genotype comparison within the ambient-temperature condition.
Condition comparisons use only genuinely paired WT mice, shown descriptively.
Do not equate exploration phenotypes with recognition, impairment or rescue.
"""
from __future__ import annotations
import argparse
from itertools import combinations,product
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Cohesive, interpretable readouts: not all 90+ AnyMaze fields have physiological meaning.
# Positive differences do not imply stronger recognition for measures of exploration.
FEATURES = [
    ('novel_preference_pct', 'Novel time preference (%)', 'memory'),
    ('novel_entry_preference_pct', 'Novel entry preference (%)', 'memory'),
    ('novel_head_time_preference_pct', 'Novel head-time preference (%)', 'memory'),
    ('novel_head_entry_preference_pct', 'Novel head-entry preference (%)', 'memory'),
    ('novel_time_s', 'Time in novel arm (s)', 'arm occupancy'),
    ('novel_entries', 'Novel arm entries', 'arm exploration'),
    ('mean_novel_visit_s', 'Mean novel visit duration (s)', 'arm exploration'),
    ('novel_latency_s', 'Latency to first novel entry (s)', 'arm exploration'),
    ('novel_entries_per_mobile_min', 'Novel entries / mobile min', 'activity normalized'),
    ('novel_entries_per_meter', 'Novel entries / m traveled', 'activity normalized'),
    ('novel_entries_share_all_arms_pct', 'Novel entries / all arm entries (%)', 'selectivity'),
    ('all_entries_per_mobile_min', 'All arm entries / mobile min', 'general exploration'),
    ('familiar_entries_per_mobile_min', 'Familiar entries / mobile min', 'general exploration'),
    ('total_arm_entries', 'All arm entries', 'general exploration'),
    ('total_distance_m', 'Total distance traveled (m)', 'locomotion'),
    ('mobile_time_s', 'Mobile time (s)', 'locomotion'),
    ('Number of line crossings', 'Line crossings', 'general exploration'),
]

def features(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    cols=['mouse_id','genotype','condition','novel_entries','familiar_entries','start_entries',
          'novel_time_s','familiar_time_s','mobile_time_s','total_distance_m']
    missing=[x for x in cols if x not in df]
    if missing: raise ValueError(f'Y-maze included CSV missing {missing}')
    df['mouse_id']=df.mouse_id.astype(str).str.strip().str.lower()
    df['condition']=df.condition.astype(str).str.strip().str.lower()
    df['genotype']=df.genotype.astype(str).str.strip().str.upper().replace({'APP':'APP/PS1'})
    if 'included' in df.columns:
        # Explicitly reject mixed included and excluded exports, rather than silently using excluded rows.
        if not df['included'].astype(str).str.lower().eq('true').all():
            raise ValueError('Input must be QC-filtered ymaze_included.csv')
    if df.duplicated(['mouse_id','condition']).any():raise ValueError('Duplicate test trials for mouse and condition')
    for k in set(x for x,_,_ in FEATURES if x in df.columns) | set(cols[3:]):
        df[k]=pd.to_numeric(df[k],errors='coerce')
    div=df.novel_time_s+df.familiar_time_s
    df['novel_preference_pct']=100*df.novel_time_s/div.replace(0,np.nan)
    div=df.novel_entries+df.familiar_entries
    df['novel_entry_preference_pct']=100*df.novel_entries/div.replace(0,np.nan)
    nh="Time the animal's head was in the novel arm zone (s)"
    fh="Time the animal's head was in the familiar arm zone (s)"
    ne="Number of entries of the animal's head to the novel arm zone"
    fe="Number of entries of the animal's head to the familiar arm zone"
    if nh in df and fh in df:
        n,f=pd.to_numeric(df[nh],errors='coerce'),pd.to_numeric(df[fh],errors='coerce')
        df['novel_head_time_preference_pct']=100*n/(n+f).replace(0,np.nan)
    if ne in df and fe in df:
        n,f=pd.to_numeric(df[ne],errors='coerce'),pd.to_numeric(df[fe],errors='coerce')
        df['novel_head_entry_preference_pct']=100*n/(n+f).replace(0,np.nan)
    df['mean_novel_visit_s']=pd.to_numeric(df.get('Average duration of visit to the novel arm zone (s)'), errors='coerce')
    df['novel_latency_s']=pd.to_numeric(df.get('novel_latency_s'), errors='coerce')
    df['total_arm_entries']=df.novel_entries+df.familiar_entries+df.start_entries
    df['novel_entries_per_mobile_min']=60*df.novel_entries/df.mobile_time_s.replace(0,np.nan)
    df['familiar_entries_per_mobile_min']=60*df.familiar_entries/df.mobile_time_s.replace(0,np.nan)
    df['all_entries_per_mobile_min']=60*df.total_arm_entries/df.mobile_time_s.replace(0,np.nan)
    df['novel_entries_per_meter']=df.novel_entries/df.total_distance_m.replace(0,np.nan)
    df['novel_entries_share_all_arms_pct']=100*df.novel_entries/df.total_arm_entries.replace(0,np.nan)
    return df

def exact_label_p(wt,app):
    wt=np.asarray(wt,dtype=float);app=np.asarray(app,dtype=float)
    if min(len(wt),len(app))<2:return np.nan
    vals=np.r_[wt,app]; n=len(wt); observed=abs(np.mean(wt)-np.mean(app));count=0;total=0
    for inds in combinations(range(len(vals)),n):
        mask=np.zeros(len(vals),bool);mask[list(inds)]=True
        d=abs(vals[mask].mean()-vals[~mask].mean());count+=d>=observed-1e-12;total+=1
    return count/total

def benjamini_hochberg(pvals):
    p=np.asarray(pvals,dtype=float);q=np.full_like(p,np.nan)
    good=np.isfinite(p);inds=np.flatnonzero(good)
    if len(inds)==0:return q
    ordinds=inds[np.argsort(p[inds])];m=len(ordinds)
    adj=np.minimum.accumulate((p[ordinds]*m/np.arange(1,m+1))[::-1])[::-1]
    q[ordinds]=np.minimum(1,adj);return q

def analyze(df):
    group=df.loc[df.condition=='ambtemp']
    wt=group.loc[group.genotype=='WT'];app=group.loc[group.genotype=='APP/PS1']
    rows=[]
    for key,label,family in FEATURES:
        if key not in df:continue
        a=pd.to_numeric(wt[key],errors='coerce').dropna().values
        b=pd.to_numeric(app[key],errors='coerce').dropna().values
        if min(len(a),len(b))<2:continue
        rows.append(dict(feature=key,label=label,domain=family,n_wt=len(a),n_app=len(b),
                         wt_mean=float(np.mean(a)),app_mean=float(np.mean(b)),
                         mean_diff_app_minus_wt=float(np.mean(b)-np.mean(a)),
                         exact_two_sided_p=exact_label_p(a,b)))
    results=pd.DataFrame(rows)
    if len(results):results['exploratory_BH_q']=benjamini_hochberg(results.exact_two_sided_p.values)
    return results.sort_values('exact_two_sided_p')

def paired_wt(df):
    w=df.loc[df.genotype=='WT']; before=w.loc[w.condition=='baseline'].set_index('mouse_id');after=w.loc[w.condition=='ambtemp'].set_index('mouse_id')
    ids=sorted(set(before.index)&set(after.index));rows=[]
    for mid in ids:
        for key,label,family in FEATURES:
            if key not in df:continue
            a,b=before.loc[mid,key],after.loc[mid,key]
            if np.isfinite(a) and np.isfinite(b):
                rows.append(dict(mouse_id=mid,feature=key,label=label,baseline=float(a),ambtemp=float(b),change=float(b-a)))
    return pd.DataFrame(rows)

def plot(df,path):
    configs=[('novel_entries_per_mobile_min','Novel entries / minute mobile'),
             ('mean_novel_visit_s','Mean novel-arm visit (s)'),
             ('novel_time_s','Total novel-arm time (s)'),
             ('novel_entry_preference_pct','Novel-arm entry preference (%)')]
    wt=df[(df.genotype=='WT')&(df.condition=='ambtemp')];app=df[(df.genotype=='APP/PS1')&(df.condition=='ambtemp')]
    fig,axs=plt.subplots(2,2,figsize=(10,7.5));rng=np.random.default_rng(14)
    for ax,(key,label) in zip(axs.flat,configs):
        for j,group in enumerate([wt,app]):
            vals=group[key].to_numpy(float)
            xx=np.full(len(vals),j)+rng.uniform(-.085,.085,len(vals))
            ax.scatter(xx,vals,s=45,zorder=3)
            mean=float(np.nanmean(vals));ax.plot([j-.13,j+.13],[mean,mean],linewidth=2.2,zorder=4)
            for x,y,mid in zip(xx,vals,group.mouse_id):ax.annotate(mid,(x,y),xytext=(3,4),textcoords='offset points',fontsize=8)
        if key=='novel_entry_preference_pct':ax.axhline(50,color='gray',ls='--',lw=1)
        ax.set_ylabel(label);ax.set_xticks([0,1],['WT (n=6)','APP/PS1 (n=3)']);ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Y-maze · Ambient temperature: exploratory genotype comparison',fontsize=15)
    fig.text(.5,.015,'Points are individual mice; horizontal lines show means. Descriptive findings; multiple testing and small n limit inference.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.035,1,.95]);fig.savefig(path,dpi=200);plt.close(fig)

def main():
    arg=argparse.ArgumentParser();arg.add_argument('--input',type=Path,required=True);arg.add_argument('--output',type=Path,default=Path('data/processed/ymaze/feature_exploration'));a=arg.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    df=features(pd.read_csv(a.input));df.to_csv(a.output/'mouse_level_features.csv',index=False)
    tab=analyze(df);tab.to_csv(a.output/'ambient_genotype_screen.csv',index=False)
    paired_wt(df).to_csv(a.output/'paired_wt_condition_changes.csv',index=False)
    plot(df,a.output/'genotype_feature_panels.png')
    print(f'Records: {len(df)}; APP/PS1 baseline: {sum((df.genotype=="APP/PS1")&(df.condition=="baseline"))}')
    print(tab[['label','n_wt','n_app','wt_mean','app_mean','exact_two_sided_p','exploratory_BH_q']].head(12).to_string(index=False))
    print('WARNING: exploratory screen; q-values control this feature family only, not the earlier unrestricted search over 96 measures.')
if __name__=='__main__':main()
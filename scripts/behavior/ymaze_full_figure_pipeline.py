#!/usr/bin/env python3
"""Reproducible Y-maze memory and exploratory feature figures, stats and Prism CSVs.

Designed for the supplied 300-s ANY-maze novel-arm recognition TEST trial.
The training arm was hidden. These data do not contain a within-session timeline.
Biological unit: mouse; baseline and warming partially repeat WT animals.

Uses the original full-measure CSV or the prior normalized Y-maze dataset.
Reports exact mouse-level permutation p-values for within-condition genotype
comparisons and exact paired sign-flip tests for within-mouse WT changes.
No genotype x condition interaction test is claimed: only one APP/PS1 baseline
record remains after requested exclusions.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import combinations, product
from pathlib import Path
import hashlib
import json
import math
import re
import textwrap

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

# These are ANALYSIS MEASURES, not independent hypotheses. The DI is a linear
# rescaling of novelty-time preference; both are shown because Prism users
# often expect both panels. Counts/occupancy also share information.
@dataclass(frozen=True)
class Measure:
    key: str
    label: str
    unit: str
    family: str
    chance: float | None = None

MEMORY = [
    Measure('novel_preference_pct','Novel-arm time preference','%', 'primary_memory',50),
    Measure('discrimination_index','Novel-arm discrimination index','(N−F)/(N+F)', 'secondary_memory',0),
    Measure('novel_entry_preference_pct','Novel-arm entry preference','%', 'secondary_memory',50),
    Measure('novel_time_s','Time in novel arm','s','secondary_memory'),
    Measure('familiar_time_s','Time in familiar arm','s','secondary_memory'),
    Measure('novel_latency_s','Latency to first novel-arm entry','s','secondary_memory'),
    Measure('novel_entries','Novel-arm entries','count','secondary_memory'),
    Measure('mean_novel_visit_s','Mean novel-arm visit duration','s','secondary_memory'),
    Measure('novel_head_time_preference_pct','Novel-arm head-time preference','%','secondary_memory',50),
]
EXPLORATORY = [
    Measure('novel_entries_per_mobile_min','Novel entries per mobile minute','min⁻¹','exploratory'),
    Measure('all_entries_per_mobile_min','All-arm entries per mobile minute','min⁻¹','exploratory'),
    Measure('novel_entries_per_meter','Novel entries per metre','m⁻¹','exploratory'),
    Measure('novel_entries_share_all_pct','Novel entries / all arm entries','%','exploratory'),
    Measure('total_arm_entries','All-arm entries','count','exploratory'),
    Measure('center_entries','Centre entries','count','exploratory'),
    Measure('line_crossings','Line crossings','count','exploratory'),
    Measure('total_distance_m','Total distance travelled','m','exploratory'),
    Measure('mobile_time_s','Total mobile time','s','exploratory'),
    Measure('average_speed_mps','Average speed','m/s','exploratory'),
    Measure('start_time_s','Time in start arm','s','exploratory'),
    Measure('center_time_s','Time in centre','s','exploratory'),
    Measure('immobile_novel_pct','Novel-arm immobility fraction','%','exploratory'),
    Measure('familiar_speed_mps','Speed in familiar arm','m/s','exploratory'),
    Measure('novel_head_entry_preference_pct','Novel head-entry preference','%','exploratory',50),
    Measure('latency_first_immobile_s','Latency to first immobility','s','exploratory'),
]
ALL = MEMORY + EXPLORATORY

# Exact requested QC decisions. m73 across both conditions; m91 automatically
# poor tracking, but retained here explicitly so the analysis is reproducible.
DEFAULT_EXCLUSIONS = [
    ('m73','*','possible duplicate: near-identical records'),
    ('m91','baseline','very low tracking (22.3%)'),
    ('m67','baseline','unexplained zone-time deficit; user excluded'),
    ('m89','baseline','unexplained zone-time deficit; user excluded'),
    ('m78','ambtemp','unexplained zone-time deficit; user excluded'),
]

RAW = {
    'tracking_pct': 'Percentage of frames tracked (%)',
    'total_distance_m':'Total distance travelled (m)',
    'average_speed_mps':'Average speed (m/s)',
    'mobile_time_s':'Total time mobile (s)',
    'novel_entries':'Number of entries to the novel arm zone',
    'familiar_entries':'Number of entries to the familiar arm zone',
    'start_entries':'Number of entries to the start arm zone',
    'novel_time_s':'Time in the novel arm zone (s)',
    'familiar_time_s':'Time in the familiar arm zone (s)',
    'start_time_s':'Time in the start arm zone (s)',
    'center_time_s':'Time in the center zone (s)',
    'center_entries':'Number of entries to the center zone',
    'line_crossings':'Number of line crossings',
    'novel_latency_s':'Latency to first entry to the novel arm zone (s)',
    'mean_novel_visit_s':'Average duration of visit to the novel arm zone (s)',
    'novel_head_time_s':"Time the animal's head was in the novel arm zone (s)",
    'familiar_head_time_s':"Time the animal's head was in the familiar arm zone (s)",
    'novel_head_entries':"Number of entries of the animal's head to the novel arm zone",
    'familiar_head_entries':"Number of entries of the animal's head to the familiar arm zone",
    'novel_immobile_time_s':'Time immobile in the novel arm zone (s)',
    'familiar_speed_mps':'Average speed in the familiar arm zone (m/s)',
    'latency_first_immobile_s':'Latency to start of first immobile episode (s)',
}

ALIASES = {
    'center_entries':'Number of entries to the center zone',
    'line_crossings':'Number of line crossings',
    'novel_head_time_s':"Time the animal's head was in the novel arm zone (s)",
    'familiar_head_time_s':"Time the animal's head was in the familiar arm zone (s)",
    'novel_head_entries':"Number of entries of the animal's head to the novel arm zone",
    'familiar_head_entries':"Number of entries of the animal's head to the familiar arm zone",
    'novel_immobile_time_s':'Time immobile in the novel arm zone (s)',
    'familiar_speed_mps':'Average speed in the familiar arm zone (m/s)',
    'latency_first_immobile_s':'Latency to start of first immobile episode (s)',
    'average_speed_mps':'Average speed (m/s)',
    'mean_novel_visit_s':'Average duration of visit to the novel arm zone (s)',
}

WT = '#69727e'
APP = '#397aa9'
CONDITIONS = ['baseline','ambtemp']
GENOTYPES = ['WT','APP/PS1']
GROUPS = [('baseline','WT'),('baseline','APP/PS1'),('ambtemp','WT'),('ambtemp','APP/PS1')]


def safe_pct(num:pd.Series, denom:pd.Series, scale:float=100.0) -> pd.Series:
    n=pd.to_numeric(num,errors='coerce')
    d=pd.to_numeric(denom,errors='coerce')
    return scale * n / d.where(d > 0)


def load_and_prepare(path: Path, exclusions:Path|None=None) -> tuple[pd.DataFrame,pd.DataFrame]:
    d=pd.read_csv(path)
    if not {'mouse_id','condition','genotype'}.issubset(d.columns):
        raise ValueError('Input must contain mouse_id, condition and genotype. Use Ymaze_AnyMaze_all_measures.csv or an existing clean table.')
    d['mouse_id']=d['mouse_id'].astype(str).str.strip().str.lower().str.replace(r'^mouse\s*','m',regex=True)
    d['condition']=d['condition'].astype(str).str.lower().str.strip().replace({'warming':'ambtemp','ambient':'ambtemp','ambient temperature':'ambtemp'})
    d['genotype']=d['genotype'].astype(str).str.strip().str.upper().replace({'APP':'APP/PS1','APPPS1':'APP/PS1'})
    if not d['condition'].isin(CONDITIONS).all() or not d['genotype'].isin(GENOTYPES).all():
        raise ValueError('Unknown condition or genotype in input; check labels before analysis.')
    for k,source in RAW.items():
        if k not in d.columns:
            if source in d.columns:d[k]=d[source]
            elif k in ALIASES and ALIASES[k] in d.columns:d[k]=d[ALIASES[k]]
            else: raise ValueError(f'Missing required ANY-maze field: {source} ({k})')
        d[k]=pd.to_numeric(d[k],errors='coerce')
    d['novel_preference_pct']=safe_pct(d.novel_time_s,d.novel_time_s+d.familiar_time_s)
    d['discrimination_index']=safe_pct(d.novel_time_s-d.familiar_time_s,d.novel_time_s+d.familiar_time_s,scale=1)
    d['novel_entry_preference_pct']=safe_pct(d.novel_entries,d.novel_entries+d.familiar_entries)
    d['novel_head_time_preference_pct']=safe_pct(d.novel_head_time_s,d.novel_head_time_s+d.familiar_head_time_s)
    d['novel_head_entry_preference_pct']=safe_pct(d.novel_head_entries,d.novel_head_entries+d.familiar_head_entries)
    d['total_arm_entries']=d.novel_entries+d.familiar_entries+d.start_entries
    d['novel_entries_per_mobile_min']=safe_pct(d.novel_entries,d.mobile_time_s,scale=60)
    d['all_entries_per_mobile_min']=safe_pct(d.total_arm_entries,d.mobile_time_s,scale=60)
    d['novel_entries_per_meter']=safe_pct(d.novel_entries,d.total_distance_m,scale=1)
    d['novel_entries_share_all_pct']=safe_pct(d.novel_entries,d.total_arm_entries)
    d['immobile_novel_pct']=safe_pct(d.novel_immobile_time_s,d.novel_time_s)
    d['zone_time_sum_s']=d[['start_time_s','novel_time_s','familiar_time_s','center_time_s']].sum(axis=1,min_count=4)
    d['zone_time_gap_s']=300-d['zone_time_sum_s']
    # Use full record key; never drop all occurrences of a mouse except wildcard exclusion.
    exclusions_table=[]
    for m,c,reason in DEFAULT_EXCLUSIONS:exclusions_table.append({'mouse_id':m,'condition':c,'reason':reason})
    if exclusions is not None:
        e=pd.read_csv(exclusions)
        if not {'mouse_id','condition'}.issubset(e.columns):
            raise ValueError('Exclusions CSV needs mouse_id and condition fields.')
        for _,r in e.iterrows():
            exclusions_table.append({'mouse_id':str(r.mouse_id).lower().replace('mouse','m').strip(), 'condition':str(r.condition).strip().lower(),'reason':str(r.get('reason','manual exclusion'))})
    exc=[]
    for _,r in d.iterrows():
        reasons=[x['reason'] for x in exclusions_table if x['mouse_id']==r.mouse_id and x['condition'] in ('*',r.condition)]
        if np.isfinite(r.tracking_pct) and r.tracking_pct < 90:
            reasons.append('tracking below 90%')
        exc.append('; '.join(dict.fromkeys(reasons)))
    d['exclude_reason']=exc
    qc=d[['mouse_id','condition','genotype','tracking_pct','zone_time_sum_s','zone_time_gap_s','exclude_reason']].copy()
    qc['included']=qc.exclude_reason.eq('')
    d=d[d.exclude_reason.eq('')].copy().reset_index(drop=True)
    if d[['mouse_id','condition']].duplicated().any():
        raise ValueError('Multiple test records with identical mouse ID and condition; review duplicates.')
    if d[ [v.key for v in ALL] ].isna().any().any():
        bad=d[[v.key for v in ALL]].isna().sum()
        raise ValueError('Required metric contains NaN: '+str(bad[bad>0].to_dict()))
    return d,qc


def exact_permutation(x, y) -> float:
    """Two-sided absolute mean difference; exhaustive label permutations."""
    x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float)
    if len(x)<2 or len(y)<2 or not np.isfinite(np.r_[x,y]).all():return np.nan
    pooled=np.r_[x,y];nx=len(x)
    observed=abs(x.mean()-y.mean())
    ntests=0;tail=0
    for indices in combinations(range(len(pooled)),nx):
        subset=pooled[list(indices)]; other=np.delete(pooled,list(indices))
        delta=abs(subset.mean()-other.mean())
        ntests+=1
        tail+=int(delta >= observed-1e-10)
    return tail/ntests


def exact_sign_flip(values) -> float:
    v=np.asarray(values,dtype=float)
    if len(v)<2 or not np.isfinite(v).all():return np.nan
    obs=abs(v.mean());tail=0
    for s in product([-1.,1.],repeat=len(v)):
        tail+=int(abs((v*np.array(s)).mean())>=obs-1e-10)
    return tail / 2**len(v)


def bh_fdr(pvalues):
    p=np.asarray(pvalues,dtype=float)
    q=np.full(p.shape,np.nan)
    mask=np.isfinite(p)
    if not mask.any():return q
    v=p[mask];order=np.argsort(v);ranked=v[order]
    adjusted=np.minimum.accumulate((ranked*len(v)/np.arange(1,len(v)+1))[::-1])[::-1]
    back=np.empty(len(v));back[order]=np.minimum(adjusted,1)
    q[mask]=back
    return q


def bootstrap_difference(x,y,seed=2026,nboot=8000):
    x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float)
    if min(len(x),len(y))<3:return (np.nan,np.nan)
    rng=np.random.default_rng(seed)
    xb=rng.choice(x,size=(nboot,len(x)),replace=True).mean(1)
    yb=rng.choice(y,size=(nboot,len(y)),replace=True).mean(1)
    return tuple(np.quantile(yb-xb,[.025,.975]))


def calculate_stats(d:pd.DataFrame):
    descriptive=[];comparisons=[];paired=[];chance=[]
    for measure in ALL:
        k=measure.key
        for c,g in GROUPS:
            a=d.loc[(d.condition==c)&(d.genotype==g),k].to_numpy(float)
            descriptive.append({'metric':k,'family':measure.family,'condition':c,'genotype':g,'n':len(a),
                'mean':np.mean(a) if len(a) else np.nan, 'sd':np.std(a,ddof=1) if len(a)>1 else np.nan,
                'median':np.median(a) if len(a) else np.nan,'min':np.min(a) if len(a) else np.nan,
                'max':np.max(a) if len(a) else np.nan,'unit':measure.unit})
        for condition in CONDITIONS:
            wt=d.loc[(d.condition==condition)&(d.genotype=='WT'),k].to_numpy(float)
            app=d.loc[(d.condition==condition)&(d.genotype=='APP/PS1'),k].to_numpy(float)
            ci=bootstrap_difference(wt,app,seed=2026+len(k))
            comparisons.append({'metric':k,'family':measure.family,'condition':condition,
                'n_WT':len(wt),'n_APP':len(app), 'mean_WT':wt.mean() if len(wt) else np.nan,
                'mean_APP':app.mean() if len(app) else np.nan,
                'delta_APP_minus_WT':app.mean()-wt.mean() if len(wt) and len(app) else np.nan,
                'bootstrap_diff_CI_low_descriptive':ci[0],'bootstrap_diff_CI_high_descriptive':ci[1],
                'exact_two_sided_p':exact_permutation(wt,app),
                'test_note':'not estimated: APP/PS1 n<2' if len(app)<2 else 'unpaired mouse-level exact permutation; exploratory'})
        for g in GENOTYPES:
            sub=d[d.genotype==g].pivot(index='mouse_id',columns='condition',values=k)
            if set(CONDITIONS).issubset(sub.columns):
                valid=sub[CONDITIONS].dropna()
                diffs=(valid.ambtemp-valid.baseline).to_numpy(float)
            else:valid=pd.DataFrame();diffs=np.array([])
            paired.append({'metric':k,'family':measure.family,'genotype':g,'n_pairs':len(diffs),
                'paired_mouse_ids':'; '.join(valid.index.astype(str).tolist()),
                'mean_change_ambtemp_minus_baseline':diffs.mean() if len(diffs) else np.nan,
                'exact_two_sided_signflip_p':exact_sign_flip(diffs),
                'test_note':'not estimated (<2 paired mice)' if len(diffs)<2 else 'paired exact sign-flip; exploratory'})
        if measure.chance is not None:
            for c,g in GROUPS:
                a=d.loc[(d.condition==c)&(d.genotype==g),k].to_numpy(float)
                chance.append({'metric':k,'condition':c,'genotype':g,'chance_value':measure.chance,
                  'n':len(a),'mean_minus_chance':a.mean()-measure.chance if len(a) else np.nan,
                  'exact_two_sided_signflip_p':exact_sign_flip(a-measure.chance)})
    desc=pd.DataFrame(descriptive); comp=pd.DataFrame(comparisons);pair=pd.DataFrame(paired);ch=pd.DataFrame(chance)
    comp['exploratory_BH_q_within_family']=np.nan
    # Distinct families and conditions: the primary outcome is intentionally
    # not mixed into the exploratory multiplicity correction.
    for c in CONDITIONS:
        for family in ['secondary_memory','exploratory']:
            idx=(comp.condition==c)&(comp.family==family)
            comp.loc[idx,'exploratory_BH_q_within_family']=bh_fdr(comp.loc[idx,'exact_two_sided_p'])
    return desc,comp,pair,ch


def _stable_jitter(mouse_id,condition):
    token=hashlib.md5((str(mouse_id)+':'+condition).encode()).hexdigest()
    return (int(token[:8],16)/0xffffffff-0.5)*0.115


def draw_panel(ax,d,metric:Measure,comp,pair,caption=''):
    k=metric.key
    xx={('baseline','WT'):-.18,('baseline','APP/PS1'):.18,('ambtemp','WT'):.82,('ambtemp','APP/PS1'):1.18}
    vals=d[k].to_numpy(float)
    lo=float(np.nanmin(vals));hi=float(np.nanmax(vals))
    if metric.chance is not None:
        lo=min(lo,metric.chance);hi=max(hi,metric.chance)
    span=max(hi-lo,1e-6)
    margin=max(span*.20,.05 if abs(hi)<2 else .3)
    bottom=lo-margin*.4;top=hi+margin*1.5
    if lo>=0 and metric.chance!=0:bottom=max(0,bottom)
    if metric.chance is not None:
        ax.axhline(metric.chance,ls=(0,(4,3)),c='#aeb6bc',lw=1.1,zorder=0)
    # Show the within-animal WT pairs only (same animals in both conditions).
    paired_wt=d[d.genotype=='WT'].pivot(index='mouse_id',columns='condition',values=k)
    if set(CONDITIONS).issubset(paired_wt.columns):
        for mid,rr in paired_wt.dropna(subset=CONDITIONS).iterrows():
            ax.plot([xx[('baseline','WT')]+_stable_jitter(mid,'baseline'),xx[('ambtemp','WT')]+_stable_jitter(mid,'ambtemp')],
              [rr.baseline,rr.ambtemp],lw=.8,c=WT,alpha=.19,zorder=1)
    for condition,genotype in GROUPS:
        dd=d[(d.condition==condition)&(d.genotype==genotype)]
        if dd.empty:continue
        xcenter=xx[(condition,genotype)]
        color=WT if genotype=='WT' else APP
        px=np.array([xcenter+_stable_jitter(mid,condition) for mid in dd.mouse_id])
        py=dd[k].to_numpy(float)
        ax.scatter(px,py,s=31,facecolor=color,edgecolor='white',linewidth=.6,zorder=4,alpha=.93)
        if len(dd)>1:
            ax.hlines(np.mean(py),xcenter-.10,xcenter+.10,lw=2.6,color='#101828',zorder=5)
        else:
            ax.scatter([xcenter],[py[0]],s=78,facecolors='none',edgecolors='#101828',lw=.9,zorder=5)
    if k=='novel_preference_pct':
        for _,r in d.iterrows():
            if r.mouse_id in {'m78','m79'}:
                x=xx[(r.condition,r.genotype)]+_stable_jitter(r.mouse_id,r.condition)
                ax.annotate(r.mouse_id,(x,r[k]),xytext=(3,4),textcoords='offset points',fontsize=7,color='#46515e')
    pval=comp.loc[(comp.condition=='ambtemp')&(comp.metric==k),'exact_two_sided_p'].iloc[0]
    qval=comp.loc[(comp.condition=='ambtemp')&(comp.metric==k),'exploratory_BH_q_within_family'].iloc[0]
    pstr=f'p={pval:.3f}' if np.isfinite(pval) else 'p=NA'
    if np.isfinite(qval):pstr+=f'  q={qval:.3f}'
    ax.text(.99,.985, 'Warm APP vs WT  '+pstr, transform=ax.transAxes,ha='right',va='top',fontsize=6.8,color='#455468',
            bbox={'facecolor':'white','alpha':.83,'edgecolor':'none','pad':1.5})
    ax.set_xlim(-.43,1.43);ax.set_ylim(bottom,top)
    ax.set_xticks([0,1]);ax.set_xticklabels(['Baseline','Warming'])
    ax.set_ylabel(metric.unit,fontsize=8)
    ax.set_title(textwrap.fill(caption or metric.label, width=35, max_lines=2),loc='left',fontsize=9.3,pad=12,fontweight='bold')
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.grid(axis='y',color='#d8dce1',alpha=.45,lw=.55)
    ax.set_axisbelow(True)


def write_figure(d,metrics,comp,pair,dest,title,subtitle,ncols=3):
    rows=math.ceil(len(metrics)/ncols)
    width=ncols*4.0; height=rows*3.08+1.05
    fig,axes=plt.subplots(rows,ncols,figsize=(width,height))
    axes=np.asarray(axes).reshape(-1)
    for i,m in enumerate(metrics):draw_panel(axes[i],d,m,comp,pair,caption=f'{chr(65+i)}   {m.label}')
    for ax in axes[len(metrics):]:ax.axis('off')
    fig.suptitle(title,fontsize=16,fontweight='bold',x=.045,ha='left',y=.995)
    fig.text(.045,.966,subtitle,fontsize=8.1,color='#5a6574',ha='left',va='top')
    legend=[Line2D([],[],marker='o',ls='none',color=WT,label='WT',markersize=7),
            Line2D([],[],marker='o',ls='none',color=APP,label='APP/PS1',markersize=7),
            Line2D([],[],color='#aeb6bc',ls='--',label='No preference (where relevant)'),
            Line2D([],[],color=WT,lw=1,alpha=.4,label='WT animals tested twice')]
    fig.legend(handles=legend,loc='upper left',bbox_to_anchor=(.045,.946),ncol=4,frameon=False,fontsize=7.4)
    fig.text(.045,.016,'Dots = individual tests; bar = mean for n≥2. Baseline APP/PS1 n=1: descriptive only. '+
             'Uncorrected p and within-family FDR q refer to genotype comparison within warming only. '+
             'Repeated WT mice are joined by faint lines.',ha='left',fontsize=7,color='#657180')
    fig.subplots_adjust(left=.065,right=.975,bottom=.07,top=.845,wspace=.31,hspace=.51)
    fig.savefig(dest.with_suffix('.png'),dpi=235,bbox_inches='tight',facecolor='white')
    fig.savefig(dest.with_suffix('.pdf'),bbox_inches='tight',facecolor='white')
    plt.close(fig)


def _safe_name(k):return re.sub(r'[^a-zA-Z0-9_]+','_',k)


def export_prism(d,outdir):
    outdir.mkdir(parents=True,exist_ok=True)
    d.to_csv(outdir/'Prism_mouse_level_all_metrics.csv',index=False)
    long=d.melt(id_vars=['mouse_id','genotype','condition'],value_vars=[v.key for v in ALL],var_name='metric',value_name='value')
    long.to_csv(outdir/'Prism_long_all_metrics.csv',index=False)
    (outdir/'column_tables').mkdir(exist_ok=True)
    (outdir/'paired_WT').mkdir(exist_ok=True)
    for m in ALL:
        groups=[]
        for c,g in GROUPS:
            dd=d[(d.condition==c)&(d.genotype==g)].sort_values('mouse_id')
            groups.append(pd.Series(dd[m.key].to_numpy(float),name=f'{g.replace("/","")}_{c}'))
        wide=pd.concat(groups,axis=1)
        # Independent replicate column tables intentionally do NOT encode pairs.
        wide.to_csv(outdir/'column_tables'/f'Prism_{_safe_name(m.key)}.csv',index=False)
        wt=d[d.genotype=='WT'].pivot(index='mouse_id',columns='condition',values=m.key)
        if set(CONDITIONS).issubset(wt.columns):
            wt=wt[CONDITIONS].dropna()
            wt.reset_index().to_csv(outdir/'paired_WT'/f'Prism_WT_paired_{m.key}.csv',index=False)
    # Single table per condition/genotype for plots preserving mouse IDs and replicates.
    d[['mouse_id','genotype','condition']+[x.key for x in ALL]].to_csv(outdir/'Prism_mouse_level_selected.csv',index=False)


def screen_all_raw_measures(d:pd.DataFrame,output:Path,comp:pd.DataFrame,pair:pd.DataFrame):
    """Independent broad post-hoc AnyMaze raw-field screen (reported separately)."""
    output.mkdir(parents=True,exist_ok=True)
    derived={m.key for m in ALL} | {'zone_time_sum_s','zone_time_gap_s','novel_head_time_s',
        'familiar_head_time_s','novel_head_entries','familiar_head_entries','novel_immobile_time_s',
        'exclude_reason','included'}
    rawkeys=[k for k in d.select_dtypes(include=[np.number]).columns if k not in derived]
    rows=[]
    for k in rawkeys:
        group={}
        for c,g in GROUPS:
            group[(c,g)]=d.loc[(d.condition==c)&(d.genotype==g),k].dropna().to_numpy(float)
        a,b=group[('ambtemp','WT')],group[('ambtemp','APP/PS1')]
        if not len(a) or not len(b):continue
        wt=d[d.genotype=='WT'].pivot(index='mouse_id',columns='condition',values=k)
        pairing=wt[CONDITIONS].dropna() if set(CONDITIONS).issubset(wt.columns) else pd.DataFrame()
        vals=pairing['ambtemp']-pairing['baseline'] if len(pairing) else np.array([])
        rows.append({'raw_feature':k,'mean_WT_baseline':np.mean(group[('baseline','WT')]),
          'mean_APP_baseline':np.mean(group[('baseline','APP/PS1')]),
          'mean_WT_warming':np.mean(a),'mean_APP_warming':np.mean(b),
          'delta_APP_minus_WT_warming':np.mean(b)-np.mean(a),
          'exact_p_warming_genotype':exact_permutation(a,b),
          'WT_paired_n':len(vals),'WT_paired_exact_p':exact_sign_flip(vals)})
    screen=pd.DataFrame(rows)
    # No meaningful p for a fully constant or degenerate metric; prevent
    # all-equal values from appearing in the top of the signal list.
    for idx,row in screen.iterrows():
        if np.nanstd(d[row.raw_feature].to_numpy(float))<1e-12:
            screen.loc[idx,'exact_p_warming_genotype']=np.nan
    screen['BH_q_ALL_raw_numeric_fields']=bh_fdr(screen.exact_p_warming_genotype)
    screen=screen.sort_values('exact_p_warming_genotype',na_position='last')
    screen.to_csv(output/'all_anymaze_raw_measure_screen.csv',index=False)
    # The figures are deliberately chosen from a post-hoc all-measures scan:
    # these are *not* confirmatory and some features reflect apparatus geometry.
    selected=screen[(screen.exact_p_warming_genotype<=.10)].head(12)
    pd.DataFrame(selected).to_csv(output/'posthoc_raw_p_le_0p10_candidates.csv',index=False)
    if len(selected):
        figures=[]
        for row in selected.itertuples():
            key=row.raw_feature
            label=str(key).replace('Number of ','').replace('animal\'s ','')
            label=label.replace(' zone','').replace(' to the ',' → ').replace(' in the ',' in ')
            label=label[:62] + ('…' if len(label)>62 else '')
            unit = '%' if '(%)' in key else ('s' if '(s)' in key else ('count' if any(word in key.lower() for word in ('number of','entries','exits')) else 'original units'))
            figures.append(Measure(key,label,unit,'raw_posthoc'))
        fake=screen.rename(columns={'raw_feature':'metric','exact_p_warming_genotype':'exact_two_sided_p',
            'BH_q_ALL_raw_numeric_fields':'exploratory_BH_q_within_family'}).assign(condition='ambtemp')
        write_figure(d,figures,fake,pair,output/'Figure_Ymaze_raw_screen_candidates',
          'Additional ANY-maze signals — raw post-hoc screen',
          'Raw AnyMaze features with uncorrected p ≤ 0.10 · some are redundant or geometry-dependent · not confirmatory',ncols=2 if len(figures)<=6 else 3)
    return screen


def build_report(d,qc,comp,pair,ch,desc,out):
    gp=d.groupby(['condition','genotype']).size().reindex(pd.MultiIndex.from_tuples(GROUPS)).fillna(0).astype(int)
    o=[]
    o.append('# Y-maze recognition and exploratory analysis — autogenerated\n')
    o.append('This is a novel-arm recognition TEST trial following familiarization with the novel arm blocked. Full test length: 300 s.')
    o.append('Primary outcome: time in novel / (novel + familiar) × 100. The discrimination index is exactly 2×preference/100 − 1, so it is NOT independent evidence.\n')
    o.append('## QC and animals')
    o.append(f'Input records = {len(qc)}, included = {len(d)}, excluded = {len(qc)-len(d)}.')
    for c,g in GROUPS:o.append(f'- {c} {g}: n={gp.loc[(c,g)]}')
    o.append('\nExcluded: '+', '.join(f'{r.mouse_id}/{r.condition}' for _,r in qc[~qc.included].iterrows())+'.')
    o.append('\n## Interpretation limits')
    o.append('- Baseline APP/PS1 n=1, hence baseline genotype inference and rescue/intervention interactions **cannot be supported**.')
    o.append('- Some WT mice appear in both conditions; use paired-only WT comparisons for formal condition p-values, and show full condition means descriptively.')
    o.append('- No APP/PS1 mice are paired after QC. Comparing APP/PS1 baseline with warming does not establish individual improvement.')
    o.append('- The 50% novelty-preference reference is a theoretical no-preference value; observed WT groups do not clearly exceed it.')
    o.append('- Exploratory family FDR q-values apply only to features in this report and do not undo earlier data-driven searching over many measurements.')
    o.append('- Visit counts, mean visit duration and movement-normalized counts share denominators/related information and are not independent discoveries.')
    o.append('- An association between exploratory switching in this cohort and spindle organization in separate EEG animals cannot be directly tested.')
    o.append('\n## Ambient-temperature genotype comparisons: selected measures')
    view=comp[(comp.condition=='ambtemp') & (comp.metric.isin(['novel_preference_pct','novel_entry_preference_pct','novel_latency_s','novel_entries_per_mobile_min','all_entries_per_mobile_min','mean_novel_visit_s','total_distance_m']))]
    for _,r in view.iterrows():
        q=('q='+f'{r.exploratory_BH_q_within_family:.3f}' if np.isfinite(r.exploratory_BH_q_within_family) else 'q=not applicable')
        o.append(f'- {r.metric}: WT={r.mean_WT:.3f}, APP/PS1={r.mean_APP:.3f}, Δ={r.delta_APP_minus_WT:+.3f}, exact p={r.exact_two_sided_p:.3f}, {q}')
    o.append('\n## Prism instructions')
    o.append('1. `prism/column_tables/Prism_<metric>.csv` has four **unpaired** columns: WT_baseline, APPPS1_baseline, WT_ambtemp, APPPS1_ambtemp. Use them for grouped individual-point figures; the row alignment is arbitrary, not paired.')
    o.append('2. `prism/paired_WT/Prism_WT_paired_<metric>.csv` retains mouse IDs and true baseline↔warming pairs; use a paired-data analysis only here.')
    o.append('3. `prism/Prism_long_all_metrics.csv` contains animal IDs for every value. Group mean tables, exact permutation p-values, FDR q-values and chance tests are in `statistics/`.')
    o.append('4. The four-condition display is descriptive. Do **not** run a standard independent two-way ANOVA: three WT animals have repeat tests and APP/PS1 baseline n=1.')
    (out/'README_results.md').write_text('\n'.join(o)+'\n',encoding='utf-8')


def run(inputfile:Path,output:Path,exclusions:Path|None=None):
    output.mkdir(parents=True,exist_ok=True)
    (output/'figures').mkdir(exist_ok=True)
    (output/'statistics').mkdir(exist_ok=True)
    data,qc=load_and_prepare(inputfile,exclusions)
    data.to_csv(output/'Ymaze_final_included_mouse_level.csv',index=False)
    qc.to_csv(output/'Ymaze_QC_and_exclusions.csv',index=False)
    desc,comp,pair,ch=calculate_stats(data)
    desc.to_csv(output/'statistics'/'group_descriptives.csv',index=False)
    comp.to_csv(output/'statistics'/'genotype_comparisons_by_condition.csv',index=False)
    pair.to_csv(output/'statistics'/'paired_condition_change.csv',index=False)
    ch.to_csv(output/'statistics'/'recognition_vs_chance.csv',index=False)
    # Show every exploratory metric; separate focus is *descriptive, not confirmatory*.
    signals=comp[(comp.condition=='ambtemp')&(comp.family=='exploratory')].sort_values('exact_two_sided_p')
    signals.to_csv(output/'statistics'/'exploratory_ranked_all_metrics.csv',index=False)
    flags=signals[signals.exact_two_sided_p<=0.10]
    flags.to_csv(output/'statistics'/'uncorrected_p_le_0p10_NOT_confirmatory.csv',index=False)
    export_prism(data,output/'prism')
    raw_screen=screen_all_raw_measures(data,output/'raw_feature_screen',comp,pair)
    subtitle=('17 tests after QC  |  Baseline: WT n=7, APP/PS1 n=1  |  Warming: WT n=6, APP/PS1 n=3  |  Mouse-level observations')
    write_figure(data,MEMORY,comp,pair,output/'figures'/'Figure_Ymaze_memory',
        'Y-maze spatial recognition memory',subtitle,ncols=3)
    write_figure(data,EXPLORATORY,comp,pair,output/'figures'/'Figure_Ymaze_exploratory_all',
        'Y-maze exploration and locomotion',subtitle,ncols=4)
    # Shortlist also includes conceptually important controls, not only low p-values.
    highlighted=[m for m in EXPLORATORY if m.key in set(flags.metric)]
    controls=[m for m in EXPLORATORY if m.key in {'total_distance_m','total_arm_entries'}]
    highlights=(highlighted+ [m for m in controls if m not in highlighted])[:9]
    write_figure(data,highlights,comp,pair,output/'figures'/'Figure_Ymaze_exploratory_focus',
        'Exploratory differences and movement controls',
        'Panels with uncorrected warming genotype p ≤ 0.10, plus locomotor controls; NOT FDR-significant; hypothesis-generating only',ncols=3)
    build_report(data,qc,comp,pair,ch,desc,output)
    with open(output/'README_results.md','a',encoding='utf-8') as f:
        f.write('\n## Broad raw-field screen\nAll numeric AnyMaze source fields are ranked in `raw_feature_screen/all_anymaze_raw_measure_screen.csv`. The companion raw-screen figure presents only uncorrected p≤0.10 candidates and may include coordinate/zone geometry artifacts or strongly redundant measures. The raw-screen q-values reflect this broader family and should be considered descriptive.\n')
    print(f'Input test records: {len(qc)}; retained: {len(data)}; excluded: {len(qc)-len(data)}')
    print(data.groupby(['condition','genotype']).size().to_string())
    print('\nAmbient genotype contrasts (sorted by exploratory p):')
    print(signals[['metric','mean_WT','mean_APP','exact_two_sided_p','exploratory_BH_q_within_family']].head(10).round(4).to_string(index=False))
    print(f'\nOutputs: {output}')


def main():
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input',required=True,type=Path,help='Original 23-record AnyMaze CSV or existing normalized full-measure CSV')
    ap.add_argument('--output',required=True,type=Path,help='Output folder for figures, stats and Prism CSV files')
    ap.add_argument('--exclusions',type=Path,default=None,help='Optional additional CSV mouse_id,condition,reason. Built-in user exclusions always apply.')
    args=ap.parse_args()
    run(args.input,args.output,args.exclusions)

if __name__=='__main__':main()
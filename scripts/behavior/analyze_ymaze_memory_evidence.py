#!/usr/bin/env python3
"""Y-maze recognition: chance performance, genotype contrast, exploration controls.

Analyzes *one row per mouse/condition* from run_ymaze.py's ymaze_included.csv.
Never treats distinct mice as paired or labels nonsignificant genotype contrasts as rescue.
"""
from __future__ import annotations

import argparse
import csv
from itertools import combinations, product
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PRIMARY = 'novel_preference_pct'
ENTRY = 'novel_entry_preference_pct'
CONDITIONS = ('baseline', 'ambtemp')
GENOTYPES = ('WT', 'APP/PS1')
FILL = {'WT': '#81858e', 'APP/PS1': '#4385b7'}


def read_rows(path):
    with open(path, newline='', encoding='utf8') as f:
        rows = list(csv.DictReader(f))
    required = {'mouse_id', 'genotype', 'condition', PRIMARY, ENTRY,
                'total_distance_m', 'total_arm_entries'}
    if not rows or not required.issubset(rows[0]):
        raise ValueError('Expected nonempty ymaze_included.csv with columns: '+', '.join(sorted(required)))
    for row in rows:
        row['mouse_id'] = row['mouse_id'].strip().lower()
        row['condition'] = row['condition'].strip().lower()
        row['genotype'] = row['genotype'].strip()
        for key in (PRIMARY, ENTRY, 'total_distance_m', 'total_arm_entries'):
            row[key] = float(row[key]) if row[key].strip() else float('nan')
        if row['condition'] not in CONDITIONS or row['genotype'] not in GENOTYPES:
            raise ValueError(f"Unknown condition/genotype: {row['condition']}/{row['genotype']}")
    if len({(r['mouse_id'],r['condition']) for r in rows}) != len(rows):
        raise ValueError('Multiple test trials for same mouse and condition; resolve first')
    return rows


def exact_signflip_mean_p(values, reference=50.0):
    """Two-sided exact mean sign-flip p, assumes symmetry of centered errors."""
    x = np.asarray(values,dtype=float)-reference
    if len(x)<2: return float('nan')
    observed = abs(x.mean()); n=0; extreme=0
    for signs in product((-1,1),repeat=len(x)):
        n+=1
        extreme+=abs(np.mean(x*np.asarray(signs)))>=observed-1e-12
    return extreme/n


def exact_unpaired_permutation_p(group_a,group_b):
    a=np.asarray(group_a,dtype=float);b=np.asarray(group_b,dtype=float)
    if min(len(a),len(b))<2: return float('nan')
    vals=np.r_[a,b]; observed=abs(b.mean()-a.mean());n=0;extreme=0
    for indices in combinations(range(len(vals)),len(b)):
        mask=np.zeros(len(vals),dtype=bool);mask[list(indices)]=True
        difference=abs(vals[mask].mean()-vals[~mask].mean())
        n+=1;extreme+=difference>=observed-1e-12
    return extreme/n


def bootstrap_mean_ci(values,seed=42,iterations=20000):
    vals=np.asarray(values,dtype=float)
    if len(vals)<2: return (float('nan'),float('nan'))
    rng=np.random.default_rng(seed)
    sampled=rng.choice(vals,size=(iterations,len(vals)),replace=True).mean(axis=1)
    return tuple(np.quantile(sampled,[.025,.975]))


def bootstrap_difference_ci(a,b,seed=42,iterations=20000):
    a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float)
    if min(len(a),len(b))<2: return (float('nan'),float('nan'))
    rng=np.random.default_rng(seed)
    ab=rng.choice(a,size=(iterations,len(a)),replace=True).mean(axis=1)
    bb=rng.choice(b,size=(iterations,len(b)),replace=True).mean(axis=1)
    return tuple(np.quantile(bb-ab,[.025,.975]))


def group(rows,condition,genotype):
    return [r for r in rows if r['condition']==condition and r['genotype']==genotype]


def fmt(x,d=1):
    return f'{x:.{d}f}' if np.isfinite(x) else 'not estimable'


def make_results(rows):
    groups=[]
    for condition in CONDITIONS:
        for genotype in GENOTYPES:
            subset=group(rows,condition,genotype)
            if not subset:continue
            x=np.array([r[PRIMARY] for r in subset],dtype=float)
            ent=np.array([r[ENTRY] for r in subset],dtype=float)
            low,high=bootstrap_mean_ci(x)
            groups.append(dict(condition=condition,genotype=genotype,n=len(x),
                mean_time_preference_pct=float(np.mean(x)),median_time_preference_pct=float(np.median(x)),
                time_above_50_count=int(np.sum(x>50)),
                mean_entry_preference_pct=float(np.mean(ent)),
                exact_signflip_p_vs_50=exact_signflip_mean_p(x),
                bootstrap_95ci_mean_low=low,bootstrap_95ci_mean_high=high,
                mean_distance_m=float(np.mean([r['total_distance_m'] for r in subset])),
                mean_arm_entries=float(np.mean([r['total_arm_entries'] for r in subset]))))
    wt=group(rows,'ambtemp','WT');app=group(rows,'ambtemp','APP/PS1')
    x=[r[PRIMARY] for r in wt];y=[r[PRIMARY] for r in app]
    contrast=dict(n_wt=len(x),n_app=len(y),mean_diff_app_minus_wt=float(np.mean(y)-np.mean(x)),
        p_exact_unpaired=exact_unpaired_permutation_p(x,y),bootstrap_95ci=bootstrap_difference_ci(x,y))
    # Complete pairs across both conditions, NEVER fill missing with cross-mouse comparisons.
    by_id={}
    for r in rows:by_id.setdefault(r['mouse_id'],{})[r['condition']]=r
    pairs=[]
    for mid,v in sorted(by_id.items()):
        if 'baseline' in v and 'ambtemp' in v:
            if v['baseline']['genotype']!=v['ambtemp']['genotype']:
                raise ValueError(f'Genotype mismatch in repeated mouse: {mid}')
            pairs.append(dict(mouse_id=mid,genotype=v['baseline']['genotype'],
                baseline=v['baseline'][PRIMARY],ambtemp=v['ambtemp'][PRIMARY],
                delta_pp=v['ambtemp'][PRIMARY]-v['baseline'][PRIMARY]))
    correlations=[]
    for condition in CONDITIONS:
        for geno in GENOTYPES:
            sub=group(rows,condition,geno)
            if len(sub)>=5:
                res=spearmanr([r[PRIMARY] for r in sub],[r['total_distance_m'] for r in sub])
                correlations.append(dict(condition=condition,genotype=geno,n=len(sub),
                    rho_distance=float(res.statistic),p_uncorrected=float(res.pvalue)))
    return groups,contrast,pairs,correlations


def save_csv(path,rows):
    if not rows: return
    with open(path,'w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def draw_panel(ax,rows,key,title,ylab):
    rng=np.random.default_rng(21)
    for xi,condition in enumerate(CONDITIONS):
        for gi,geno in enumerate(GENOTYPES):
            sub=group(rows,condition,geno);vals=[r[key] for r in sub]
            if not sub:continue
            x=xi+(-.17 if gi==0 else .17)
            jitter=rng.uniform(-.048,.048,size=len(sub))
            ax.scatter(np.array([x]*len(sub))+jitter,vals,s=61,zorder=3,
                color=FILL[geno],edgecolors='white',linewidths=.65)
            ax.plot([x-.08,x+.08],[np.mean(vals)]*2,color='#191a1d',lw=2.2,zorder=4)
            if key==PRIMARY and condition=='baseline' and geno=='WT':
                target=next((v for v in sub if v['mouse_id']=='m79'),None)
                if target:
                    j=sub.index(target)
                    ax.annotate('m79',xy=(x+jitter[j],target[key]),xytext=(14,11),textcoords='offset points',fontsize=8,color='#656565')
    ax.axhline(50,ls='--',color='#aa6a41',lw=1.2,alpha=.8)
    ax.set_xticks([0,1],['Baseline','Ambient warming'])
    ax.set_xlim(-.42,1.42);ax.set_ylim(0,100)
    ax.set_ylabel(ylab);ax.set_title(title,loc='left',fontweight='bold',fontsize=11)
    ax.grid(axis='y',alpha=.14)
    ax.spines[['top','right']].set_visible(False)


def plot_figure(rows,pairs,correlations,path):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9.5,'axes.labelcolor':'#202329','text.color':'#202329'})
    fig,axes=plt.subplots(2,2,figsize=(12.2,9.1))
    draw_panel(axes[0,0],rows,PRIMARY,'A  Novel-arm time preference (primary)','Novel / (novel + familiar), %')
    draw_panel(axes[0,1],rows,ENTRY,'B  Novel-arm entry preference','Novel / (novel + familiar) entries, %')
    ax=axes[1,0]
    shapes={'baseline':'o','ambtemp':'^'}
    for condition in CONDITIONS:
        for geno in GENOTYPES:
            sub=group(rows,condition,geno)
            ax.scatter([r['total_distance_m'] for r in sub],[r[PRIMARY] for r in sub],s=70,
                       marker=shapes[condition],color=FILL[geno],edgecolors='white',linewidth=.7,
                       label=f"{geno} · {'baseline' if condition=='baseline' else 'warming'}")
    ax.axhline(50,ls='--',color='#aa6a41',lw=1.2,alpha=.8)
    ax.set(xlabel='Total distance travelled (m)',ylabel='Novel-arm time preference (%)',ylim=(0,100))
    ax.set_title('C  Exploration and novel-arm preference',loc='left',fontweight='bold',fontsize=11)
    ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.13)
    ax.legend(frameon=False,fontsize=8,loc='upper left',ncol=2)
    ax=axes[1,1]
    wtpairs=[p for p in pairs if p['genotype']=='WT']
    for p in wtpairs:
        ax.plot([0,1],[p['baseline'],p['ambtemp']],marker='o',color='#7d838e',lw=1.65,alpha=.86)
        ax.annotate(p['mouse_id'],xy=(1,p['ambtemp']),xytext=(5,-1),textcoords='offset points',fontsize=8,color='#555')
    ax.axhline(50,ls='--',color='#aa6a41',lw=1.2,alpha=.8)
    ax.set(xlim=(-.15,1.3),ylim=(0,100),ylabel='Novel-arm time preference (%)')
    ax.set_xticks([0,1],['Baseline','Ambient warming'])
    ax.set_title('D  Within-mouse change (WT only)',loc='left',fontweight='bold',fontsize=11)
    ax.grid(axis='y',alpha=.13);ax.spines[['top','right']].set_visible(False)
    import matplotlib.lines as mlines
    handles=[mlines.Line2D([0],[0],marker='o',color='none',markerfacecolor=FILL[g],markeredgecolor='white',markersize=8,label=g) for g in GENOTYPES]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.975),ncol=2,frameon=False)
    fig.suptitle('Novel-arm Y-maze | separate behavioral cohort',x=.06,y=.995,ha='left',fontsize=17,fontweight='bold')
    fig.text(.06,.015,'Each point = one test record. Black ticks = group means. Dashed line = theoretical no-preference (50%).\n'
                     'Exclusions: m73 (both conditions), m91 and m67/m89 (baseline), m78 (warming). These data alone do not demonstrate memory rescue.',fontsize=8,color='#555')
    fig.tight_layout(rect=(0,.05,1,.96),h_pad=2.2,w_pad=2.2)
    fig.savefig(path,dpi=220,facecolor='white');plt.close(fig)


def report(groups,contrast,pairs,corr):
    g={(a['condition'],a['genotype']):a for a in groups}
    rows=['# Y-maze spatial recognition — statistical evidence and limitations','',
    'Dataset: 17 retained test records after six documented record exclusions; independent behavioral cohort.',
    'Novel-arm time preference is primary. Chance = 50% only if familiar and novel arms are equally preferred without learning.',
    'Two-sided exact sign-flip probabilities are exploratory and assume symmetric centered errors; the small n limits interpretation.',
    'Bootstrap intervals are descriptive percentile CIs and unstable for very small samples.','',
    '## Recognition versus chance','',
    '| Condition | Genotype | n | Mean time preference | Median | Above 50% | Exact sign-flip p | Bootstrap 95% CI of mean |','|---|---|---:|---:|---:|---:|---:|---|']
    for condition in CONDITIONS:
        for geno in GENOTYPES:
            a=g.get((condition,geno))
            if a:
                ci_label = (f"{fmt(a['bootstrap_95ci_mean_low'])}% to {fmt(a['bootstrap_95ci_mean_high'])}%" if a['n'] >= 2 else 'not estimable')
                rows.append(f"| {condition} | {geno} | {a['n']} | {fmt(a['mean_time_preference_pct'])}% | {fmt(a['median_time_preference_pct'])}% | {a['time_above_50_count']}/{a['n']} | {fmt(a['exact_signflip_p_vs_50'],3)} | {ci_label} |")
    lo,hi=contrast['bootstrap_95ci'];rows.extend(['','## Genotype contrast under ambient temperature','',
        f"APP/PS1 minus WT mean novelty preference = {fmt(contrast['mean_diff_app_minus_wt'])} percentage points (APP n={contrast['n_app']}, WT n={contrast['n_wt']}).",
        f"Exploratory exact label-permutation p = {fmt(contrast['p_exact_unpaired'],3)}; descriptive bootstrap 95% CI: {fmt(lo)} to {fmt(hi)} percentage points.",
        'This is a genotype comparison *within warming*, not a warming effect, memory rescue, or a genotype × temperature interaction.',
        '', '## Within-animal temperature comparisons','',
        '| Mouse | Genotype | Baseline % | Warming % | Difference (pp) |','|---|---|---:|---:|---:|'])
    for p in pairs:rows.append(f"| {p['mouse_id']} | {p['genotype']} | {fmt(p['baseline'])} | {fmt(p['ambtemp'])} | {fmt(p['delta_pp'])} |")
    rows.extend(['', f"No APP/PS1 animals have both included test records. WT paired n={sum(p['genotype']=='WT' for p in pairs)}, without consistent direction.",
                 '', '## Checks affecting validity','',
                 '- WT novel preference is not convincingly above 50% in either condition, so a reliable task positive control has not been demonstrated.',
                 '- A single APP/PS1 baseline test (m78) cannot establish a baseline deficit; a rescue requires a deficit and an appropriate improvement relative to a comparator.',
                 '- Baseline/ambient-temperature groups overlap only partly. Treat condition means as partly independent, not a fully paired treatment design.',
                 '- Total exploration, novelty seeking, arm-position bias, anxiety and locomotor activity can affect zone occupancy.',
                 '- Time and entry preference diverge for individual mice (notably m79); retain these subjects but inspect time-vs-entry interpretation.',
                 '- Confirm novel-arm assignment was counterbalanced, both accessible arms were explored during training, retention interval, repeated-test order, temperature timing, and animal sex/age.',
                 '- Unlike sleep-correlated studies, these behavioral animals are not mouse-matched to spindle recordings. Do not report individual spindle–memory correlations.',
                 '', '## Descriptive locomotor associations',''])
    for c in corr: rows.append(f"- {c['genotype']} {c['condition']}: Spearman ρ = {fmt(c['rho_distance'],2)} for distance versus novelty preference (n={c['n']}, uncorrected p = {fmt(c['p_uncorrected'],3)}). Exploratory only.")
    rows.extend(['','## Conservative results statement','',
                 'Following exclusion of recordings with tracking or zone-time QC concerns, the remaining behavioral cohort showed heterogeneous novel-arm preference. WT mice did not show convincing preference above the theoretical 50% reference at baseline or under ambient warming. A baseline genotype deficit could not be assessed reliably because only one APP/PS1 record remained. Under warming, APP/PS1 mice had a numerically higher novel-arm time preference than WT mice, but the difference was uncertain and does not establish an improvement or rescue. The data therefore do not support a claim of memory impairment or memory rescue in this cohort.'])
    return '\n'.join(rows)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=Path('data/processed/ymaze/ymaze_included.csv'))
    parser.add_argument('--output',type=Path,default=Path('data/processed/ymaze/memory_evidence'))
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    rows=read_rows(args.input);g,contrast,pairs,corr=make_results(rows)
    save_csv(args.output/'recognition_chance_tests.csv',g)
    save_csv(args.output/'paired_test_records.csv',pairs)
    save_csv(args.output/'locomotion_correlations.csv',corr)
    plot_figure(rows,pairs,corr,args.output/'ymaze_memory_figure.png')
    (args.output/'ymaze_memory_interpretation.md').write_text(report(g,contrast,pairs,corr),encoding='utf8')
    print(f'Processed {len(rows)} records, {len({r["mouse_id"] for r in rows})} unique mice')
    print(f'APP/PS1 minus WT in warming: {contrast["mean_diff_app_minus_wt"]:.1f} pp; exact p={contrast["p_exact_unpaired"]:.3f}')
    print(f'Wrote results to {args.output}')


if __name__=='__main__':main()
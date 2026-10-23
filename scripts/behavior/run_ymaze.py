#!/usr/bin/env python3
"""Reproducible two-trial Y-maze novel-arm recognition analysis.

Input is an AnyMaze *mouse-level test-trial* CSV, NOT the familiarization
trial. Never exclude a mouse solely because its novel preference is low.
"""
from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

COLUMNS = {
    'Percentage of frames tracked (%)': 'tracking_pct',
    'Total distance travelled (m)': 'total_distance_m',
    'Total time mobile (s)': 'mobile_time_s',
    'Number of entries to the novel arm zone': 'novel_entries',
    'Number of entries to the familiar arm zone': 'familiar_entries',
    'Number of entries to the start arm zone': 'start_entries',
    'Time in the novel arm zone (s)': 'novel_time_s',
    'Time in the familiar arm zone (s)': 'familiar_time_s',
    'Time in the start arm zone (s)': 'start_time_s',
    'Time in the center zone (s)': 'center_time_s',
    'Latency to first entry to the novel arm zone (s)': 'novel_latency_s',
}
REQUIRED = ('mouse_id', 'condition', 'genotype', 'tracking_pct',
            'novel_time_s', 'familiar_time_s')
BEHAVIOR = ['novel_preference_pct', 'discrimination_index',
            'novel_entry_preference_pct', 'total_distance_m', 'total_arm_entries']
GENOTYPES = ('WT', 'APP/PS1')


def load_behavior(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path).rename(columns=COLUMNS)
    missing = set(REQUIRED) - set(data.columns)
    if missing:
        raise ValueError(f'Missing columns: {sorted(missing)}')
    data['mouse_id'] = data.mouse_id.astype('string').str.strip().str.lower()
    data['condition'] = data.condition.astype('string').str.strip().str.lower().replace({
        'ambient warming': 'ambtemp', 'warming': 'ambtemp', 'ambt': 'ambtemp',
        'bas': 'baseline', 'bl': 'baseline'})
    data['genotype'] = data.genotype.astype('string').str.strip().str.upper().replace({
        'APP': 'APP/PS1', 'APPPS1': 'APP/PS1'})
    if data[['mouse_id', 'condition']].isna().any().any() or (data.mouse_id == '').any():
        raise ValueError('Every test trial needs mouse_id and condition')
    if data.duplicated(['mouse_id', 'condition']).any():
        raise ValueError('Duplicate mouse_id + condition; resolve trial identity first')
    if not data.genotype.isin(GENOTYPES).all():
        raise ValueError('Unexpected genotype label; expected WT or APP/PS1')
    for col in ('tracking_pct', 'novel_time_s', 'familiar_time_s',
                'novel_entries', 'familiar_entries', 'start_entries',
                'total_distance_m'):
        if col in data:
            data[col] = pd.to_numeric(data[col], errors='coerce')
    denom = data.novel_time_s + data.familiar_time_s
    data['novel_preference_pct'] = np.where(
        denom > 0, 100 * data.novel_time_s / denom, np.nan)
    data['discrimination_index'] = np.where(
        denom > 0, (data.novel_time_s - data.familiar_time_s) / denom, np.nan)
    if {'novel_entries', 'familiar_entries'} <= set(data.columns):
        entries = data.novel_entries + data.familiar_entries
        data['novel_entry_preference_pct'] = np.where(
            entries > 0, 100 * data.novel_entries / entries, np.nan)
    if {'novel_entries', 'familiar_entries', 'start_entries'} <= set(data.columns):
        data['total_arm_entries'] = (
            data.novel_entries + data.familiar_entries + data.start_entries)
    return data


def apply_qc(data: pd.DataFrame, exclusions: Path | None,
             min_tracking: float = 90.0, min_exploration_s: float = 30.,
             expected_duration_s: float = 300.,
             max_zone_gap_fraction: float = .05) -> pd.DataFrame:
    data = data.copy()
    data['qc_reason'] = ''
    # Zone-time completeness is a *review warning*, not an automatic exclusion.
    # The four zones should approximately partition the nominal 300-s trial.
    # This assumption needs checking against the AnyMaze zone layout and video.
    zone_cols = ['start_time_s', 'novel_time_s', 'familiar_time_s', 'center_time_s']
    if all(c in data.columns for c in zone_cols):
        zones = data[zone_cols].apply(pd.to_numeric, errors='coerce')
        data['zone_time_sum_s'] = zones.sum(axis=1, min_count=len(zone_cols))
        data['zone_gap_s'] = expected_duration_s - data.zone_time_sum_s
        data['zone_coverage_pct'] = 100 * data.zone_time_sum_s / expected_duration_s
        data['zone_time_review_needed'] = (data.zone_time_sum_s.isna() |
            (data.zone_gap_s.abs() > expected_duration_s * max_zone_gap_fraction))
    else:
        data['zone_time_sum_s'] = np.nan
        data['zone_gap_s'] = np.nan
        data['zone_coverage_pct'] = np.nan
        data['zone_time_review_needed'] = False
    def flag(mask, reason):
        data.loc[mask, 'qc_reason'] += reason + ';'
    flag(data.tracking_pct.isna() | (data.tracking_pct < min_tracking), 'low_tracking')
    flag(data.novel_time_s.isna() | data.familiar_time_s.isna() |
         (data.novel_time_s < 0) | (data.familiar_time_s < 0), 'invalid_arm_time')
    flag((data.novel_time_s + data.familiar_time_s) < min_exploration_s,
         'insufficient_arm_exploration')
    # Explicit exclusions are recorded separately from automatic, objective QC.
    if exclusions and exclusions.exists():
        manual = pd.read_csv(exclusions).fillna('')
        if not {'mouse_id', 'condition', 'reason'} <= set(manual.columns):
            raise ValueError('Exclusions CSV needs mouse_id,condition,reason')
        for r in manual.itertuples(index=False):
            mid = str(r.mouse_id).strip().lower()
            cond = str(r.condition).strip().lower()
            mask = (data.mouse_id == mid) & ((data.condition == cond) | (cond == '*'))
            if not mask.any():
                raise ValueError(f'Exclusion does not match input: {mid}/{cond}')
            flag(mask, f'manual:{r.reason}')
    elif exclusions:
        raise FileNotFoundError(exclusions)
    data['included'] = data.qc_reason.eq('')
    return data


def exact_permutation_p(wt: np.ndarray, app: np.ndarray) -> float:
    """Two-sided exact label permutation test for difference in group means."""
    values = np.r_[wt, app]
    k = len(app)
    observed = abs(app.mean() - wt.mean())
    if not 0 < k < len(values):
        return float('nan')
    num, tail = 0, 0
    for ids in combinations(range(len(values)), k):
        other = np.ones(len(values), dtype=bool)
        other[list(ids)] = False
        diff = abs(values[list(ids)].mean() - values[other].mean())
        tail += diff >= observed - 1e-10
        num += 1
    return tail / num


def group_tables(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, tests = [], []
    for (condition, genotype), g in data.groupby(['condition', 'genotype']):
        for metric in BEHAVIOR:
            if metric not in g:
                continue
            v = g[metric].dropna()
            rows.append(dict(condition=condition, genotype=genotype, metric=metric,
                             n=len(v), mean=v.mean(), sd=v.std(ddof=1),
                             median=v.median()))
    for condition, g in data.groupby('condition'):
        for metric in BEHAVIOR:
            if metric not in g:
                continue
            wt = g.loc[g.genotype == 'WT', metric].dropna().to_numpy(float)
            app = g.loc[g.genotype == 'APP/PS1', metric].dropna().to_numpy(float)
            if min(len(wt), len(app)) < 2:
                continue
            tests.append(dict(condition=condition, metric=metric, n_wt=len(wt),
                              n_app=len(app), mean_diff_app_minus_wt=app.mean()-wt.mean(),
                              exact_permutation_p=exact_permutation_p(wt, app)))
    return pd.DataFrame(rows), pd.DataFrame(tests)


def plot_preferences(data: pd.DataFrame, path: Path) -> None:
    conditions = [x for x in ('baseline', 'ambtemp') if x in set(data.condition)]
    conditions += sorted(set(data.condition) - set(conditions))
    fig, ax = plt.subplots(figsize=(max(5, 1.8 * len(conditions) + 2), 4.6))
    rng = np.random.default_rng(42)
    for idx, condition in enumerate(conditions):
        for offset, geno, color in [(-.16, 'WT', 'gray'), (.16, 'APP/PS1', 'steelblue')]:
            vals = data[(data.condition == condition) & (data.genotype == geno)]
            x = np.full(len(vals), idx + offset) + rng.uniform(-.065, .065, len(vals))
            ax.scatter(x, vals.novel_preference_pct, s=42, alpha=.9, color=color,
                       label=geno if idx == 0 else None)
            for xx, y, mid in zip(x, vals.novel_preference_pct, vals.mouse_id):
                ax.annotate(mid, (xx, y), xytext=(3, 3),
                            textcoords='offset points', fontsize=7, alpha=.7)
    ax.axhline(50, linestyle='--', color='black', linewidth=1, alpha=.55,
               label='No preference (50%)')
    ax.set_xticks(range(len(conditions)), conditions)
    ax.set_ylabel('Novel-arm preference (% of novel + familiar time)')
    ax.set_ylim(0, 100)
    ax.legend(frameon=False, fontsize=9)
    ax.set_title('Y-maze recognition test — included mice')
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_behavior_diagnostics(data: pd.DataFrame, out: Path) -> None:
    """Show time- versus entry-based preference and locomotor activity.

    Diagnostic scatter plots only; no statistical outlier removal.
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for genotype, color in [('WT', 'gray'), ('APP/PS1', 'steelblue')]:
        g = data[data.genotype == genotype]
        axes[0].scatter(g.novel_preference_pct, g.novel_entry_preference_pct,
                        label=genotype, color=color, s=46)
        axes[1].scatter(g.total_distance_m, g.novel_preference_pct,
                        label=genotype, color=color, s=46)
        for _, row in g.iterrows():
            axes[0].annotate(f"{row.mouse_id} ({row.condition})",
                             (row.novel_preference_pct, row.novel_entry_preference_pct),
                             xytext=(3, 3), textcoords='offset points', fontsize=7)
            axes[1].annotate(f"{row.mouse_id} ({row.condition})",
                             (row.total_distance_m, row.novel_preference_pct),
                             xytext=(3, 3), textcoords='offset points', fontsize=7)
    axes[0].axhline(50, ls='--', color='black', alpha=.3)
    axes[0].axvline(50, ls='--', color='black', alpha=.3)
    axes[0].set(xlabel='Novel preference by time (%)',
                ylabel='Novel preference by arm entries (%)', xlim=(0, 100), ylim=(0, 100))
    axes[1].axhline(50, ls='--', color='black', alpha=.3)
    axes[1].set(xlabel='Total distance travelled (m)',
                ylabel='Novel preference by time (%)', ylim=(0, 100))
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle('Y-maze exploratory activity diagnostics (included records)')
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def link_spindles(data: pd.DataFrame, spindle_file: Path,
                  out: Path, min_n: int = 5) -> None:
    sp = pd.read_csv(spindle_file)
    keys = ['mouse_id', 'condition']
    if not set(keys) <= set(sp.columns):
        raise ValueError('Spindle table needs mouse_id and condition columns')
    sp['mouse_id'] = sp.mouse_id.astype('string').str.strip().str.lower()
    sp['condition'] = sp.condition.astype('string').str.strip().str.lower()
    if sp.duplicated(keys).any():
        raise ValueError('Spindle table must contain ONE row per mouse + condition')
    features = [c for c in sp.columns if c.startswith(('spindle_', 'pac_'))]
    if not features:
        raise ValueError('Expected spindle_* or pac_* numeric features')
    for col in features:
        sp[col] = pd.to_numeric(sp[col], errors='coerce')
    joined = data.merge(sp, on=keys, how='left', validate='one_to_one', indicator=True)
    joined.to_csv(out / 'spindle_join_audit.csv', index=False)
    # Correlation is calculated separately for each condition/genotype to avoid
    # repeated-animal pseudoreplication and genotype/condition confounding.
    rows = []
    for (condition, genotype), g in joined.groupby(['condition', 'genotype']):
        for feat in features:
            pair = g[['novel_preference_pct', feat]].dropna()
            n = len(pair)
            rho, p = (float('nan'), float('nan'))
            if n >= min_n and pair[feat].nunique() > 1 and pair.novel_preference_pct.nunique() > 1:
                rho, p = spearmanr(pair.novel_preference_pct, pair[feat])
            rows.append(dict(condition=condition, genotype=genotype, feature=feat,
                             matched_n=n, rho=rho, p_exploratory_uncorrected=p,
                             status='exploratory' if np.isfinite(rho) else 'insufficient_n_or_variation'))
    pd.DataFrame(rows).to_csv(out / 'spindle_correlations.csv', index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True, help='Mouse-level AnyMaze test-trial CSV')
    parser.add_argument('--output', type=Path, default=Path('data/processed/ymaze'))
    parser.add_argument('--exclusions', type=Path,
                        default=Path('config/ymaze_exclusions.csv'))
    parser.add_argument('--spindles', type=Path, help='Optional matched mouse/condition spindle features CSV')
    parser.add_argument('--trial-duration', type=float, default=300.,
                        help='Nominal test duration in seconds (default: 300)')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    raw = load_behavior(args.input)
    qc = apply_qc(raw, args.exclusions, expected_duration_s=args.trial_duration)
    qc.to_csv(args.output / 'qc_all_trials.csv', index=False)
    qc.loc[qc.zone_time_review_needed].to_csv(args.output / 'qc_zone_review.csv', index=False)
    clean = qc[qc.included].copy()
    clean.to_csv(args.output / 'ymaze_included.csv', index=False)
    qc.loc[~qc.included].to_csv(args.output / 'ymaze_excluded.csv', index=False)
    group, comparisons = group_tables(clean)
    group.to_csv(args.output / 'group_descriptives.csv', index=False)
    comparisons.to_csv(args.output / 'genotype_exploratory_permutation.csv', index=False)
    plot_preferences(clean, args.output / 'novel_preference.png')
    if 'novel_entry_preference_pct' in clean and 'total_distance_m' in clean:
        plot_behavior_diagnostics(clean, args.output / 'behavior_diagnostics.png')
    if args.spindles:
        link_spindles(clean, args.spindles, args.output)
    print(f'Test records {len(qc)}; included {len(clean)}; excluded {len(qc) - len(clean)}')
    print(qc.loc[~qc.included, ['mouse_id', 'condition', 'qc_reason']].to_string(index=False))
    n_review = int(qc.zone_time_review_needed.sum())
    print(f'Zone-time consistency: {n_review} record(s) flagged for manual review (not automatically excluded)')
    if n_review:
        print(qc.loc[qc.zone_time_review_needed,
                     ['mouse_id', 'condition', 'zone_time_sum_s', 'zone_gap_s', 'included']]
              .to_string(index=False))
    print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
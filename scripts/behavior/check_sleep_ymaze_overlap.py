#!/usr/bin/env python3
"""Audit mouse-level Y-maze vs spindle-probability cohort overlap.

This does NOT estimate correlations by itself. In particular, a shared mouse ID
is necessary but not sufficient: recording and behavioral timing must also be
verified before interpreting a memory/sleep association.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path
import pandas as pd

SPINDLE_COLUMNS = (
    'spindle_pre_rem_occupancy_25_15',
    'spindle_pre_wake_occupancy_25_15',
    'spindle_prerem_minus_prewake_25_15',
)
PHASES = {('baseline', 'BASELINE_ZT3_9'): 'baseline',
          ('warming', 'WARM_ACTIVE'): 'ambtemp'}
TRANSITIONS = {'NREM_to_REM': 'spindle_pre_rem_occupancy_25_15',
               'NREM_to_WAKE': 'spindle_pre_wake_occupancy_25_15'}


def normalize_mouse_id(value: object) -> str:
    """Represent `m69`, `mouse69`, `M-069` as `m69`. Do not guess other IDs."""
    s = str(value).strip().lower()
    match = re.fullmatch(r'(?:mouse|m)[\s_-]*0*(\d+)', s)
    return f'm{int(match.group(1))}' if match else s


def normalize_genotype(value: object) -> str:
    s = str(value).strip().upper().replace('/', '').replace('-', '')
    aliases = {'APP': 'APP/PS1', 'APPPS1': 'APP/PS1', 'WT': 'WT', 'WILDTYPE': 'WT'}
    if s not in aliases:
        raise ValueError(f'Unknown genotype: {value!r}')
    return aliases[s]


def load_behavior(path: Path) -> pd.DataFrame:
    b = pd.read_csv(path)
    required = {'mouse_id', 'condition', 'genotype', 'novel_preference_pct'}
    if not required <= set(b.columns):
        raise ValueError(f'Behavior missing: {sorted(required - set(b.columns))}')
    if 'included' in b.columns:
        b = b[b['included'].astype(str).str.lower().eq('true')].copy()
    b['mouse_id'] = b.mouse_id.map(normalize_mouse_id)
    b['condition'] = b.condition.astype(str).str.strip().str.lower()
    b['genotype'] = b.genotype.map(normalize_genotype)
    if b.duplicated(['mouse_id', 'condition']).any():
        raise ValueError('Behavior has repeated mouse+condition; fix trial identity')
    return b[['mouse_id', 'condition', 'genotype', 'novel_preference_pct']].copy()


def load_spindle_probability(path: Path) -> pd.DataFrame:
    s = pd.read_csv(path)
    required = {'mouse_id', 'genotype', 'condition', 'phase', 'transition_type',
                'occupancy_probability_25_15'}
    if not required <= set(s.columns):
        raise ValueError(f'Spindle export missing: {sorted(required - set(s.columns))}')
    s['mouse_id'] = s.mouse_id.map(normalize_mouse_id)
    s['genotype'] = s.genotype.map(normalize_genotype)
    s['condition'] = s.condition.astype(str).str.strip().str.lower()
    s['phase'] = s.phase.astype(str).str.strip()
    s['transition_type'] = s.transition_type.astype(str).str.strip()
    s['behavior_condition'] = [PHASES.get((c, p)) for c, p in zip(s.condition, s.phase)]
    s = s[s.behavior_condition.notna() & s.transition_type.isin(TRANSITIONS)].copy()
    keys = ['mouse_id', 'behavior_condition', 'transition_type']
    if s.duplicated(keys).any():
        raise ValueError('Spindle summary has multiple rows per mouse/phase/transition')
    s['feature'] = s.transition_type.map(TRANSITIONS)
    s['occupancy_probability_25_15'] = pd.to_numeric(
        s.occupancy_probability_25_15, errors='coerce')
    # The source values are fractions, not percentages.
    if s.occupancy_probability_25_15.dropna().between(0, 1).eq(False).any():
        raise ValueError('Unexpected spindle occupancy outside [0, 1]')
    grouped = (s.pivot(index=['mouse_id', 'behavior_condition'],
                       columns='feature', values='occupancy_probability_25_15')
                 .reset_index().rename(columns={'behavior_condition': 'condition'}))
    for c in SPINDLE_COLUMNS[:2]:
        if c not in grouped:
            grouped[c] = float('nan')
    grouped[SPINDLE_COLUMNS[2]] = (
        grouped[SPINDLE_COLUMNS[0]] - grouped[SPINDLE_COLUMNS[1]])
    geno = s.groupby(['mouse_id', 'behavior_condition']).genotype.nunique()
    if (geno > 1).any():
        raise ValueError('Conflicting spindle genotype for a mouse+condition')
    one_geno = (s.drop_duplicates(['mouse_id', 'behavior_condition'])
                [['mouse_id', 'behavior_condition', 'genotype']]
                .rename(columns={'behavior_condition': 'condition', 'genotype': 'genotype_sleep'}))
    grouped = grouped.merge(one_geno, on=['mouse_id', 'condition'], validate='one_to_one')
    return grouped


def build_audit(behavior: pd.DataFrame, spindle: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    audit = behavior.merge(spindle, on=['mouse_id', 'condition'], how='left',
                           indicator=True, validate='one_to_one')
    mismatch = audit.genotype_sleep.notna() & (audit.genotype != audit.genotype_sleep)
    if mismatch.any():
        raise ValueError('Genotype mismatch: ' + ', '.join(audit.loc[mismatch, 'mouse_id']))
    audit['spindle_match_status'] = audit['_merge'].map(
        {'both': 'matched_id_and_condition', 'left_only': 'no_matching_sleep_measurement',
         'right_only': 'unexpected'}).astype(str)
    audit = audit.drop(columns=['_merge'])
    matched = audit[audit.spindle_match_status == 'matched_id_and_condition'].copy()
    return audit, matched


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--behavior', required=True, type=Path,
                    help='data/processed/ymaze/ymaze_included.csv')
    ap.add_argument('--spindles', required=True, type=Path,
                    help='mouse_probability_summary_25s.csv from sleep_app_ach_warm')
    ap.add_argument('--output', type=Path, default=Path('data/processed/ymaze/sleep_matching'))
    args = ap.parse_args()
    b = load_behavior(args.behavior)
    s = load_spindle_probability(args.spindles)
    audit, matched = build_audit(b, s)
    args.output.mkdir(parents=True, exist_ok=True)
    audit.to_csv(args.output / 'ymaze_sleep_matching_audit.csv', index=False)
    feature_cols = ['mouse_id', 'condition', *SPINDLE_COLUMNS]
    matched[feature_cols].to_csv(args.output / 'spindle_features_matched.csv', index=False)
    print(f'Y-maze records: {len(b)}, spindle mouse/condition records: {len(s)}')
    print(f'Y-maze records with matching spindle mouse/condition: {len(matched)}')
    print(f'Unique behavior mice: {b.mouse_id.nunique()}; unique sleep mice: {s.mouse_id.nunique()}')
    if matched.empty:
        print('NO OVERLAP: no spindle–memory correlations can be performed yet.')
    else:
        print('Matched IDs must still be verified for recording/test dates and timing.')
    print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
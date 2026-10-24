#!/usr/bin/env python3
"""Analyze 60-second ANY-maze Data-page segment exports, one row per test segment.

NOT a parser for proprietary .szd/.szk/.szv files. Uses exported CSV only.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def canonical(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(s).lower()).strip()


def choose_col(columns: list[str], phrases: tuple[str, ...], required: bool = True) -> str | None:
    names = {canonical(x): x for x in columns}
    for phrase in phrases:
        key = canonical(phrase)
        if key in names:
            return names[key]
    # Allow ANY-maze versions to add bracketed units or punctuation, but not unrelated measures
    for phrase in phrases:
        key = canonical(phrase)
        hits = [orig for norm,orig in names.items() if norm.startswith(key + ' ') or norm == key]
        if len(hits) == 1:
            return hits[0]
    if required:
        raise ValueError(f"Missing column for {phrases}; available columns: {columns}")
    return None


def detect_csv(path: Path) -> pd.DataFrame:
    text = None
    for encoding in ('utf-8-sig', 'utf-16', 'cp1252'):
        try:
            candidate = path.read_text(encoding=encoding)
            if '\x00' in candidate[:200]:
                continue
            text = candidate
            break
        except (UnicodeError, UnicodeDecodeError):
            continue
    if text is None:
        raise ValueError(f"Cannot decode {path}. Export a UTF-8 CSV from ANY-maze.")
    try:
        dialect = csv.Sniffer().sniff(text[:8192], delimiters=',;\t')
        sep = dialect.delimiter
    except csv.Error:
        sep = ','
    df = pd.read_csv(path, sep=sep, encoding=encoding)
    if df.shape[1] <= 1:
        raise ValueError(f"Only {df.shape[1]} column found. Export the Data page as CSV or tab-separated text.")
    return df


def mouse_id(value: object) -> str:
    match = re.search(r"(?:mouse|m)\s*[-_ ]*0*(\d+)\b", str(value), re.I)
    if match:
        return f"m{int(match.group(1))}"
    # Numeric animal IDs are ambiguous with ANY-maze's automatic animal number.
    raise ValueError(f"Unrecognized animal ID: {value!r}; use explicit IDs such as m78")


def condition(value: object) -> str:
    s = canonical(str(value))
    if s in ('baseline', 'base', 'control') or 'baseline' in s:
        return 'baseline'
    if s in ('ambtemp', 'amb temp', 'ambient temperature', 'warming', 'warm') or any(q in s for q in ['ambtemp','ambient temp','warming']):
        return 'ambtemp'
    raise ValueError(f"Unrecognized condition: {value!r}. Use baseline or ambtemp.")


def parse_segment(value: object) -> tuple[float,float]:
    s = str(value).lower().replace(',', '.')
    numbers = re.findall(r"\d+(?:\.\d+)?", s)
    if len(numbers) < 2:
        raise ValueError(f"Segment {value!r} does not state start and end times; expected '0-60 seconds'.")
    a,b = float(numbers[0]), float(numbers[1])
    if b <= a:
        raise ValueError(f"Invalid segment {value!r}")
    return a,b


def float_series(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values.astype(str).str.replace(',', '.', regex=False), errors='coerce')


def parse_export(df: pd.DataFrame) -> pd.DataFrame:
    cols = list(df.columns)
    animal = choose_col(cols, ('Animal ID','Animal','Mouse ID'))
    treatment = choose_col(cols, ('Treatment','Condition','Group'))
    segment = choose_col(cols, ('Segment of test','Test segment','Segment'))
    novel_t = choose_col(cols, ('Time in the novel arm zone (s)','Time in the novel arm zone','Time in novel arm zone','Novel time s'))
    familiar_t = choose_col(cols, ('Time in the familiar arm zone (s)','Time in the familiar arm zone','Time in familiar arm zone','Familiar time s'))
    entry_n = choose_col(cols, ('Number of entries to the novel arm zone','Number of entries to novel arm zone','Novel entries'), required=False)
    entry_f = choose_col(cols, ('Number of entries to the familiar arm zone','Number of entries to familiar arm zone','Familiar entries'), required=False)
    dist = choose_col(cols, ('Total distance travelled (m)','Total distance travelled','Distance travelled m'), required=False)
    mobile = choose_col(cols, ('Total time mobile (s)','Total time mobile','Time mobile s'), required=False)
    track = choose_col(cols, ('Percentage of frames tracked (%)','Percentage of frames tracked'), required=False)
    r = pd.DataFrame()
    r['mouse_id'] = df[animal].map(mouse_id)
    r['condition'] = df[treatment].map(condition)
    ranges = df[segment].map(parse_segment)
    r['start_s'] = ranges.map(lambda x: x[0])
    r['end_s'] = ranges.map(lambda x: x[1])
    for name, c in [('novel_s',novel_t),('familiar_s',familiar_t),('novel_entries',entry_n),('familiar_entries',entry_f),('distance_m',dist),('mobile_s',mobile),('tracked_pct',track)]:
        r[name] = float_series(df[c]) if c else np.nan
    required = ['novel_s','familiar_s']
    if r[required].isna().any().any():
        raise ValueError('Novel/familiar zone time has missing/non-numeric values. Check CSV format and decimal separators.')
    return r


def exclude(df: pd.DataFrame, exclusions_path: Path) -> pd.DataFrame:
    if not exclusions_path.is_file():
        raise FileNotFoundError(exclusions_path)
    x = pd.read_csv(exclusions_path, dtype=str).fillna('')
    keep = pd.Series(True,index=df.index)
    for _,row in x.iterrows():
        try:
            mid=mouse_id(row['mouse_id'])
        except ValueError:
            continue
        cond=canonical(row.get('condition','*'))
        mask=(df['mouse_id']==mid)
        if cond not in ('*','all',''):
            mask &= (df['condition']==condition(cond))
        keep &= ~mask
    return df.loc[keep].copy()


def with_metadata(df: pd.DataFrame, metadata_path: Path | None) -> pd.DataFrame:
    if metadata_path is None:
        df = df.copy()
        df['genotype'] = 'unknown'
        return df
    meta = pd.read_csv(metadata_path)
    for key in ('mouse_id', 'condition', 'genotype'):
        if key not in meta.columns:
            raise ValueError(f"Metadata {metadata_path} lacks {key}")
    meta = meta[['mouse_id','condition','genotype']].drop_duplicates()
    meta['mouse_id'] = meta['mouse_id'].map(mouse_id)
    meta['condition'] = meta['condition'].map(condition)
    meta['genotype'] = meta['genotype'].astype(str).str.strip().str.upper().replace({'APP':'APP/PS1', 'APPPS1':'APP/PS1', 'APP/PS1':'APP/PS1'})
    if meta.duplicated(['mouse_id','condition']).any():
        raise ValueError('Conflicting genotype assignments in metadata for the same mouse-condition')
    result = df.merge(meta,on=['mouse_id','condition'],how='left',validate='many_to_one')
    if result.genotype.isna().any():
        missing=result[result.genotype.isna()][['mouse_id','condition']].drop_duplicates().to_dict('records')
        raise ValueError(f"Unmatched animal-condition IDs in metadata: {missing}")
    return result


def build_outputs(df: pd.DataFrame) -> tuple[pd.DataFrame,pd.DataFrame,pd.DataFrame]:
    if df.duplicated(['mouse_id','condition','start_s','end_s']).any():
        raise ValueError('Duplicate mouse-condition segment. Separate repeat trials before analysis.')
    out = df.sort_values(['mouse_id','condition','start_s']).copy()
    out['duration_s']=out['end_s']-out['start_s']
    out['preference_pct']=100*out['novel_s']/(out['novel_s']+out['familiar_s']).replace(0,np.nan)
    out['entry_preference_pct']=100*out['novel_entries']/(out['novel_entries']+out['familiar_entries']).replace(0,np.nan)
    out['zone_time_overrun']=out['novel_s']+out['familiar_s'] > out['duration_s']+0.5
    out['any_time_outside_preference_arms']=(out['novel_s']+out['familiar_s'] < out['duration_s']-0.5)
    out['valid_segment']=~out['zone_time_overrun']
    out['expected_60s_bin']=np.isclose(out['duration_s'],60,atol=.1) & np.isclose(out['start_s']/60,np.round(out['start_s']/60),atol=.002)
    # Guard against partial/missing bins. A mouse with incomplete bins is not included in cumulative 300 s.
    if 'genotype' not in out:
        out['genotype']='unknown'
    cumulative=[]
    for (mid,cond),g in out.groupby(['mouse_id','condition']):
        for end in (60,120,300):
            q=g[(g.start_s >= 0) & (g.end_s <= end+0.01)].sort_values('start_s')
            contiguous=len(q)>0 and abs(q.iloc[0].start_s)<.01 and abs(q.iloc[-1].end_s-end)<.01 and np.allclose(q.start_s.iloc[1:],q.end_s.iloc[:-1].to_numpy(),atol=.01)
            enough= q['valid_segment'].all() and q['expected_60s_bin'].all()
            if not contiguous or not enough:
                cumulative.append({'mouse_id':mid,'condition':cond,'window_end_s':end,'genotype':g.genotype.iloc[0],'novel_preference_pct':np.nan,'novel_entry_preference_pct':np.nan,'novel_time_s':np.nan,'familiar_time_s':np.nan,'complete':False})
                continue
            tnov=q.novel_s.sum(); tfam=q.familiar_s.sum()
            en=q.novel_entries.sum(min_count=1); ef=q.familiar_entries.sum(min_count=1)
            cumulative.append({'mouse_id':mid,'condition':cond,'window_end_s':end,'genotype':g.genotype.iloc[0],
                              'novel_preference_pct':100*tnov/(tnov+tfam) if tnov+tfam else np.nan,
                              'novel_entry_preference_pct':100*en/(en+ef) if np.isfinite(en+ef) and en+ef else np.nan,
                              'novel_time_s':tnov,'familiar_time_s':tfam,'complete':True})
    cum=pd.DataFrame(cumulative)
    summary=(cum.dropna(subset=['novel_preference_pct']).groupby(['condition','genotype','window_end_s']).agg(
        n=('mouse_id','size'),mean_preference_pct=('novel_preference_pct','mean'),
        sd_preference_pct=('novel_preference_pct','std')).reset_index())
    return out,cum,summary


def plot(out: pd.DataFrame, cumulative: pd.DataFrame, dst: Path):
    fig,axs=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for c,ax in zip(('baseline','ambtemp'),axs):
        dd=out[out.condition==c]
        for mid,g in dd.groupby('mouse_id'):
            xx=(g.start_s+g.end_s)/2/60
            ax.plot(xx,g.preference_pct,marker='o',alpha=.55,linewidth=1,linestyle='--' if g.genotype.iloc[0]=='APP/PS1' else '-',label=f'{mid} ({g.genotype.iloc[0]})')
        ax.axhline(50,color='0.35',ls='--',lw=1)
        ax.set(title=f'{c}: per-minute novel preference',xlabel='Minute of test',ylabel='Novel preference (%)',ylim=(0,100),xticks=(.5,1.5,2.5,3.5,4.5))
        ax.grid(alpha=.15)
    fig.savefig(dst,dpi=170);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for c,ax in zip(('baseline','ambtemp'),axs):
        dd=cumulative[(cumulative.condition==c)&(cumulative.complete)]
        for mid,g in dd.groupby('mouse_id'):
            ax.plot(g.window_end_s,g.novel_preference_pct,marker='o',alpha=.7,linewidth=1,linestyle='--' if g.genotype.iloc[0]=='APP/PS1' else '-',label=f'{mid} ({g.genotype.iloc[0]})')
        ax.axhline(50,color='0.35',ls='--',lw=1)
        ax.set(title=f'{c}: cumulative novel preference',xlabel='Elapsed time (s)',ylabel='Novel preference (%)',ylim=(0,100),xticks=(60,120,300))
        ax.grid(alpha=.15)
    fig.savefig(dst.with_name('cumulative_preference.png'),dpi=170);plt.close(fig)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True,type=Path,help='ANY-maze Data page export with one row per 60-second segment')
    p.add_argument('--exclusions',required=True,type=Path)
    p.add_argument('--metadata',type=Path,help='Validated included Y-maze mouse-level CSV with mouse_id, condition, genotype')
    p.add_argument('--output',required=True,type=Path)
    a=p.parse_args(argv)
    df=parse_export(detect_csv(a.input))
    cleaned=exclude(df,a.exclusions)
    cleaned=with_metadata(cleaned,a.metadata)
    out,cum,summary=build_outputs(cleaned)
    a.output.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.output/'time_bins_mouse_level.csv',index=False)
    cum.to_csv(a.output/'early_60_120_300s_mouse_level.csv',index=False)
    summary.to_csv(a.output/'timecourse_descriptives.csv',index=False)
    plot(out,cum,a.output/'minute_by_minute_preference.png')
    print(f'Imported {len(df)} rows, retained {len(out)} after exclusions, {len(cum)//3} mouse-condition records; output {a.output}')
    missing=cum[~cum.complete]
    if len(missing):
        print(f'WARNING: {len(missing)} cumulative windows incomplete; see early_60_120_300s_mouse_level.csv')
    if out.zone_time_overrun.any():
        print(f'WARNING: {(out.zone_time_overrun).sum()} segments have impossible overlapping zone occupancy. Validate AnyMaze zones.')

if __name__ == '__main__':
    main()
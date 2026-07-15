# Example input audit

This audit documents the example files used to bootstrap the loader. It is not a biological inclusion decision.

## Baseline recording

`300526-m63-bas-1.mat`

- EEG: 22,284,863 samples at 610.3515625 Hz
- EMG: 22,284,863 samples at 610.3515625 Hz
- Processed photometry variable: `ne`
- Photometry: 37,141,375 samples at 1017.2526245 Hz
- EEG/EMG duration: approximately 36,511.520 s
- Photometry duration: approximately 36,511.457 s
- Maximum signal-duration mismatch: approximately 0.063 s
- Complete common one-second epochs: 36,511

## Scoring export

`final_scoring.csv`

- Rows: 36,532
- Epoch resolution: 1 s
- Codes: Undefined `-1`, Wake `0`, NREM `1`, REM `2`
- The scoring extends approximately 20.5 s beyond the continuous signals.
- The pipeline therefore retains 36,511 complete epochs and records 21 trimmed scoring rows in QC.
- Microarousals are derived later; they are not present in the original scoring export.

Before group analysis, confirm that this scoring file belongs to this recording and complete genotype, sex, and age metadata in the manifest.

## Expected scoring-derived QC for this pair

After retaining 36,511 complete epochs and applying the strict `<15 s` microarousal rule:

- Original Wake seconds: 11,672
- NREM seconds: 20,766
- REM seconds: 4,072
- Undefined seconds: 1
- Derived microarousal bouts: 44
- Derived microarousal seconds: 325
- Exactly 15-second Wake bouts retained as Wake: 12

These values provide a useful first regression check after the initial MATLAB run.

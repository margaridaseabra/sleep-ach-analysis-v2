# Data dictionary

## Session manifest

| Field | Meaning |
|---|---|
| `session_id` | Stable standardized identifier |
| `date` | Recording date in ISO format |
| `mouse_id` | Animal identifier |
| `genotype` | `WT`, `APPPS1`, or temporarily `UNKNOWN` |
| `sex` | `F`, `M`, or temporarily `UNKNOWN` |
| `age_days` | Age on recording date |
| `condition` | `baseline`, `ambtemp`, or `donepezil` |
| `replicate` | Session number for the same mouse/condition |
| `recording_file` | Relative or absolute MAT-file path |
| `scoring_file` | Relative or absolute scoring CSV path |
| `photometry_variable` | MAT variable containing the processed photometry trace |
| `photometry_is_preprocessed` | Whether the trace already represents a corrected signal |
| `temperature_file` | Optional temperature CSV |
| `temperature_offset_s` | Value added to logger time to align with recording time |
| `zt_start` | Zeitgeber time at recording start |
| `include_session` | Whether the pipeline should process the row |
| `notes` | Free-text metadata or warnings |

## Standardized MAT output

```text
eeg_200                  single, 200 Hz
emg_200                  single, 200 Hz
ach_1hz                  double, 1 Hz
scores_original          int16, 1 Hz
scores_analysis          int16, 1 Hz
scores_nrem_like         int16, 1 Hz
metadata                 struct
qc                       one-row table
temperature_aligned      table, optional
```

## Scoring CSV output

The exported `scores_analysis.csv` contains all aligned original scoring columns plus:

| Field | Meaning |
|---|---|
| `original_code` | Duration-aligned original score |
| `analysis_code` | Score after microarousal derivation |
| `nrem_like_code` | Score with microarousals merged into NREM |

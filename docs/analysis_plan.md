# Analysis plan

## Phase 1 — preprocessing and QC

- Validate MAT variables and sampling rates.
- Align EEG, EMG, ACh, scoring, and temperature.
- Preserve review provenance and undefined epochs.
- Generate one QC report per session.
- Freeze one baseline and one warming session as reference outputs.

## Phase 2 — sleep architecture

- State percentage by 3-hour bin and by hour.
- Bout count and duration.
- Microarousal burden and NREM fragmentation index.
- Kaplan–Meier and mouse-clustered Cox models for NREM and REM bouts.

## Phase 3 — state-resolved EEG

- State PSD and band power.
- Spindle detection and features.
- Pre-transition spindle probability.

## Phase 4 — ACh

- State and transition QC.
- NREM low-frequency PSD split by subsequent outcome.
- Transition-aligned peak latency, AUC, and rate of rise.
- NREM variability, bout range, and stability features.

## Phase 5 — exploratory models

- ACh feature survival models.
- Phase–amplitude coupling.
- Cohort/batch sensitivity analyses.

The new recordings should first be treated as a replication cohort. Combine them with the thesis cohort only after preprocessing comparability and batch effects have been assessed.

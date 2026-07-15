# Sleep–ACh analysis pipeline v2

Clean, reproducible preprocessing and analysis scaffold for EEG, EMG, sleep scoring, cortical ACh photometry, and ambient-temperature recordings from the APP/PS1 sleep project.

The repository is designed to reproduce the core processing decisions used in the MSc thesis while keeping raw data, analysis-ready data, and source code separate.

## What the first preprocessing stage does

For every session listed in `config/sessions.csv`, the pipeline:

1. inspects the raw MAT-file without modifying it;
2. checks EEG, EMG, and photometry duration and sampling rates;
3. imports the full scoring export and preserves review provenance;
4. trims scoring and signals to the common complete-second duration;
5. derives microarousals from Wake runs shorter than 15 seconds;
6. creates a separate NREM-like score vector in which microarousals are merged into NREM;
7. applies a 50 Hz notch to EEG/EMG and resamples them to 200 Hz;
8. converts the already-processed photometry trace to a robust 1 Hz median trace;
9. aligns and summarizes a temperature log when present;
10. saves standardized outputs and session-level QC reports.

Raw data are never edited.

## Repository layout

```text
sleep-ach-analysis-v2/
├── config/                 session manifest and analysis settings
├── data/
│   ├── raw/                local raw files; ignored by Git
│   ├── interim/            temporary outputs; ignored by Git
│   └── processed/          standardized outputs; ignored by Git
├── docs/                   decisions, data dictionary, and analysis plan
├── examples/               synthetic demo-data generator
├── scripts/                Git and data-staging helpers
├── src/                    reusable MATLAB functions
├── tests/                  lightweight MATLAB tests
├── project_setup.m         adds project folders to the MATLAB path
└── run_preprocessing.m     main entry point
```

## Requirements

- MATLAB R2025b recommended
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox is recommended for later analyses

The initial preprocessing code uses `iirnotch`, `filtfilt`, and `resample` from Signal Processing Toolbox.

## 1. Put the repository on your computer

Unzip the starter repository, then open a terminal inside it.

```bash
bash scripts/initialize_git.sh
```

This initializes a local Git repository and creates the first commit using your own configured Git identity.

## 2. Create the GitHub repository

With GitHub CLI installed and authenticated:

```bash
gh repo create sleep-ach-analysis-v2 --private --source=. --remote=origin --push
```

Alternatively, create an empty private repository on GitHub and run:

```bash
git remote add origin git@github.com:YOUR_USERNAME/sleep-ach-analysis-v2.git
git push -u origin main
```

Do not upload raw animal recordings to GitHub.

## 3. Stage the example recording locally

For the uploaded baseline example:

```bash
bash scripts/stage_baseline_example.sh \
  /path/to/300526-m63-bas-1.mat \
  /path/to/final_scoring.csv
```

For the uploaded temperature example:

```bash
bash scripts/stage_temperature_example.sh \
  /path/to/temperature_log_20260530_185347_m61-ambt-2.csv
```

The scripts copy data under `data/raw/`, which is excluded by `.gitignore`.

## 4. Complete the session manifest

Edit `config/sessions.csv` before cohort analysis. In particular, fill in genotype, sex, age, and confirm that each scoring file belongs to the listed recording.

The first example row is configured for:

```text
300526-m63-bas-1.mat
final_scoring.csv
```

## 5. Run preprocessing

Open MATLAB in the repository root and run:

```matlab
project_setup
run_preprocessing
```

Outputs appear in:

```text
data/processed/<session_id>/
```

Expected outputs:

```text
preprocessed_session.mat
scores_analysis.csv
qc_summary.csv
qc_report.md
temperature_aligned.csv   % warming sessions only
```

## Important scoring convention

The authoritative scoring remains the imported `final_code` vector:

```text
-1 Undefined
 0 Wake
 1 NREM
 2 REM
```

The pipeline produces separate derived vectors:

```text
analysis_code   Wake runs <15 s become microarousal code 15
nrem_like_code  microarousal code 15 is temporarily merged into NREM code 1
```

Exactly 15-second Wake bouts remain Wake.

## Important photometry convention

The supplied MAT-files contain a processed photometry variable named `ne`, rather than separate raw 465 nm and 405 nm channels. The pipeline therefore preserves it as a processed signal and creates a 1 Hz median trace. It does not claim to repeat 405-to-465 fitting or motion/hemodynamic correction.

## First development milestones

1. Run the preprocessing pipeline on one baseline session.
2. Review the generated QC report and scoring warnings.
3. Add one warming session and verify temperature alignment.
4. Freeze reference outputs from one thesis-era session.
5. Rebuild sleep architecture analyses against standardized inputs.
6. Add ACh, spindle, survival, and PAC modules only after the preprocessing outputs are stable.

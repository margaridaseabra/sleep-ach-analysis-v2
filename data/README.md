# Local data layout

Raw and derived data are intentionally excluded from Git.

Use one folder per session:

```text
data/raw/20260530_m063_baseline_r01/
├── 300526-m63-bas-1.mat
└── final_scoring.csv
```

Warming sessions can additionally contain:

```text
temperature_log_YYYYMMDD_HHMMSS_mXX-ambt-X.csv
```

The manifest in `config/sessions.csv` is the authoritative map between files and biological metadata.

# Cross-cohort matching audit: Y-maze recognition and spindles

The AnyMaze behavior data concern a novel-arm recognition protocol (training with the novel arm blocked, then test with it opened). Behavioral condition `ambtemp` is not proof the test followed learning and sleep manipulation in an experimental sequence appropriate for a memory-consolidation causal claim.

The mouse-level spindle probability export lives in the older `sleep_app_ach_warm` repository, on `paper/pac-comodulograms-2026-10-23`, at `results/prism/spindles/probability/mouse_probability_summary_25s.csv`.

## Usage

From `~/Downloads/sleep-ach-analysis-v2`:

```bash
conda activate ymaze-analysis
python scripts/behavior/check_sleep_ymaze_overlap.py \
  --behavior data/processed/ymaze/ymaze_included.csv \
  --spindles ~/sleep_app_ach_warm/results/prism/spindles/probability/mouse_probability_summary_25s.csv \
  --output data/processed/ymaze/sleep_matching
```

Outputs: `ymaze_sleep_matching_audit.csv` includes every included behavioral test and whether there is a same-mouse and same-condition spindle metric. `spindle_features_matched.csv` is an analysis-ready one-row-per-mouse-and-condition file **only for matching records**.

These are explicitly selected spindle windows:

- behavioral `baseline` ↔ sleep `condition=baseline`, `phase=BASELINE_ZT3_9`
- behavioral `ambtemp` ↔ sleep `condition=warming`, `phase=WARM_ACTIVE`
- events: `NREM_to_REM` and `NREM_to_WAKE`, using *fractional* `occupancy_probability_25_15` (not percent). The contrast is preREM minus preWake.

Do not interpret these as trial-timed sleep recordings without checking the actual test date, sleep date, and learning/retention intervals. **Even if mouse IDs overlap, same-animal sleep data alone do not prove consolidation.**

## Known current limitation, October 2026

The available spindle-probability and PAC exports each contain 16 unique mice (`mouse1,2,3,4,8,10,11,12,59,61,63,64,66,68,70,72`). The retained 17 behavioral test records involve 14 other unique mice (`m65,69,74,76,77,78,79,80,81,82,83,84,85,87`). Their overlap is zero. The audit must not quietly calculate a correlation or substitute a group mean.

Next: obtain / process sleep EEG for animals actually in the Y-maze cohort (if recorded), verify genotype and test timing, then rerun matching. If those recordings don't exist, treat Y-maze and spindle results as **independent cohorts** and report only parallel group-level phenotypes, not within-animal correlations.
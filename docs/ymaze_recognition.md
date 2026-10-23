# Y-maze novel-arm recognition: behavioral and spindle association pipeline

This pipeline analyzes the **test trial** of a two-trial Y-maze task. During familiarization the novel arm was blocked; it was then made accessible for the recognition test. It does **not** compute spontaneous alternation.

## Requirements

Python 3.10+; install `numpy`, `pandas`, `scipy`, `matplotlib` (e.g. `python -m pip install numpy pandas scipy matplotlib`). The rest of the repository's MATLAB processing is unaffected.

## First analysis (no spindle file required)

Put `Ymaze_AnyMaze_all_measures.csv` or `Ymaze_AnyMaze_preliminary_summary.csv` **locally** in `data/raw/ymaze/`. These are per-animal AnyMaze test-trial exports with `mouse_id`, `condition`, and `genotype` columns. The script accepts either the original exported column names or the normalized short names.

From the repository root:

```bash
python scripts/behavior/run_ymaze.py \
  --input data/raw/ymaze/Ymaze_AnyMaze_all_measures.csv \
  --output data/processed/ymaze \
  --exclusions config/ymaze_exclusions.csv
```

Automatically exclude tracking `<90%` and total novel+familiar exploration `<30 s`, or invalid arm times; record all reasons in `qc_all_trials.csv`. The explicit exclusions CSV removes **both** m73 records pending original-video verification because their measurements are virtually identical across conditions. The low-tracking m91 baseline record is removed automatically. Low preference is **never** a reason to exclude an animal. Manual exclusions must be documented, not chosen based on resulting genotype effects or correlations.

### Zone-time review (added after examining the 23 AnyMaze records)

With a nominal 300-second test, the sum of **novel + familiar + start + center** time should approximately account for the trial *if* the four AnyMaze zones form a nonoverlapping partition. The pipeline now computes `zone_time_sum_s`, `zone_gap_s`, `zone_coverage_pct`, and a `zone_time_review_needed` flag; the default trigger is more than 5% of expected duration (15 seconds). You can change the assumed trial duration with `--trial-duration 300`. **Zone-time flags are warnings only; they do not exclude mice.** Inspect the zone layout, the original video, actual trial length, and source reports before deciding whether a missing/excess zone sum is a real tracking artifact.

In the supplied test set, the warnings affect: m89 baseline APP/PS1 (104.1 s accounted for), m67 baseline APP/PS1 (186.9 s), m78 ambient APP/PS1 (263.6 s), and both already-excluded m73 records (~319 s; likely overlapping zones or duplicated processing). Their high frame-tracking percentages do not establish that arm-zone measures are valid. Removing m67 and m89 immediately would leave just one APP/PS1 mouse at baseline; do not drop them based only on the sums.

The additional `behavior_diagnostics.png` contrasts time-based versus entry-based novelty preference, and distance travelled versus time preference. Long visits or low movement can explain discrepancies between the two preference measures (e.g. m79 baseline WT: about 9% of arm time but 71% of novel+familiar entries), without constituting an automatic exclusion criterion.

Output files: `qc_all_trials.csv`, `qc_zone_review.csv`, `ymaze_included.csv`, `ymaze_excluded.csv`, `group_descriptives.csv`, `genotype_exploratory_permutation.csv`, `novel_preference.png`, and `behavior_diagnostics.png`. Re-running replaces derived files, not source data. The statistical permutation test compares genotypes **within** each condition only. Different conditions share some mice, so these cross-condition records cannot be treated as independent.

## Outcome definitions

Primary: `100 * novel_time_s / (novel_time_s + familiar_time_s)` (50% is no preference). Secondary: `(novel - familiar) / (novel + familiar)`; entry-based novelty preference; distance traveled, total entries, tracking completeness. Inspect video and zone definitions for anomalous movement before drawing conclusions. No post-hoc exclusion of poor learners.

## Connecting to spindle exports

Create a **mouse-level** spindle table `data/raw/ymaze/spindle_features.csv`, one row per `mouse_id,condition`, with numeric feature columns beginning `spindle_` or `pac_`, for example:

```csv
mouse_id,condition,spindle_density_nrem,spindle_pre_rem_probability,pac_nrem_prerem
m65,baseline,4.2,0.19,0.37
```

The numbers above are **examples**, not measured values. Export actual mouse-level results from the spindle and PAC analysis and record which sleep episode they came from and its timing relative to Y-maze training/testing. Run:

```bash
python scripts/behavior/run_ymaze.py \
  --input data/raw/ymaze/Ymaze_AnyMaze_all_measures.csv \
  --output data/processed/ymaze \
  --exclusions config/ymaze_exclusions.csv \
  --spindles data/raw/ymaze/spindle_features.csv
```

The script validates uniqueness, joins only on `(mouse_id, condition)`, writes `spindle_join_audit.csv` including unmatched behavior rows, and produces `spindle_correlations.csv` with exploratory Spearman associations separately by genotype **and** condition when at least 5 matched mice have nonconstant values. Smaller groups are marked insufficient rather than producing misleading p-values. No pooled-across-condition pseudo-replication or causal claim.

## Interpretation and follow-up

1. Validate protocol metadata: training duration; delay to test; testing order; when ambient temperature was applied; sleep recording relative to learning; age/sex; video identity for m73; and behavior-to-EEG mouse IDs. Record session dates and counterbalancing in a separate manifest.
2. Check baseline WT novel preference against the 50% reference. Lack of demonstrable learning in controls limits interpretation of genotype deficits or rescue.
3. Primary spindle candidate: NREM slow-wave–sigma PAC preceding REM. Secondary: pre-transition spindle probability (pre-REM, separately pre-wake), spindle density, and sleep stability; prioritize a small number to limit multiple comparisons.
4. Exploratory correlations are **associations**; only properly aligned post-training sleep recording can support a sleep-dependent consolidation hypothesis. Even that cannot alone establish causality.
5. After QC and adequate sample size, use mixed-effects regression or an appropriately paired design for genotype × condition. The current partially repeated dataset is better treated descriptively for condition comparisons.

## Tests

```bash
python -m unittest discover -s tests -p 'test_ymaze_pipeline.py'
```

Never commit animal-level raw data or generated analysis outputs unless approved by the lab; keep them in ignored `data/raw/` and `data/processed/` directories.
# Y-maze: exploration beyond novelty preference

This analysis is **exploratory**, using the 17 test records included after the user-requested QC exclusions: m73 in both conditions, m91 baseline (low tracking), m67 and m89 baseline (zone-time gaps), and m78 ambtemp (zone-time gap). The conditions refer to the behavioral export labels, not to contemporaneous recordings of sleep.

## Scientific question

Could differences in genotype or ambient-temperature condition affect **how** mice explore a novel-arm recognition Y-maze (entries, returns, visit duration, movement) even without establishing a difference in spatial recognition memory?

## Results from ambient-temperature condition

Independent animals: 6 WT versus 3 APP/PS1. Two-sided **exact permutation of genotype labels** was used for mouse-level mean differences (84 possible labelings). Values and p-values below are descriptive/post hoc.

| Outcome | WT mean | APP/PS1 mean | p (uncorrected) |
|---|---:|---:|---:|
| Novel entries / mobile min | 2.53 | 4.59 | 0.0119 |
| Novel entries / meter traveled | 1.05 | 1.86 | 0.0238 |
| All arm entries / mobile min | 7.28 | 11.62 | 0.0238 |
| Mean duration of a novel-arm visit (s) | 11.25 | 6.17 | 0.0357 |
| Number of novel arm entries | 7.33 | 14.00 | 0.0476 |
| Novel entries as share of *all* arm entries (%) | 35.64 | 40.26 | 0.3810 |
| Novel-entry preference (novel / [novel + familiar], %) | 49.37 | 62.48 | 0.2024 |
| Novel time preference (%) | 45.29 | 57.86 | 0.3095 |
| Absolute time in novel arm (s) | 86.5 | 83.37 | 0.9405 |
| Total distance moved (m) | 7.28 | 8.33 | 0.6786 |

The primary picture is **more, briefer visits**, rather than longer overall novel-arm residence. Normalization by mobile time is interesting because each APP/PS1 animal exceeds all WT animals for this endpoint. Yet *all-arm* visit rate is also higher. The novel arm is **not conclusively preferred relative to general arm switching**.

Across a deliberately restricted collection of 17 endpoints in `ambient_genotype_screen.csv`, the smallest Benjamini–Hochberg adjusted q-value is approximately **0.135**; none survives FDR 0.05. These features were selected *after looking at the data*, so even these q-values should **not** be taken as confirmatory multiplicity control across the entire exploration effort. Furthermore, mean visit duration is closely mathematically linked to entry number and total arm time; it is **not an independent biological discovery**.

## Condition effect

Only three genuinely paired WT mice remain (m76, m80, m82). All three show a decline in novel-arm **entry** preference between baseline and ambtemp (changes -27.9, -10.3, -6.2 percentage points). Their **time** preference changes are mixed (+5.0, -12.6, +14.7 points). With n=3, no reliable condition inference is possible; for all three changes in the same direction, a two-sided sign-flip test gives p=0.25. There is only one remaining baseline APP/PS1 mouse and no paired APP/PS1 data. Accordingly, **no genotype-by-temperature effect and no rescue test** can be reliably estimated.

## Interpretation to carry forward

The dataset may reveal an *exploratory strategy / arm-switching* difference under ambient-temperature conditions in APP/PS1 versus WT mice. It does **not** establish memory impairment or rescue. Memory-specific support would require control animals showing valid recognition, sufficient APP/PS1 baseline observations, consistent novel selectivity, and properly matched condition exposures/repeated-test order.

Importantly, zone entry counts depend on AnyMaze ROI boundary definitions and may be inflated by repeated back-and-forth threshold crossings. Validate the detected entries on example videos before treating these measures as behavioral phenotypes. Consider entries with a minimum dwell threshold or requiring center entry between successive arm visits if the tracking time series is available.

## Run

From `sleep-ach-analysis-v2` project root (after running `run_ymaze.py` and activating `ymaze-analysis`):

```bash
python scripts/behavior/explore_ymaze_features.py \
  --input data/processed/ymaze/ymaze_included.csv \
  --output data/processed/ymaze/feature_exploration
```

Outputs: `mouse_level_features.csv`, `ambient_genotype_screen.csv`, `paired_wt_condition_changes.csv` and `genotype_feature_panels.png`.

### Unit tests

```bash
python -m unittest discover -s tests -p 'test_explore_ymaze_features.py' -v
```
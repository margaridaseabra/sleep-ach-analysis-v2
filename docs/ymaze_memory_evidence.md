# Interpreting two-trial Y-maze recognition and potential memory rescue

The behavioral cohort is independent of the EEG/spindle/PAC cohort: cross-cohort mechanistic comparison can only be conceptual, never within-mouse correlation.

## Primary endpoint

The novel arm was blocked during familiarization, and then opened for the recognition test. For the classical protocol, primary outcome: `100*novel_time/(novel_time+familiar_time)`; start arm excluded. The theoretical null reference is 50% *if* arm allocation was counterbalanced and the arms were equally preferred without training.

Use the existing `run_ymaze.py` first, honoring the exclusions in `config/ymaze_exclusions.csv` (m73 both, m91 baseline, m67/m89 baseline, m78 warming). Then run:

```bash
python scripts/behavior/analyze_ymaze_memory_evidence.py \
    --input data/processed/ymaze/ymaze_included.csv \
    --output data/processed/ymaze/memory_evidence
```

Outputs: `ymaze_memory_figure.png`, `recognition_chance_tests.csv`, `paired_test_records.csv`, `locomotion_correlations.csv`, and `ymaze_memory_interpretation.md`.

## Logic of evidence

1. **Assay validity:** WT animals should show credible novelty preference. Also check familiarization exposure, test delay, maze cues, side counterbalancing and repeat tests.
2. **Genotype impairment:** requires a robust baseline WT-vs-APP/PS1 comparison; with only one baseline APP/PS1 mouse it cannot be estimated.
3. **Warming rescue:** requires impairment plus a suitable warming-associated change in APP/PS1, ideally a treatment × genotype interaction with sufficient independent biological replicates and appropriate pairing/controls. A WT-vs-APP comparison within the warming condition alone cannot show rescue.
4. **Locomotion and choice behavior:** report novel-arm time and entry preferences, exploration and distance; avoid treating low novelty preference alone as a QC failure.

Two-sided exact sign-flip tests against the theoretical 50% reference and label-permutation genotype comparison under warming are exploratory; bootstrap intervals in n=3 groups are descriptive and highly unstable. No mass hypothesis testing or covariate-heavy regression is warranted at these sample sizes.

## Evidence-based interpretation

The available 17 test records do not establish a memory deficit, or rescue of such a deficit. This is not equivalent to demonstrating intact memory. The clearest next action is validation of the task and more baseline APP/PS1/WT animals under a protocol that can test genotype and warming effects.
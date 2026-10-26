# Full Y-maze figure, statistics, and Prism export

Reproducible independent analysis of **novel-arm recognition** in the Y-maze. Generates all-memory-outcomes figure, full exploratory-behavior figure, a short exploratory figure, a post-hoc raw AnyMaze screen, statistics, and Prism-ready CSVs. All diagrams compare **WT / APP/PS1** at **baseline / ambient warming** with every retained animal shown.

## Input and QC

Runs on the original `Ymaze_AnyMaze_all_measures.csv` (23 test records) or the full-measure normalized CSV. By default excludes: m73 both conditions (possible duplicate), m91 baseline (<90% tracking), m67 baseline / m89 baseline / m78 ambient (zone-time discrepancies; user-requested). Expected sample sizes after QC: baseline WT 7, baseline APP/PS1 1, ambient WT 6, ambient APP/PS1 3. An optional exclusions CSV adds more restrictions and never overrides the built-in excluded animals.

The data are from a two-trial **novel-arm recognition** test where the novel arm was blocked during familiarization. Do not compute spontaneous alternation from these exports.

## Run on Mac

```bash
cd ~/Downloads/sleep-ach-analysis-v2
conda activate ymaze-analysis
python -m pip install numpy pandas matplotlib

python scripts/behavior/ymaze_full_figure_pipeline.py \
 --input data/raw/ymaze/Ymaze_AnyMaze_all_measures.csv \
 --output data/processed/ymaze/full_figures \
 --exclusions config/ymaze_exclusions.csv

open data/processed/ymaze/full_figures/figures/
open data/processed/ymaze/full_figures/prism/
```

The Python script can also be downloaded and run from any folder with `--input` pointing to the exported AnyMaze CSV. No Windows, AnyMaze installation, video extraction, or supplementary export is required.

## Files produced

- `figures/Figure_Ymaze_memory.png/.pdf`: nine standard time-based, discrimination, entry-based, latency and arm-exploration outcomes.
- `figures/Figure_Ymaze_exploratory_all.png/.pdf`: all 16 selected nonmemory/locomotion and switching outcomes.
- `figures/Figure_Ymaze_exploratory_focus.png/.pdf`: uncorrected p ≤ 0.10 exploratory outcomes plus negative controls; **not confirmatory**.
- `raw_feature_screen/Figure_Ymaze_raw_screen_candidates.png/.pdf`: other original AnyMaze numerical columns with p ≤ 0.10, regardless of biological relevance; these may be redundant or geometry-dependent.
- `statistics/group_descriptives.csv`: mouse count, mean, SD, median and range for every endpoint and condition-by-genotype cell.
- `statistics/genotype_comparisons_by_condition.csv`: exact two-sided permutation tests APP vs WT *within* baseline or warming. Also shows mean difference and exploratory mouse-bootstrap bounds (n≥3 each genotype).
- `statistics/paired_condition_change.csv`: exact paired sign-flip tests for same mouse in both conditions; only three WT pairs remain.
- `statistics/recognition_vs_chance.csv`: two-sided mouse-level exact tests against no-preference 50% (DI=0); baseline APP/PS1 n=1 has no inference.
- `statistics/exploratory_ranked_all_metrics.csv`: selected exploratory panel features with FDR q-values within this defined family.
- `raw_feature_screen/all_anymaze_raw_measure_screen.csv`: ranked broad source-feature screen, with a **separate FDR correction across all tested raw numeric fields**.
- `prism/column_tables/Prism_<metric>.csv`: four columns of independent-looking observations for plotting; **row alignment does not represent within-mouse pairing**.
- `prism/paired_WT/Prism_WT_paired_<metric>.csv`: three validated WT pairs with mouse IDs, baseline and warming columns.
- `prism/Prism_mouse_level_all_metrics.csv` and `Prism_long_all_metrics.csv`: animal IDs for every endpoint.
- `Ymaze_final_included_mouse_level.csv`, `Ymaze_QC_and_exclusions.csv`, `README_results.md`: complete audit and interpretation notes.

## Statistical design and interpretation

The **primary** recognition endpoint is `100 × novel_time/(novel_time+familiar_time)`. The discrimination index is `(novel−familiar)/(novel+familiar)` and is the same information expressed on a −1 to +1 scale. The chance reference for these indices is theoretical and assumes an unbiased test.

Cross-genotype exact permutation tests are run separately within each condition. The baseline APP/PS1 n=1 comparison has no p-value. Condition comparisons display all individual animals, but inferential condition p-values are reported **only for the three genuinely paired WT mice**. **Do not run an independent two-way ANOVA or make a genotype-by-warming interaction/rescue claim** from this design. APP/PS1 has no retained matched pairs.

FDR q-values are exploratory and calculated separately within `secondary_memory`, `exploratory`, and all raw AnyMaze-field families. The broad raw screen is post hoc and should never be presented as an independently validated finding. Significant-looking p-values from these analyses are hypothesis generating, especially given prior feature searching. These mice are separate from the spindle cohort—no direct spindle–memory correlation or mechanism can be inferred.

Tests can be run from the repo root with `python -m unittest discover -s tests -p 'test_ymaze_full_figure_pipeline.py' -v`. The raw-data tests run when the original CSV exists at `data/raw/ymaze/Ymaze_AnyMaze_all_measures.csv`; otherwise those tests are skipped, while tests of the permutation/FDR helpers still run.
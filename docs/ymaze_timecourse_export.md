# Y-maze novel-arm recognition: time-resolved ANY-maze export

This is a **separate exploratory analysis** of recognition over the 300-s test, not a substitute for the full-test primary outcome. The source `.szd` and `.szv` files are proprietary; export the experiment from ANY-maze on Windows. Do not try to rename `.szd` or `.szv` to CSV.

## Export 60-s bins from ANY-maze

1. Back up the original experiment. Open `Ymaze.szd` in ANY-maze; if necessary use `File > Open experiment`. `.szk` files are ANY-maze backups that can be opened using the backup file type.
2. In `Protocol > Analysis across time`, choose **60-second segments** for the 300-second test.
3. Open the **Data** page, then use **Select data** to include the following measures, using the precise zone names configured in the apparatus:
   - Animal ID, Treatment, **Segment of test**
   - Time in the **novel arm** zone, time in the **familiar arm** zone
   - Number of entries to each of these two zones
   - Total distance travelled and Total time mobile (optional)
   - Percentage of frames tracked (optional)
4. With **Segment of test** included, expect 5 rows per 300-s test, labelled `0-60 seconds`, `60-120 seconds`, ... `240-300 seconds`.
5. Save the Data-page spreadsheet as CSV, e.g. `Ymaze_AnyMaze_60s_segments.csv`. If the experiment includes training trials, filter to the **recognition-test stage only**; never combine training and testing.

Alternative: `File > Export > Export test data > CSV` exports raw position/time/zone event data for each test, which can be used for more flexible analyses (including first novel-arm choice and event-resolved habituation) in a *future* raw-data pipeline. This script is for the 60-s Data-page segments, **not** the per-test raw event export.

## Run

```bash
cd ~/Downloads/sleep-ach-analysis-v2
conda activate ymaze-analysis
git pull origin feature/ymaze-novel-arm-recognition
mkdir -p data/raw/ymaze
# Copy the exported CSV into data/raw/ymaze/ first.
python scripts/behavior/analyze_ymaze_timecourse.py \
  --input data/raw/ymaze/Ymaze_AnyMaze_60s_segments.csv \
  --exclusions config/ymaze_exclusions.csv \
  --metadata data/processed/ymaze/ymaze_included.csv \
  --output data/processed/ymaze/timecourse
open data/processed/ymaze/timecourse/
```

## Outputs

- `time_bins_mouse_level.csv`: novel preference and entry preference for each 60-second segment, including QC flags.
- `early_60_120_300s_mouse_level.csv`: per-mouse cumulative first 60, first 120, full 300 s values (incomplete windows are *not* estimated).
- `timecourse_descriptives.csv`: descriptive sample sizes, means and SD by condition, genotype and window.
- `minute_by_minute_preference.png`: individual-mouse preference trajectories by condition.
- `cumulative_preference.png`: cumulative 60/120/300-s preference trajectories.

**Important**: This pipeline is descriptive, and does not claim genotype or treatment effects. Group analysis must use animals (not 60-s bins) as units; the partially repeated design and tiny APP/PS1 baseline sample must be respected. The `config/ymaze_exclusions.csv` exclusions are applied automatically (m73 both conditions, m91 baseline by original tracking QC, m67 and m89 baseline, m78 ambtemp). The `--metadata` file is the validated included dataset; the script requires exactly one genotype for each mouse-condition and never infers genotype from filenames. Recreate `ymaze_included.csv` using `run_ymaze.py` if needed. Look for learned novelty bias at 0–60 s, then assess adaptation; activity confounds and starting-arm location must be considered. A 50% reference assumes novel and familiar arms are equivalent and correctly mapped.

For details on ANY-maze Data page segment exports, see the ANY-maze documentation at https://www.anymaze.co.uk/help/any-maze_data_settings_segment.htm and https://www.anymaze.co.uk/help/any-maze_exporttestdata.htm.
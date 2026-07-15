# Preprocessing decisions

## Raw data immutability

Raw MAT, scoring, temperature, and video files are never edited by the pipeline. Every transformation is written to `data/processed/<session_id>/`.

## Common duration

The usable session duration is the minimum duration across EEG, EMG, and photometry. Only complete one-second epochs are retained:

```text
n_complete_epochs = floor(min(signal durations in seconds))
```

Scoring rows beyond that point are excluded from physiological analyses and counted in QC.

## Scoring representations

Three vectors are kept:

1. `original_code`: imported scoring after duration alignment.
2. `analysis_code`: short Wake bouts are reclassified as microarousals.
3. `nrem_like_code`: microarousals are merged into NREM only for analyses that require uninterrupted NREM-like bouts.

The state codes are:

```text
-1 Undefined
 0 Wake
 1 NREM
 2 REM
15 Microarousal
```

The microarousal definition is a contiguous Wake run strictly shorter than 15 seconds. A 15-second Wake run remains Wake.

## EEG and EMG

The first standardized derivative applies a 50 Hz notch filter and resamples EEG and EMG to 200 Hz. No additional digital bandpass is applied because acquisition hardware already applied the experiment-specific filtering described in the thesis.

## Photometry

The current example MAT-file stores only the processed variable `ne`. The pipeline treats this as an already-processed signal and creates a robust one-second median trace. Raw 465 nm and 405 nm processing cannot be reconstructed from this file.

## Review provenance

All scoring columns, including `final_source`, `review_status`, and `review_notes`, are preserved. Unreviewed epochs trigger a QC warning but do not currently stop preprocessing.

## Temperature

Temperature time is aligned using:

```text
aligned_time_s = Time_s + temperature_offset_s
```

If the temperature logger started 120 seconds before the physiological recording, enter `-120` in the manifest.

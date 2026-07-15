function qc = build_qc_summary(sessionRow, recordingInfo, scoring, maInfo, cfg, temperatureSummary)
%BUILD_QC_SUMMARY Build a one-row table of session-level QC metrics.

if nargin < 7
    temperatureSummary = struct();
end

orig = double(scoring.original_code);
analysis = double(scoring.analysis_code);
n = numel(analysis);

qc = table();
qc.session_id = string(sessionRow.session_id);
qc.pipeline_version = string(cfg.pipeline_version);
qc.condition = string(sessionRow.condition);
qc.mouse_id = string(sessionRow.mouse_id);
qc.genotype = string(sessionRow.genotype);
qc.sex = string(sessionRow.sex);
qc.eeg_input_fs = recordingInfo.eeg_fs;
qc.photometry_input_fs = recordingInfo.photometry_fs;
qc.eeg_output_fs = cfg.target_eeg_fs;
qc.eeg_duration_s = recordingInfo.eeg_duration_s;
qc.emg_duration_s = recordingInfo.emg_duration_s;
qc.photometry_duration_s = recordingInfo.photometry_duration_s;
qc.common_duration_s = recordingInfo.common_duration_s;
qc.max_signal_duration_difference_s = recordingInfo.max_duration_difference_s;
qc.n_complete_epochs = recordingInfo.n_complete_epochs;
qc.scoring_rows_original = scoring.n_rows_original;
qc.scoring_rows_used = scoring.n_rows_used;
qc.scoring_rows_trimmed = scoring.n_rows_trimmed;
qc.reviewed_fraction = scoring.reviewed_fraction;
qc.unreviewed_epochs = scoring.unreviewed_count;
qc.undefined_epochs = sum(orig == cfg.code_undefined);
qc.wake_seconds_original = sum(orig == cfg.code_wake);
qc.nrem_seconds = sum(orig == cfg.code_nrem);
qc.rem_seconds = sum(orig == cfg.code_rem);
qc.wake_seconds_analysis = sum(analysis == cfg.code_wake);
qc.microarousal_seconds = sum(analysis == cfg.code_microarousal);
qc.microarousal_bouts = maInfo.n_microarousal_bouts;
qc.exactly_15s_wake_bouts = maInfo.exact_threshold_wake_bouts;
qc.wake_percent_analysis = 100 * qc.wake_seconds_analysis / n;
qc.nrem_percent = 100 * qc.nrem_seconds / n;
qc.rem_percent = 100 * qc.rem_seconds / n;
qc.microarousal_percent = 100 * qc.microarousal_seconds / n;

qc.temperature_rows = 0;
qc.warming_detected = false;
qc.warming_on_s = NaN;
qc.warming_off_s = NaN;
qc.target_temperature_c = NaN;
qc.active_temperature_mean_c = NaN;
qc.active_temperature_sd_c = NaN;
qc.active_within_0_5c_fraction = NaN;

if ~isempty(fieldnames(temperatureSummary))
    qc.temperature_rows = temperatureSummary.n_rows;
    qc.warming_detected = temperatureSummary.warming_detected;
    qc.warming_on_s = temperatureSummary.warming_on_s;
    qc.warming_off_s = temperatureSummary.warming_off_s;
    qc.target_temperature_c = temperatureSummary.target_temperature_c;
    qc.active_temperature_mean_c = temperatureSummary.active_temperature_mean_c;
    qc.active_temperature_sd_c = temperatureSummary.active_temperature_sd_c;
    qc.active_within_0_5c_fraction = temperatureSummary.active_within_0_5c_fraction;
end
end

function preprocess_session(sessionRow, cfg)
%PREPROCESS_SESSION Run the standardized preprocessing for one manifest row.

sessionId = string(sessionRow.session_id);
outDir = fullfile(cfg.processed_root, char(sessionId));
ensure_dir(outDir);
outMat = fullfile(outDir, 'preprocessed_session.mat');

if isfile(outMat) && ~cfg.overwrite_existing
    error('Output exists for %s. Set cfg.overwrite_existing=true to replace it.', sessionId);
end

recordingPath = resolve_project_path(sessionRow.recording_file, cfg.repo_root);
scoringPath = resolve_project_path(sessionRow.scoring_file, cfg.repo_root);
temperaturePath = resolve_project_path(sessionRow.temperature_file, cfg.repo_root);
photoVariable = string(sessionRow.photometry_variable);

recordingInfo = inspect_recording_mat(recordingPath, photoVariable);
scoring = load_scoring_table(scoringPath, recordingInfo.n_complete_epochs, cfg);

[scoresAnalysis, maInfo] = derive_microarousals(scoring.original_code, cfg.score_epoch_s, cfg);
scoresNremLike = make_nrem_like_scores(scoresAnalysis, cfg);
scoring.analysis_code = scoresAnalysis;
scoring.nrem_like_code = scoresNremLike;

fprintf('  Loading and preprocessing EEG...\n');
D = load(recordingPath, 'eeg');
[eeg200, eegOutputFs, eegDetails] = preprocess_continuous_signal(D.eeg, recordingInfo.eeg_fs, cfg);
clear D;

fprintf('  Loading and preprocessing EMG...\n');
D = load(recordingPath, 'emg');
[emg200, emgOutputFs, emgDetails] = preprocess_continuous_signal(D.emg, recordingInfo.emg_fs, cfg);
clear D;

nTargetSamples = recordingInfo.n_complete_epochs * cfg.target_eeg_fs;
eeg200 = eeg200(1:min(numel(eeg200), nTargetSamples));
emg200 = emg200(1:min(numel(emg200), nTargetSamples));
if numel(eeg200) < nTargetSamples || numel(emg200) < nTargetSamples
    error('Resampled EEG/EMG is shorter than the expected complete-epoch duration.');
end

fprintf('  Loading and reducing photometry to 1 Hz...\n');
photoLoad = load(recordingPath, char(recordingInfo.photometry_variable));
photoRaw = photoLoad.(char(recordingInfo.photometry_variable));
ach1Hz = downsample_median_1hz(photoRaw, recordingInfo.photometry_fs, recordingInfo.n_complete_epochs);
clear photoLoad photoRaw;

scoringOut = scoring.table;
scoringOut.original_code = scoring.original_code;
scoringOut.analysis_code = scoresAnalysis;
scoringOut.nrem_like_code = scoresNremLike;
writetable(scoringOut, fullfile(outDir, 'scores_analysis.csv'));

temperatureAligned = table();
temperatureSummary = struct();
if strlength(temperaturePath) > 0
    offset = value_to_double(sessionRow.temperature_offset_s, 0);
    temperature = load_temperature_log(temperaturePath, offset);
    temperatureAligned = temperature.table;
    temperatureSummary = summarize_temperature(temperatureAligned, cfg);
    writetable(temperatureAligned, fullfile(outDir, 'temperature_aligned.csv'));
end

qc = build_qc_summary(sessionRow, recordingInfo, scoring, maInfo, cfg, temperatureSummary);
writetable(qc, fullfile(outDir, 'qc_summary.csv'));

warningsList = strings(0, 1);
if scoring.n_rows_trimmed > 0
    warningsList(end+1) = sprintf('%d scoring rows were beyond the complete signal duration and were trimmed.', scoring.n_rows_trimmed); %#ok<AGROW>
end
if ~isnan(scoring.reviewed_fraction) && scoring.reviewed_fraction < 1
    warningsList(end+1) = sprintf('Only %.2f%% of aligned scoring epochs are reviewed or auto-filled.', 100 * scoring.reviewed_fraction); %#ok<AGROW>
end
if scoring.undefined_count > 0
    warningsList(end+1) = sprintf('%d aligned scoring epochs are Undefined.', scoring.undefined_count); %#ok<AGROW>
end
if strcmpi(string(sessionRow.genotype), 'UNKNOWN') || strcmpi(string(sessionRow.sex), 'UNKNOWN')
    warningsList(end+1) = 'Genotype and/or sex is UNKNOWN in the manifest.'; %#ok<AGROW>
end
if ~logical(sessionRow.photometry_is_preprocessed)
    warningsList(end+1) = 'Manifest marks photometry as not preprocessed, but this pipeline currently assumes a processed single trace.'; %#ok<AGROW>
end
if strcmpi(string(sessionRow.condition), 'ambtemp') && isempty(fieldnames(temperatureSummary))
    warningsList(end+1) = 'Ambient-warming session has no temperature log.'; %#ok<AGROW>
end
write_qc_report(fullfile(outDir, 'qc_report.md'), qc, warningsList);

metadata = struct();
metadata.session_id = sessionId;
metadata.date = string(sessionRow.date);
metadata.mouse_id = string(sessionRow.mouse_id);
metadata.genotype = string(sessionRow.genotype);
metadata.sex = string(sessionRow.sex);
metadata.age_days = value_to_double(sessionRow.age_days, NaN);
metadata.condition = string(sessionRow.condition);
metadata.replicate = value_to_double(sessionRow.replicate, NaN);
metadata.zt_start = value_to_double(sessionRow.zt_start, NaN);
metadata.recording_file = string(recordingPath);
metadata.scoring_file = string(scoringPath);
metadata.photometry_variable_original = recordingInfo.photometry_variable;
metadata.photometry_is_preprocessed = logical(sessionRow.photometry_is_preprocessed);
metadata.eeg_input_fs = recordingInfo.eeg_fs;
metadata.eeg_output_fs = eegOutputFs;
metadata.emg_output_fs = emgOutputFs;
metadata.photometry_input_fs = recordingInfo.photometry_fs;
metadata.photometry_output_fs = cfg.photometry_output_fs;
metadata.n_complete_epochs = recordingInfo.n_complete_epochs;
metadata.pipeline_version = string(cfg.pipeline_version);
metadata.eeg_preprocessing = eegDetails;
metadata.emg_preprocessing = emgDetails;

scores_original = scoring.original_code; %#ok<NASGU>
scores_analysis = scoresAnalysis; %#ok<NASGU>
scores_nrem_like = scoresNremLike; %#ok<NASGU>
eeg_200 = eeg200; %#ok<NASGU>
emg_200 = emg200; %#ok<NASGU>
ach_1hz = ach1Hz; %#ok<NASGU>

if isempty(temperatureAligned)
    save(outMat, 'eeg_200', 'emg_200', 'ach_1hz', 'scores_original', ...
        'scores_analysis', 'scores_nrem_like', 'metadata', 'qc', cfg.save_mat_version);
else
    save(outMat, 'eeg_200', 'emg_200', 'ach_1hz', 'scores_original', ...
        'scores_analysis', 'scores_nrem_like', 'temperatureAligned', ...
        'metadata', 'qc', cfg.save_mat_version);
end

fprintf('  Saved %s\n', outMat);
end

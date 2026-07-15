function cfg = analysis_defaults(repoRoot)
%ANALYSIS_DEFAULTS Central preprocessing and scoring parameters.

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
end

cfg = struct();
cfg.repo_root = char(repoRoot);
cfg.raw_root = fullfile(cfg.repo_root, 'data', 'raw');
cfg.interim_root = fullfile(cfg.repo_root, 'data', 'interim');
cfg.processed_root = fullfile(cfg.repo_root, 'data', 'processed');

cfg.pipeline_version = '0.1.0';

% Scoring
cfg.score_epoch_s = 1;
cfg.code_undefined = -1;
cfg.code_wake = 0;
cfg.code_nrem = 1;
cfg.code_rem = 2;
cfg.code_microarousal = 15;
cfg.microarousal_threshold_s = 15;
cfg.microarousal_rule = 'strictly_less_than';
cfg.require_all_epochs_reviewed = false;

% EEG/EMG preprocessing
cfg.target_eeg_fs = 200;
cfg.apply_notch = true;
cfg.notch_frequency_hz = 50;
cfg.notch_quality_factor = 35;
cfg.output_signal_type = 'single';

% Photometry
cfg.photometry_output_fs = 1;
cfg.photometry_downsample_method = 'median_per_second';

% Temperature alignment
% aligned_time_s = logger Time_s + temperature_offset_s
cfg.temperature_setpoint_active_threshold_c = 1;

% Output
cfg.overwrite_existing = false;
cfg.save_mat_version = '-v7.3';
end

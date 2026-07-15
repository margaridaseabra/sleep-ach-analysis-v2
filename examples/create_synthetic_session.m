function create_synthetic_session()
%CREATE_SYNTHETIC_SESSION Create a tiny test session under data/raw/demo.

repoRoot = project_setup();
outDir = fullfile(repoRoot, 'data', 'raw', 'demo');
ensure_dir(outDir);

rng(1);
durationS = 120;
eeg_frequency = 610.3515625;
ne_frequency = 1017.2526245;
tEeg = (0:1/eeg_frequency:durationS-1/eeg_frequency)';
tPhoto = (0:1/ne_frequency:durationS-1/ne_frequency)';
eeg = single(0.4*sin(2*pi*2*tEeg) + 0.1*sin(2*pi*50*tEeg) + 0.2*randn(size(tEeg)));
emg = single(0.2*randn(size(tEeg)) + 0.08*sin(2*pi*50*tEeg));
ne = 0.2*sin(2*pi*0.02*tPhoto) + 0.05*randn(size(tPhoto));
num_class = 3;
start_time = 0;
kept_intervals = [0, 1e11];
video_name = '';
video_path = '';
video_start_time = 0;
save(fullfile(outDir, 'demo.mat'), 'eeg', 'emg', 'ne', 'num_class', ...
    'eeg_frequency', 'ne_frequency', 'start_time', 'kept_intervals', ...
    'video_name', 'video_path', 'video_start_time');

code = ones(durationS, 1);
code(1:20) = 0;
code(41:50) = 2;
code(61:70) = 0;      % 10 s microarousal candidate
code(91:105) = 0;     % exactly 15 s, must remain Wake
T = table(strings(durationS,1), (0:durationS-1)', (0:durationS-1)', (1:durationS)', ...
    strings(durationS,1), code, repmat("synthetic", durationS, 1), ...
    repmat("reviewed", durationS, 1), strings(durationS,1), ...
    'VariableNames', {'recording_id','epoch_id','t0_s','t1_s','final_state', ...
    'final_code','final_source','review_status','review_notes'});
writetable(T, fullfile(outDir, 'final_scoring.csv'));

fprintf('Synthetic session written to %s\n', outDir);
end

%% RUN_PREPROCESSING Preprocess every included session in config/sessions.csv.
clear;
clc;

repoRoot = project_setup();
cfg = analysis_defaults(repoRoot);
manifestPath = fullfile(repoRoot, 'config', 'sessions.csv');
manifest = load_manifest(manifestPath);

included = manifest.include_session;
if ~any(included)
    error('No sessions are marked include_session=true in %s.', manifestPath);
end

fprintf('Found %d included session(s).\n', nnz(included));

rows = find(included);
failures = strings(0, 1);
for k = 1:numel(rows)
    idx = rows(k);
    sessionRow = manifest(idx, :);
    sessionId = string(sessionRow.session_id);
    fprintf('\n[%d/%d] Processing %s\n', k, numel(rows), sessionId);
    try
        preprocess_session(sessionRow, cfg);
    catch ME
        warning('Session %s failed: %s', sessionId, ME.message);
        failures(end + 1, 1) = sessionId + ": " + string(ME.message); %#ok<SAGROW>
    end
end

if isempty(failures)
    fprintf('\nAll included sessions completed successfully.\n');
else
    fprintf('\n%d session(s) failed:\n', numel(failures));
    fprintf('  - %s\n', failures);
    error('Preprocessing completed with failures.');
end

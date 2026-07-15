function scoring = load_scoring_table(scoringPath, nCompleteEpochs, cfg)
%LOAD_SCORING_TABLE Import a full QC scoring export or a legacy two-column file.

if ~isfile(scoringPath)
    error('Scoring file not found: %s', scoringPath);
end

T = readtable(scoringPath, 'TextType', 'string');
vars = string(T.Properties.VariableNames);

if ismember("final_code", vars)
    code = double(T.final_code);
elseif ismember("score", vars)
    code = double(T.score);
else
    error('Scoring file must contain final_code or score.');
end

if ismember("t0_s", vars)
    t0 = double(T.t0_s);
elseif ismember("time_s", vars)
    t0 = double(T.time_s);
else
    t0 = (0:height(T)-1)';
    T.t0_s = t0;
end

if ismember("t1_s", vars)
    t1 = double(T.t1_s);
else
    t1 = t0 + cfg.score_epoch_s;
    T.t1_s = t1;
end

if any(abs(diff(t0) - cfg.score_epoch_s) > 1e-6)
    error('Scoring t0_s is not a regular %.3f-second grid.', cfg.score_epoch_s);
end

allowed = [cfg.code_undefined, cfg.code_wake, cfg.code_nrem, cfg.code_rem];
nonMissingCode = code(~isnan(code));
if any(~ismember(nonMissingCode, allowed))
    bad = unique(nonMissingCode(~ismember(nonMissingCode, allowed)));
    error('Scoring contains unsupported code(s): %s', strjoin(string(bad), ', '));
end
code(isnan(code)) = cfg.code_undefined;

nOriginal = height(T);
nUse = min(nOriginal, nCompleteEpochs);
T = T(1:nUse, :);
code = code(1:nUse);

if nUse < nCompleteEpochs
    error('Scoring has %d rows but %d complete signal epochs are available.', nUse, nCompleteEpochs);
end

reviewedFraction = NaN;
unreviewedCount = NaN;
if ismember("review_status", string(T.Properties.VariableNames))
    status = lower(strtrim(string(T.review_status)));
    reviewedMask = status == "reviewed" | status == "auto_filled";
    reviewedFraction = mean(reviewedMask);
    unreviewedCount = sum(~reviewedMask);
end

scoring = struct();
scoring.table = T;
scoring.original_code = int16(code);
scoring.n_rows_original = nOriginal;
scoring.n_rows_used = nUse;
scoring.n_rows_trimmed = max(0, nOriginal - nUse);
scoring.reviewed_fraction = reviewedFraction;
scoring.unreviewed_count = unreviewedCount;
scoring.undefined_count = sum(code == cfg.code_undefined);
end

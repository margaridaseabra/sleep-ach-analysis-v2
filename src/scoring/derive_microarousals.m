function [analysisCode, info] = derive_microarousals(originalCode, epochSeconds, cfg)
%DERIVE_MICROAROUSALS Reclassify Wake runs strictly shorter than 15 seconds.

code = int16(originalCode(:));
analysisCode = code;

isWake = code == cfg.code_wake;
transitions = diff([false; isWake; false]);
starts = find(transitions == 1);
stops = find(transitions == -1) - 1;

durations = (stops - starts + 1) * epochSeconds;
shortMask = durations < cfg.microarousal_threshold_s;

for i = find(shortMask(:))'
    analysisCode(starts(i):stops(i)) = int16(cfg.code_microarousal);
end

info = struct();
info.n_microarousal_bouts = sum(shortMask);
info.microarousal_seconds = sum(durations(shortMask));
info.exact_threshold_wake_bouts = sum(abs(durations - cfg.microarousal_threshold_s) < 1e-9);
info.threshold_s = cfg.microarousal_threshold_s;
info.rule = cfg.microarousal_rule;
end

function summary = summarize_temperature(temperatureTable, cfg)
%SUMMARIZE_TEMPERATURE Compute logger and warming-window QC values.

summary = struct();
summary.n_rows = height(temperatureTable);
summary.start_s = min(double(temperatureTable.aligned_time_s));
summary.end_s = max(double(temperatureTable.aligned_time_s));
summary.duration_s = summary.end_s - summary.start_s;

t = double(temperatureTable.aligned_time_s);
setpoint = double(temperatureTable.Setpoint_C);
temp = double(temperatureTable.Temperature_C);
active = setpoint > cfg.temperature_setpoint_active_threshold_c;

summary.warming_detected = any(active);
summary.warming_on_s = NaN;
summary.warming_off_s = NaN;
summary.target_temperature_c = NaN;
summary.active_temperature_mean_c = NaN;
summary.active_temperature_sd_c = NaN;
summary.active_temperature_min_c = NaN;
summary.active_temperature_max_c = NaN;
summary.active_within_0_5c_fraction = NaN;

if any(active)
    activeIdx = find(active);
    summary.warming_on_s = t(activeIdx(1));
    summary.warming_off_s = t(activeIdx(end));
    summary.target_temperature_c = median(setpoint(active), 'omitnan');
    summary.active_temperature_mean_c = mean(temp(active), 'omitnan');
    summary.active_temperature_sd_c = std(temp(active), 'omitnan');
    summary.active_temperature_min_c = min(temp(active));
    summary.active_temperature_max_c = max(temp(active));
    summary.active_within_0_5c_fraction = mean(abs(temp(active) - summary.target_temperature_c) <= 0.5);
end
end

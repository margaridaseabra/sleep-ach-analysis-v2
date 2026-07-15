function write_qc_report(reportPath, qc, warningsList)
%WRITE_QC_REPORT Write a human-readable Markdown QC report.

fid = fopen(reportPath, 'w');
if fid < 0
    error('Could not open QC report for writing: %s', reportPath);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# Session QC: %s\n\n', char(qc.session_id));
fprintf(fid, '- Condition: `%s`\n', char(qc.condition));
fprintf(fid, '- Mouse: `%s`\n', char(qc.mouse_id));
fprintf(fid, '- Genotype: `%s`\n', char(qc.genotype));
fprintf(fid, '- Sex: `%s`\n', char(qc.sex));
fprintf(fid, '- Complete aligned duration: %.3f s\n', qc.common_duration_s);
fprintf(fid, '- Complete scoring epochs used: %d\n', qc.n_complete_epochs);
fprintf(fid, '- Scoring rows trimmed: %d\n', qc.scoring_rows_trimmed);
fprintf(fid, '- Reviewed fraction: %.4f\n', qc.reviewed_fraction);
fprintf(fid, '- Undefined epochs: %d\n', qc.undefined_epochs);
fprintf(fid, '- Microarousal bouts: %d\n', qc.microarousal_bouts);
fprintf(fid, '- Microarousal seconds: %d\n', qc.microarousal_seconds);
fprintf(fid, '- Exactly 15 s Wake bouts retained as Wake: %d\n', qc.exactly_15s_wake_bouts);

fprintf(fid, '\n## State composition after MA derivation\n\n');
fprintf(fid, '- Wake: %.3f %%\n', qc.wake_percent_analysis);
fprintf(fid, '- NREM: %.3f %%\n', qc.nrem_percent);
fprintf(fid, '- REM: %.3f %%\n', qc.rem_percent);
fprintf(fid, '- Microarousal: %.3f %%\n', qc.microarousal_percent);

if qc.temperature_rows > 0
    fprintf(fid, '\n## Temperature\n\n');
    fprintf(fid, '- Logger rows: %d\n', qc.temperature_rows);
    fprintf(fid, '- Warming detected: %s\n', char(string(qc.warming_detected)));
    fprintf(fid, '- Warming on: %.3f s\n', qc.warming_on_s);
    fprintf(fid, '- Warming off: %.3f s\n', qc.warming_off_s);
    fprintf(fid, '- Target: %.3f C\n', qc.target_temperature_c);
    fprintf(fid, '- Active mean: %.3f C\n', qc.active_temperature_mean_c);
    fprintf(fid, '- Active SD: %.3f C\n', qc.active_temperature_sd_c);
    fprintf(fid, '- Fraction within +/-0.5 C: %.4f\n', qc.active_within_0_5c_fraction);
end

fprintf(fid, '\n## Warnings\n\n');
if isempty(warningsList)
    fprintf(fid, 'No pipeline warnings.\n');
else
    for i = 1:numel(warningsList)
        fprintf(fid, '- %s\n', warningsList(i));
    end
end
end

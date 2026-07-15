function temperature = load_temperature_log(temperaturePath, offsetSeconds)
%LOAD_TEMPERATURE_LOG Read and align an Arduino temperature log.

if ~isfile(temperaturePath)
    error('Temperature file not found: %s', temperaturePath);
end

T = readtable(temperaturePath, 'TextType', 'string');
required = ["Time_s", "Setpoint_C", "Temperature_C"];
missing = setdiff(required, string(T.Properties.VariableNames));
if ~isempty(missing)
    error('Temperature log is missing column(s): %s', strjoin(missing, ', '));
end

T.aligned_time_s = double(T.Time_s) + double(offsetSeconds);
T = sortrows(T, 'aligned_time_s');

temperature = struct();
temperature.table = T;
temperature.offset_s = double(offsetSeconds);
end

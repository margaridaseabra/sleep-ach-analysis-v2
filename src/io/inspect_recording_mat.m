function info = inspect_recording_mat(recordingPath, photometryVariable)
%INSPECT_RECORDING_MAT Inspect variables and signal durations without loading signals.

if ~isfile(recordingPath)
    error('Recording MAT-file not found: %s', recordingPath);
end

vars = whos('-file', recordingPath);
varNames = string({vars.name});
required = ["eeg", "emg", "eeg_frequency", "ne_frequency"];
missing = setdiff(required, varNames);
if ~isempty(missing)
    error('Recording is missing required variable(s): %s', strjoin(missing, ', '));
end

photometryVariable = string(photometryVariable);
if strlength(strtrim(photometryVariable)) == 0
    candidates = ["ne", "ach", "ACh", "dff", "DFF"];
    found = candidates(ismember(candidates, varNames));
    if isempty(found)
        error('No supported photometry variable found.');
    end
    photometryVariable = found(1);
elseif ~ismember(photometryVariable, varNames)
    error('Photometry variable %s was not found.', photometryVariable);
end

meta = load(recordingPath, 'eeg_frequency', 'ne_frequency', 'start_time', ...
    'kept_intervals', 'num_class', 'video_name', 'video_path', 'video_start_time');

eegVar = vars(varNames == "eeg");
emgVar = vars(varNames == "emg");
photoVar = vars(varNames == photometryVariable);

info = struct();
info.path = string(recordingPath);
info.photometry_variable = photometryVariable;
info.eeg_n_samples = prod(eegVar.size);
info.emg_n_samples = prod(emgVar.size);
info.photometry_n_samples = prod(photoVar.size);
info.eeg_fs = double(meta.eeg_frequency);
info.emg_fs = info.eeg_fs;
info.photometry_fs = double(meta.ne_frequency);
info.eeg_duration_s = info.eeg_n_samples / info.eeg_fs;
info.emg_duration_s = info.emg_n_samples / info.emg_fs;
info.photometry_duration_s = info.photometry_n_samples / info.photometry_fs;
info.common_duration_s = min([info.eeg_duration_s, info.emg_duration_s, info.photometry_duration_s]);
info.n_complete_epochs = floor(info.common_duration_s);
info.max_duration_difference_s = max([info.eeg_duration_s, info.emg_duration_s, info.photometry_duration_s]) - info.common_duration_s;
info.metadata = meta;
end

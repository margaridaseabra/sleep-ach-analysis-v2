function [y, outputFs, details] = preprocess_continuous_signal(x, inputFs, cfg)
%PREPROCESS_CONTINUOUS_SIGNAL Apply the standardized notch and resampling.

validateattributes(x, {'numeric'}, {'vector', 'nonempty', 'real', 'finite'});
validateattributes(inputFs, {'numeric'}, {'scalar', 'positive', 'finite'});

x = x(:);
details = struct();
details.input_fs = inputFs;
details.notch_applied = false;

if cfg.apply_notch
    w0 = cfg.notch_frequency_hz / (inputFs / 2);
    if w0 >= 1
        error('Notch frequency %.3f Hz is invalid for input fs %.3f Hz.', cfg.notch_frequency_hz, inputFs);
    end
    bw = w0 / cfg.notch_quality_factor;
    [b, a] = iirnotch(w0, bw);
    x = filtfilt(b, a, x);
    details.notch_applied = true;
end

[p, q] = rat(cfg.target_eeg_fs / inputFs, 1e-12);
y = resample(x, p, q);
outputFs = inputFs * p / q;

details.resample_p = p;
details.resample_q = q;
details.output_fs = outputFs;
details.target_fs = cfg.target_eeg_fs;

switch lower(cfg.output_signal_type)
    case 'single'
        y = single(y);
    case 'double'
        y = double(y);
    otherwise
        error('Unsupported cfg.output_signal_type: %s', cfg.output_signal_type);
end
end

function y = downsample_median_1hz(x, inputFs, nSeconds)
%DOWNSAMPLE_MEDIAN_1HZ Robust per-second median with bounded memory use.

validateattributes(x, {'numeric'}, {'vector', 'nonempty', 'real'});
validateattributes(inputFs, {'numeric'}, {'scalar', 'positive', 'finite'});
validateattributes(nSeconds, {'numeric'}, {'scalar', 'integer', 'positive'});

x = x(:);
y = NaN(nSeconds, 1);
starts = floor((0:nSeconds-1)' * inputFs) + 1;
stops = floor((1:nSeconds)' * inputFs);
stops = min(stops, numel(x));

for second = 1:nSeconds
    firstSample = starts(second);
    lastSample = stops(second);
    if firstSample <= lastSample && firstSample <= numel(x)
        y(second) = median(double(x(firstSample:lastSample)), 'omitnan');
    end
end
end

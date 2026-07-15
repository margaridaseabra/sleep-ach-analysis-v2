function test_downsample_median_1hz()
x = (1:8)';
y = downsample_median_1hz(x, 4, 2);
assert(isequal(y, [2.5; 6.5]));
fprintf('test_downsample_median_1hz passed.\n');
end

function run_all_tests()
%RUN_ALL_TESTS Run lightweight function tests.
project_setup();

tests = {@test_derive_microarousals, @test_downsample_median_1hz, @test_filename_parser};
for i = 1:numel(tests)
    feval(tests{i});
end
fprintf('All %d tests passed.\n', numel(tests));
end

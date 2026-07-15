function test_derive_microarousals()
cfg = analysis_defaults(project_setup());
code = int16([1; zeros(14,1); 1; zeros(15,1); 2]);
[out, info] = derive_microarousals(code, 1, cfg);
assert(all(out(2:15) == 15));
assert(all(out(17:31) == 0));
assert(info.n_microarousal_bouts == 1);
assert(info.microarousal_seconds == 14);
assert(info.exact_threshold_wake_bouts == 1);
fprintf('test_derive_microarousals passed.\n');
end

function test_filename_parser()
meta = parse_session_filename('300526-m63-bas-1.mat');
assert(meta.mouse_id == "m63");
assert(meta.condition == "baseline");
assert(meta.replicate == 1);
assert(year(meta.date) == 2026 && month(meta.date) == 5 && day(meta.date) == 30);
fprintf('test_filename_parser passed.\n');
end

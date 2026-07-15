function meta = parse_session_filename(filename)
%PARSE_SESSION_FILENAME Parse names such as 300526-m63-bas-1.mat.

[~, stem, ~] = fileparts(filename);
expr = '^(?<date>\d{6})-m(?<mouse>\d+)-(?<condition>[A-Za-z0-9_]+)-(?<replicate>\d+)$';
tokens = regexp(stem, expr, 'names', 'once');
if isempty(tokens)
    error('Filename does not match DDMMYY-mXX-condition-replicate: %s', filename);
end

dt = datetime(tokens.date, 'InputFormat', 'ddMMyy');
meta = struct();
meta.date = dt;
meta.mouse_id = "m" + string(tokens.mouse);
meta.condition_raw = string(tokens.condition);
meta.replicate = str2double(tokens.replicate);

switch lower(meta.condition_raw)
    case {'bas', 'baseline'}
        meta.condition = "baseline";
    case {'ambt', 'ambtemp', 'warming', 'warm'}
        meta.condition = "ambtemp";
    case {'drug', 'donepezil', 'chei'}
        meta.condition = "donepezil";
    otherwise
        meta.condition = lower(meta.condition_raw);
end
end

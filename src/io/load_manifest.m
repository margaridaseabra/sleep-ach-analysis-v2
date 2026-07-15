function manifest = load_manifest(manifestPath)
%LOAD_MANIFEST Read and validate the session manifest.

if ~isfile(manifestPath)
    error('Manifest not found: %s', manifestPath);
end

manifest = readtable(manifestPath, 'TextType', 'string');
required = ["session_id", "recording_file", "scoring_file", "condition", "include_session"];
missing = setdiff(required, string(manifest.Properties.VariableNames));
if ~isempty(missing)
    error('Manifest is missing required column(s): %s', strjoin(missing, ', '));
end

manifest.include_session = normalize_logical(manifest.include_session);
if ismember('photometry_is_preprocessed', manifest.Properties.VariableNames)
    manifest.photometry_is_preprocessed = normalize_logical(manifest.photometry_is_preprocessed);
end

if numel(unique(manifest.session_id)) ~= height(manifest)
    error('session_id values must be unique.');
end
end

function values = normalize_logical(values)
if islogical(values)
    return;
end
if isnumeric(values)
    values = values ~= 0;
    return;
end
s = lower(strtrim(string(values)));
validTrue = ismember(s, ["true", "1", "yes", "y"]);
validFalse = ismember(s, ["false", "0", "no", "n", ""]);
if any(~(validTrue | validFalse))
    error('include_session contains an unrecognized value.');
end
values = validTrue;
end

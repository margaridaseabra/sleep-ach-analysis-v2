function resolved = resolve_project_path(pathValue, repoRoot)
%RESOLVE_PROJECT_PATH Resolve a manifest path relative to the repository root.

pathValue = string(pathValue);
if strlength(strtrim(pathValue)) == 0 || ismissing(pathValue)
    resolved = "";
    return;
end

candidate = char(pathValue);
if is_absolute_path(candidate)
    resolved = string(candidate);
else
    resolved = string(fullfile(repoRoot, candidate));
end
end

function tf = is_absolute_path(pathValue)
if ispc
    tf = ~isempty(regexp(pathValue, '^[A-Za-z]:[\\/]', 'once')) || startsWith(pathValue, '\\');
else
    tf = startsWith(pathValue, '/');
end
end

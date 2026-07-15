function repoRoot = project_setup()
%PROJECT_SETUP Add project code and configuration folders to the MATLAB path.

repoRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(repoRoot, 'src')));
addpath(fullfile(repoRoot, 'config'));
addpath(fullfile(repoRoot, 'examples'));
addpath(fullfile(repoRoot, 'tests'));

fprintf('Sleep-ACh project configured at:\n  %s\n', repoRoot);
end

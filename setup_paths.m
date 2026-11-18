function root = setup_paths()
% Add only this package's directories to the current MATLAB session.
root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root,'model'), fullfile(root,'calibration'), ...
    fullfile(root,'sensitivity'), fullfile(root,'simulations'), fullfile(root,'third_party'));
end

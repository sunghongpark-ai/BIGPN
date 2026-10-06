function [Risk, Results] = main(varargin)

OriginalPath = path;
Cleanup = onCleanup(@() path(OriginalPath));
addpath(fullfile(fileparts(mfilename('fullpath')), 'model'));
[Risk, Results] = RunBIGPN(varargin{:});

end

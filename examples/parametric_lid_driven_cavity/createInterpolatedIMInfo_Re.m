function createInterpolatedIMInfo_Re(Re_new, config)
%CREATEINTERPOLATEDIMINFO_RE Interpolate and save IMInfo at Re_new
%   Uses config.interpMethod = 'linear' or 'spline'
%   Reads from IMInfoRDInfoForTraining/<folder>/IMInfo.mat
%   Writes to IMInfoRDInfoInterpolated/<folder>/IMInfo.mat

baseTrainDir = config.baseTrainDir;
Re_all      = config.Re_values(:).';
targetF     = sprintf(config.folderPattern, Re_new);
outDir      = fullfile(config.baseOutputDir, targetF);

% If exact match, just copy
if any(Re_all == Re_new)
    src = fullfile(baseTrainDir, targetF, 'IMInfo.mat');
    if ~exist(src,'file')
        error('Source IMInfo not found: %s', src);
    end
    if ~exist(outDir,'dir'), mkdir(outDir); end
    copyfile(src, fullfile(outDir,'IMInfo.mat'));
    return;
end

% Find bracketing values
Re_lo = max(Re_all(Re_all < Re_new));
Re_hi = min(Re_all(Re_all > Re_new));
if isempty(Re_lo) || isempty(Re_hi)
    error('Re_new = %.3f outside [%g, %g]', Re_new, min(Re_all), max(Re_all));
end
alpha = (Re_new - Re_lo) / (Re_hi - Re_lo);

% Load endpoints
folderLo = fullfile(baseTrainDir, sprintf(config.folderPattern, Re_lo));
folderHi = fullfile(baseTrainDir, sprintf(config.folderPattern, Re_hi));
fileLo   = fullfile(folderLo, 'IMInfo.mat');
fileHi   = fullfile(folderHi, 'IMInfo.mat');
if ~exist(fileLo,'file') || ~exist(fileHi,'file')
    error('Missing IMInfo.mat in %s or %s', folderLo, folderHi);
end
Lo = load(fileLo, 'IMInfo').IMInfo;
Hi = load(fileHi, 'IMInfo').IMInfo;

fields = {
    'parametrization.tangentSpaceAtOrigin'
    'parametrization.nonlinearCoefficients'
    'dual_basis'
};

% Prepare for interp1
switch config.interpMethod
    case 'linear'
        Re_vals = [Re_lo, Re_hi];
        TSdisc   = cat(3, Lo.parametrization.tangentSpaceAtOrigin, ...
                           Hi.parametrization.tangentSpaceAtOrigin);
        Re_query = Re_new;
        method   = 'linear';
    case 'spline'
        Re_vals = Re_all;
        TSdisc   = [];
        for r = Re_vals
            fn = fullfile(baseTrainDir, sprintf(config.folderPattern, r), 'IMInfo.mat');
            D  = load(fn, 'IMInfo').IMInfo.parametrization.tangentSpaceAtOrigin;
            TSdisc = cat(3, TSdisc, D);
        end
        Re_query = Re_new;
        method   = 'spline';
    otherwise
        error('Unknown interpMethod: %s', config.interpMethod);
end

IMInfo = Lo;

for i = 1:numel(fields)
    parts = strsplit(fields{i}, '.');
    % Gather data for interpolation
    dataArr = [];
    for r = Re_vals
        fn = fullfile(baseTrainDir, sprintf(config.folderPattern, r), 'IMInfo.mat');
        A = getfield(load(fn, 'IMInfo').IMInfo, parts{:});
        dataArr = cat(ndims(A)+1, dataArr, A);
    end
    sz      = size(dataArr);
    flat    = reshape(dataArr, [], numel(Re_vals));
    out2d   = interp1(Re_vals, flat.', Re_query, method).';
    B       = reshape(out2d, sz(1:end-1));
    IMInfo = setfield(IMInfo, parts{:}, B);
end

% Save result
V = IMInfo.parametrization.tangentSpaceAtOrigin;
% IMInfo.chart.map = @(x) V' * x;
% IMInfo.chart.map = @(x) IMInfo.dual_basis * x;
dual_basis = IMInfo.dual_basis;
IMInfo.chart.map = @(x) dual_basis * x;
IMInfo.parametrization.map = @(q) [V, IMInfo.parametrization.nonlinearCoefficients] * IMInfo.parametrization.phi(q);
if ~exist(outDir,'dir'), mkdir(outDir); end
save(fullfile(outDir,'IMInfo.mat'), 'IMInfo', '-v7.3');
end
function createInterpolatedRDInfo_Re(Re_new, config)
%CREATEINTERPOLATEDRDINFO_RE Interpolate & save RDInfo at Re_val
%   Uses config.interpMethod = 'linear' or 'spline'
%   config.style = 'polynomial' | 'modal' | 'normalform'

baseTrainDir = config.baseTrainDir;
targetF     = sprintf(config.folderPattern, Re_new);
outDir      = fullfile(config.baseOutputDir, targetF);
Re_all      = config.Re_values(:).';

% --- Find bracket [R0,R1] ---
if Re_new <= Re_all(1)
    R0 = Re_all(1); R1 = Re_all(2);
elseif Re_new >= Re_all(end)
    R0 = Re_all(end-1); R1 = Re_all(end);
else
    k  = find(Re_all <= Re_new, 1, 'last');
    if k == numel(Re_all), k = k-1; end
    R0 = Re_all(k); R1 = Re_all(k+1);
end

% --- Fields to interpolate ---
switch lower(config.style)
    case 'polynomial'
        fields = {
            'reducedDynamics.coefficients'
            'inverseTransformation.coefficients'
            % 'inverseTransformation.lintransf'
            'transformation.coefficients'
            % 'transformation.lintransf'
            'conjugateDynamics.coefficients'
            'eigenvaluesLinPartFlow'
            % 'eigenvectorsLinPart'
        };
    case 'modal'
        fields = {
            'reducedDynamics.coefficients'
            'inverseTransformation.coefficients'
            % 'inverseTransformation.lintransf'
            'transformation.coefficients'
            % 'transformation.lintransf'
            'conjugateDynamics.coefficients'
            'eigenvaluesLinPartFlow'
            'eigenvectorsLinPart'
        };
    % case 'schur'
    %     fields = {
    %         'reducedDynamics.coefficients'
    %         'inverseTransformation.coefficients'
    %         'inverseTransformation.lintransf'
    %         'transformation.coefficients'
    %         'transformation.lintransf'
    %         'conjugateDynamics.coefficients'
    %         'eigenvaluesLinPartFlow'
    %         'eigenvectorsLinPart'
    %     };
    case 'normalform'
        fields = {
            'reducedDynamics.coefficients'
            'inverseTransformation.coefficients'
            % 'inverseTransformation.lintransf'
            'transformation.coefficients'
            % 'transformation.lintransf'
            'conjugateDynamics.coefficients'
            'eigenvaluesLinPartFlow'
            'eigenvectorsLinPart'
        };
    otherwise
        error('Unknown style: %s', config.style);
end

% Load template for interpolation
tmplPath = fullfile(baseTrainDir, sprintf(config.folderPattern, R0), 'RDInfo.mat');
R0_data = load(tmplPath, 'RDInfo').RDInfo;
Rint = R0_data;

% Interpolation setup
switch config.interpMethod
    case 'linear'
        Re_vals  = [R0, R1];
        Re_query = Re_new;
        method   = 'linear';
    case 'spline'
        Re_vals  = Re_all;
        Re_query = Re_new;
        method   = 'spline';
    otherwise
        error('Unknown interpMethod: %s', config.interpMethod);
end

% --- Interpolate each field ---

for i = 1:numel(fields)
    parts = strsplit(fields{i}, '.');
    % gather data
    dataArr = [];
    for r = Re_vals
        R = load(fullfile(baseTrainDir, sprintf(config.folderPattern, r), 'RDInfo.mat'), 'RDInfo').RDInfo;
        A = getfield(R, parts{:});
        dataArr = cat(ndims(A)+1, dataArr, A);
    end
    % reshape & interp1
    sz      = size(dataArr);
    flat    = reshape(dataArr, [], numel(Re_vals));
    if strcmp(fields{i}, 'transformation.coefficients') || ...
       strcmp(fields{i}, 'inverseTransformation.coefficients') || ...
       strcmp(fields{i}, 'eigenvectorsLinPart')
        interpF = interp1(Re_vals, flat.', Re_query, 'nearest').'; % for normalform transformation do nearest
    else
        interpF = interp1(Re_vals, flat.', Re_query, method).';
    end
    % interpF = interp1(Re_vals, flat.', Re_query, method).';
    B       = reshape(interpF, sz(1:end-1));
    % assign back
    Rint = setfield(Rint, parts{:}, B);
end


% --- Model-specific map assignments ---
switch lower(config.style)
    case 'polynomial'
        RDInfo = Rint;
        RDInfo.transformation.map = @(z) RDInfo.transformation.coefficients * RDInfo.transformation.phi(z);
        RDInfo.inverseTransformation.map = @(eta) RDInfo.inverseTransformation.coefficients * RDInfo.inverseTransformation.phi(eta);
        RDInfo.reducedDynamics.map = @(x) RDInfo.reducedDynamics.coefficients * RDInfo.reducedDynamics.phi(x);
        RDInfo.conjugateDynamics.map = @(z) RDInfo.conjugateDynamics.coefficients * RDInfo.conjugateDynamics.phi(z);
    
    % case 'schur'
    %     RDInfo = Rint;
    %     V = RDInfo.eigenvectorsLinPart;
    %     RDInfo.transformation.map = @(y)V*y;
    %     RDInfo.inverseTransformation.map = @(x)V\x;
    %     RDInfo.reducedDynamics.map = @(x) RDInfo.reducedDynamics.coefficients * RDInfo.reducedDynamics.phi(x);
    %     RDInfo.conjugateDynamics.map = @(y) RDInfo.conjugateDynamics.coefficients * RDInfo.reducedDynamics.phi(y);

    case 'modal'
        RDInfo = Rint;
        V = RDInfo.eigenvectorsLinPart;
        % [V,~,~] = eigSorted(RDInfo.reducedDynamics.coefficients(:,1:2));
        % RDInfo.eigenvectorsLinPart = V;
        RDInfo.transformation.map = @(y)V*y;
        RDInfo.inverseTransformation.map = @(x)V\x;
        RDInfo.reducedDynamics.map = @(x) RDInfo.reducedDynamics.coefficients * RDInfo.reducedDynamics.phi(x);
        % RDInfo.conjugateDynamics.map = @(y) V\(RDInfo.reducedDynamics.coefficients * RDInfo.reducedDynamics.phi(V * y));
        RDInfo.conjugateDynamics.map = @(y) RDInfo.conjugateDynamics.coefficients * RDInfo.reducedDynamics.phi(y);

    case 'normalform'
        RDInfo = Rint;
        % conjugate dynamics (map takes [z; conj(z)])
        V = RDInfo.eigenvectorsLinPart;
        % [V,D,~] = eigSorted(RDInfo.reducedDynamics.coefficients(:,1:2));
        % RDInfo.eigenvectorsLinPart = V;
        cdyn = RDInfo.conjugateDynamics;
        RDInfo.conjugateDynamics.map = @(z) transformationComplexConj(cdyn.coefficients * cdyn.phi(transformationComplexConj(z)));
        % --- Forward transformation ---
        T = RDInfo.transformation;
        RDInfo.transformation.map = @(z) V * transformationComplexConj(T.coefficients * T.phi(transformationComplexConj(z)));
        % --- Inverse transformation ---
        iT = RDInfo.inverseTransformation;
        RDInfo.inverseTransformation.map = @(eta) transformationComplexConj(iT.coefficients * iT.phi(V\eta));
        % --- Reduced dynamics ---
        RDInfo.reducedDynamics.map = @(x) RDInfo.reducedDynamics.coefficients * RDInfo.reducedDynamics.phi(x);
end

% --- Save interpolated RDInfo ---
% outDir = fullfile('IMInfoRDInfoInterpolated', sprintf(config.folderPattern, Re_val));
if ~exist(outDir,'dir'), mkdir(outDir); end
save(fullfile(outDir, 'RDInfo.mat'), 'RDInfo');
end
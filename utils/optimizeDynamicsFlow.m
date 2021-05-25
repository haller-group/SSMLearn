function [R,Tinv,N,T,Info,bestOrder,errors] = optimizeDynamicsFlow(etaData, varargin)

p = inputParser;
addParameter(p, 'MaxOrder', 13);
addParameter(p, 'errortype', 'integration', @(x)(strcmp(x,'integration')||strcmp(x,'conjugacy')))
addParameter(p, 'l_vals', [0]);
addParameter(p, 'test_split', 0.2)
parse(p, varargin{:});
l_vals = p.Results.l_vals;

nTraj = size(etaData,1);

if strcmp(p.Results.errortype, 'conjugacy') || p.Results.test_split == 0
    indTest = 1:nTraj;
else
    indTest = round(1:1/p.Results.test_split:nTraj);
end
indTrain = setdiff(1:nTraj, indTest);
if isempty(indTrain)
    indTrain = indTest;
end

ROMOrders = 3:2:p.Results.MaxOrder;
bestResult = Inf;
errors = -ones(size(ROMOrders));

for iOrder = 1:length(ROMOrders)
    ROMOrder = ROMOrders(iOrder);
    [~,Tinvj,Nj,Tj,Infoj] = IMdynamics_flow(etaData(indTrain,:), 'R_PolyOrd', ROMOrder,...
        'style', 'normalform', 'l_vals', l_vals, 'n_folds', 5*(length(l_vals)>1),...
        'fig_disp_nf', 0);
    
    if strcmp(p.Results.errortype, 'integration')
        zData = transformComplex(Tinvj, etaData(indTest,:));
        zRec = integrateFlows(Nj, zData, @(q) q);
        etaRec = transformComplex(Tj, zRec);
        [redDist, ~, ampDist] = computeRecDynErrors(etaRec, etaRec, ...
            etaData(indTest,:), etaData(indTest,:));
        errors(iOrder) = mean(ampDist);
    elseif strcmp(p.Results.errortype, 'conjugacy')
        errors(iOrder) = Infoj.error;
    end
    
    if errors(iOrder) < bestResult
        bestResult = errors(iOrder);
        bestOrder = ROMOrder;
    end
end

[R,Tinv,N,T,Info] = IMdynamics_flow(etaData, 'R_PolyOrd', bestOrder,...
        'style', 'normalform', 'l_vals', l_vals, 'n_folds', 5*(length(l_vals)>1));
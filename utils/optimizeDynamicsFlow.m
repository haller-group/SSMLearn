function [R,Tinv,N,T,Info,bestOrder,errors] = optimizeDynamicsFlow(etaData, varargin)

bestResult = Inf;
ROMOrders = 3:2:13;
% l_vals = [0,1e-4,1e-2,1e0];
l_vals = [0];
errors = -ones(size(ROMOrders));
errortype = 'integration';
% errortype = 'conjugacy';
nTraj = size(etaData,1);
indTest = 1:5:nTraj;
indTrain = setdiff(1:nTraj, indTest);
if isempty(indTrain)
    indTrain = indTest;
end
if strcmp(errortype, 'conjugacy')
    indTrain = 1:nTraj;
end

for iOrder = 1:length(ROMOrders)
    ROMOrder = ROMOrders(iOrder);
    [Rj,Tinvj,Nj,Tj,Infoj] = IMdynamics_flow(etaData(indTrain,:), 'R_PolyOrd', ROMOrder,...
        'style', 'normalform', 'l_vals', l_vals, 'n_folds', 5*(length(l_vals)>1),...
        'fig_disp_nf', 0);
    
    if strcmp(errortype, 'integration')
        zData = transformComplex(Tinvj, etaData(indTest,:));
        zRec = integrateFlows(Nj, zData, @(q) q);
        etaRec = transformComplex(Tj, zRec);
        [redDist, ~, ampDist] = computeRecDynErrors(etaRec, etaRec, ...
            etaData(indTest,:), etaData(indTest,:));
        errors(iOrder) = mean(ampDist);
    elseif strcmp(errortype, 'conjugacy')
        errors(iOrder) = Infoj.error;
    end
    
    if errors(iOrder) < bestResult
        bestResult = errors(iOrder);
        bestOrder = ROMOrder;
    end
end

[R,Tinv,N,T,Info] = IMdynamics_flow(etaData, 'R_PolyOrd', bestOrder,...
        'style', 'normalform', 'l_vals', l_vals, 'n_folds', 5*(length(l_vals)>1));
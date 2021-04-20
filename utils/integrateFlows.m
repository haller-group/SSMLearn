function [etaRec, xRec] = integrateFlows(R, etaData, SSMFunction)

opts = odeset('AbsTol',1e-7);
nTraj = size(etaData,1);
etaRec = cell(nTraj,2); xRec = cell(nTraj,2);
for iTraj = 1:nTraj
    tStart = etaData{iTraj,1}(1);
    tEnd = etaData{iTraj,1}(end);
    nSamp = length(etaData{iTraj,1});
    [t, x] = ode23tb(@(t,y) R(y), linspace(tStart, tEnd, nSamp), etaData{iTraj,2}(:,1), opts);
    etaRec{iTraj,1} = t.';
    etaRec{iTraj,2} = x.';
    xRec{iTraj,1} = t.';
    xRec{iTraj,2} = SSMFunction(etaRec{iTraj,2});
end
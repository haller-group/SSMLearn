function [yRec, xRec] = integrateFlows(R, yData, SSMFunction)

opts = odeset('AbsTol',1e-7);
nTraj = size(yData,1);
yRec = cell(nTraj,2); xRec = cell(nTraj,2);
for iTraj = 1:nTraj
    tStart = yData{iTraj,1}(1);
    tEnd = yData{iTraj,1}(end);
    nSamp = length(yData{iTraj,1});
    [t, x] = ode23tb(@(t,y) R(y), linspace(tStart, tEnd, nSamp), yData{iTraj,2}(:,1), opts);
    yRec{iTraj,1} = t.';
    yRec{iTraj,2} = x.';
    xRec{iTraj,1} = t.';
    xRec{iTraj,2} = SSMFunction(yRec{iTraj,2});
end
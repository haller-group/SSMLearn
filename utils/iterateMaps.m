function [etaRec, yRec] = iterateMaps(R, etaData, SSMFunction)

nTraj = size(etaData,1);
etaRec = cell(nTraj,2); yRec = cell(nTraj,2);
for iTraj = 1:nTraj
    nPoints = length(etaData{iTraj,1});
    etaRec{iTraj,1} = etaData{iTraj,1};
    yRec{iTraj,1} = etaData{iTraj,1};
    etaRec{iTraj,2} = iterate_map(R, nPoints, etaData{iTraj,2}(:,1));
    yRec{iTraj,2} = SSMFunction(etaRec{iTraj,2});
end
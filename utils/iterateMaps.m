function [yRec, xRec] = iterateMaps(R, yData, SSMFunction)

nTraj = size(yData,1);
yRec = cell(nTraj,2); xRec = cell(nTraj,2);
for iTraj = 1:nTraj
    nPoints = length(yData{iTraj,1});
    yRec{iTraj,1} = yData{iTraj,1};
    xRec{iTraj,1} = yData{iTraj,1};
    yRec{iTraj,2} = iterate_map(R, nPoints, yData{iTraj,2}(:,1));
    xRec{iTraj,2} = SSMFunction(yRec{iTraj,2});
end
function [reducedTrajDist, fullTrajDist] = computeRecDynErrors(yRec, xRec, yData, xData)

nTraj = size(yData,1);
reducedTrajDist = zeros(nTraj,1);
fullTrajDist = zeros(nTraj,1); 
for iTraj = 1:nTraj
    nPoints = length(yData{iTraj,1});
    reducedTrajDist(iTraj) = sum(vecnorm(yRec{iTraj,2} - yData{iTraj,2}, 2)) / nPoints / max(vecnorm(yData{iTraj,2},2));
    fullTrajDist(iTraj) = sum(vecnorm(xRec{iTraj,2} - xData{iTraj,2}, 2)) / nPoints / max(vecnorm(xData{iTraj,2},2));
end
function [reducedTrajDist, fullTrajDist] = computeRecDynErrors(yRec, xRec, yData, xData)

nTraj = size(yData,1);
reducedTrajDist = zeros(nTraj,1);
fullTrajDist = zeros(nTraj,1); 
for iTraj = 1:nTraj
    reducedTrajDist(iTraj) = mean(vecnorm(yRec{iTraj,2} - yData{iTraj,2}, 2)) / max(vecnorm(yData{iTraj,2},2));
    fullTrajDist(iTraj) = mean(vecnorm(xRec{iTraj,2} - xData{iTraj,2}, 2)) / max(vecnorm(xData{iTraj,2},2));
end
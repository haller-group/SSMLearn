function [reducedTrajDist, fullTrajDist] = computeRecDynErrors(etaRec, yRec, etaData, yData)

nTraj = size(etaData,1);
reducedTrajDist = zeros(nTraj,1);
fullTrajDist = zeros(nTraj,1); 
for iTraj = 1:nTraj
    reducedTrajDist(iTraj) = mean(vecnorm(etaRec{iTraj,2} - etaData{iTraj,2}, 2)) / max(vecnorm(etaData{iTraj,2},2));
    fullTrajDist(iTraj) = mean(vecnorm(yRec{iTraj,2} - yData{iTraj,2}, 2)) / max(vecnorm(yData{iTraj,2},2));
end
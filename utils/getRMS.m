function RMS = getRMS(xData, SSMFunction, V)
% Compute the mean distance of trajectories in xData and the same
% trajectories projected onto the manifold defined by V and SSMFunction
%
% INPUT
% xData         (n x 2) cell array containing trajectories
% SSMFunction   Anonymous function for the SSM
% V             (d x 2) 2D subspace tangent to the SSM at the fixed point

nTraj = size(xData,1);
meanDist = zeros(nTraj,1);
for iTraj = 1:nTraj
    meanDist(iTraj) = mean(vecnorm(xData{iTraj,2} - SSMFunction(V.'*xData{iTraj,2})));
end
RMS = mean(meanDist);
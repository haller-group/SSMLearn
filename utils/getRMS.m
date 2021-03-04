function RMS = getRMS(xData, SSMFunction, V)
% Compute the mean distance of trajectories in xData and the same
% trajectories projected onto the manifold defined by V and SSMFunction
%
% INPUT
% xData         (n x 1) cell array containing trajectories
% SSMFunction   Anonymous function for the SSM
% V             (d x 2) 2D subspace tangent to the SSM at the fixed point

meanDist = zeros(length(xData),1);
for i = 1:length(xData)
    meanDist(i) = mean(vecnorm(xData{i} - SSMFunction(V'*xData{i})));
end
RMS = mean(meanDist);
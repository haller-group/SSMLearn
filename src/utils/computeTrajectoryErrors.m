function [trajErrors, ampErrors] = computeTrajectoryErrors(yData1, yData2, varargin)
% [trajErrors, ampErrors] = computeTrajectoryErrors(yData1, yData2)
% Compute the normed mean distance of trajectories in yData1 and yData2. The
% total error for each trajectory is normalized by the mean norm of yData2.
% Optionally compute the difference in norm between each pair of points.
% This is useful when the amplitude difference between the trajectory sets
% is more important than the phase difference
%
% INPUT
% yData1       {nTraj x 2}   cell array containing trajectories
% yData2       {nTraj x 2}   cell array containing trajectories
% compinds     (m x 1)       vector of indices to include in comparison
%                               (default 1:end)
%
% OUTPUT
% trajErrors   (nTraj x 1)   vector of distance errors
% ampErrors    (nTraj x 1)   vector of amplitude errors

nTraj = size(yData2,1);
trajErrors = zeros(nTraj,1);
ampErrors = zeros(nTraj,1);
if isempty(varargin)
    inds = 1:size(yData1{1,2},1);
    assert(size(yData1{1,2},1) == size(yData2{1,2},1), 'Dimensionalities of trajectories are inconsistent');
else
    inds = varargin{1};
end
    
for iTraj = 1:nTraj
    if size(yData1{iTraj,2},2) == size(yData2{iTraj,2},2)
        trajErrors(iTraj) = mean(vecnorm(yData1{iTraj,2}(inds,:) - yData2{iTraj,2}(inds,:), 2, 1)) / max(vecnorm(yData2{iTraj,2}(inds,:), 2, 1));
        ampErrors(iTraj) = mean(abs(vecnorm(yData1{iTraj,2}(inds,:), 2, 1) - vecnorm(yData2{iTraj,2}(inds,:), 2, 1))) / mean(vecnorm(yData2{iTraj,2}(inds,:), 2, 1));
    else % integration has failed
        trajErrors(iTraj) = Inf;
        ampErrors(iTraj) = Inf;
    end
end
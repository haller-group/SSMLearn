function outData = sliceTrajectories(inData, interval)
% outData = sliceTrajectories(inData, interval)
% Returns the states in inData that occur between times interval(1) and interval(2).
%
% INPUT
% inData        cell (nTraj x 2)    trajectories
% interval      vector (1 x 2)      time interval
%
% OUTPUT
% outData       cell (nTraj x 2)    trajectories with times in the interval

outData = cell(size(inData));
for iTraj = 1:size(inData,1)
    outData{iTraj,1} = inData{iTraj,1}(:,inData{iTraj,1} >= interval(1) & inData{iTraj,1} <= interval(2));
    outData{iTraj,2} = inData{iTraj,2}(:,inData{iTraj,1} >= interval(1) & inData{iTraj,1} <= interval(2));
end
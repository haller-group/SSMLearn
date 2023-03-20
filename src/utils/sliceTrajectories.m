function outData = sliceTrajectories(inData, interval)
% outData = sliceTrajectories(inData, interval)
% Returns the states in inData that occur between times interval(1) and interval(2).
%
% INPUT
% inData        cell   {nTraj x 2}    trajectories
% interval      vector (1 x 2)        time interval for all trajectories,
%            or matrix (nTraj x 2)    or pass a matrix to apply a different
%                                     interval for each trajectory
% OUTPUT
% outData       cell   {nTraj x 2}    trajectories with times in the
%                                     interval(s)

if size(interval, 1) == 1
    interval = repmat(interval, size(inData, 1), 1);
end
outData = cell(size(inData));
for iTraj = 1:size(inData,1)
    startInd = find(inData{iTraj,1} >= interval(iTraj,1), 1, 'first');
    endInd = find(inData{iTraj,1} <= interval(iTraj,2), 1, 'last');
    for iCol = 1:size(inData,2)
        outData{iTraj,iCol} = inData{iTraj,iCol}(:,startInd:endInd);
    end
end
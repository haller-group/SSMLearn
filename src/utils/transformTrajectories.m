function outData = transformTrajectories(map, inData)
% outData = transformTrajectories(map, inData)
% 
% Apply a map to each trajectory in inData.
% Example: transform reduced coordinates to normal form:
%     zData = transformTrajectories(RDInfo.inverseTransformation.map, etaData);
%   
% INPUT
% map        function handle  function from the space of inData to outData.
% inData     {nTraj x 2}      Cell array of trajectories. First column 
%                             contains time, second column contains state.
%
% OUTPUT
% outData    {nTraj x 2}      Trajectories in transformed coordinates.

nTraj = size(inData,1);
outData = cell(size(inData));
for iTraj = 1:nTraj
    outData{iTraj,1} = inData{iTraj,1};
    outData{iTraj,2} = map(inData{iTraj,2});
end
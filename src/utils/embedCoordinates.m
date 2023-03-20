function yData = embedCoordinates(xData, embeddingDimension, delaySteps)
% yData = embedCoordinates(xData, embeddingDimension, delaySteps)
% Embeds the n-dim. time series xData into a coordinate system of dimension
% p = n*embeddingDimension. 
%
% INPUT
% xData        {nTraj,2}  cell array, the first column contains time 
%                           instances (1 x mi each) and the second column 
%                           the trajectories (n x mi each)
% embeddingDimension int  number of delayed measurements to stack in yData
% delaySteps         int  number of timesteps to delay each measurement
%
%
% OUTPUT
% yData         {nTraj,2}  cell array where the first column contains
%                            time instances (1 x mi each) and the second 
%                            column the trajectories (p x mi each)
% 
% EXAMPLE
% xData{1,1} = [0,0.1,0.2,0.3,0.4];
% xData{1,2} = [1,2,3,4,5;
%               6,7,8,9,10];
% yData = embedCoordinates(xData, 2, 2)
%  yields:
% yData{1,1} = [0,0.1,0.2]
% yData{1,2} = [1 2 3
%               6 7 8
%               3 4 5
%               8 9 10]

yData = cell(size(xData, 1), 2);

for iTraj = 1:size(xData,1)
    yData{iTraj,2} = reshape(xData{iTraj,2}(:, ...
        (0:delaySteps:embeddingDimension*delaySteps-1)'+ ...
        (1:end-(embeddingDimension-1)*delaySteps)), ...
        size(xData{iTraj,2},1)*embeddingDimension, []);
    yData{iTraj,1} = xData{iTraj,1}(1:size(yData{iTraj,2}, 2));
end
function yData = getProjectedTrajs(xData, V)

nTraj = size(xData, 1);
yData = cell(size(xData));
for iTraj = 1:nTraj
    yData{iTraj,1} = xData{iTraj,1};
    yData{iTraj,2} = V.'*xData{iTraj,2};
end
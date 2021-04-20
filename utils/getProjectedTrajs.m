function etaData = getProjectedTrajs(yData, V)

nTraj = size(yData, 1);
etaData = cell(size(yData));
for iTraj = 1:nTraj
    etaData{iTraj,1} = yData{iTraj,1};
    etaData{iTraj,2} = V.'*yData{iTraj,2};
end
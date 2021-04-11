function zData = transformComplex(iT, yData)

nTraj = size(yData,1);
zData = cell(size(yData));
for iTraj = 1:nTraj
    zData{iTraj,1} = yData{iTraj,1};
    zData{iTraj,2} = iT(yData{iTraj,2});
end
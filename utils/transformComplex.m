function zData = transformComplex(iT, etaData)

nTraj = size(etaData,1);
zData = cell(size(etaData));
for iTraj = 1:nTraj
    zData{iTraj,1} = etaData{iTraj,1};
    zData{iTraj,2} = iT(etaData{iTraj,2});
end
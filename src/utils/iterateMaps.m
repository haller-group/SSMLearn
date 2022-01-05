function etaRec = iterateMaps(map, etaData)

nTraj = size(etaData,1);
etaRec = cell(nTraj,2);
for iTraj = 1:nTraj
    nPoints = length(etaData{iTraj,1});
    etaRec{iTraj,1} = etaData{iTraj,1};
    etaRec{iTraj,2} = iterateMap(map, nPoints, etaData{iTraj,2}(:,1));
end
function etaRec = iterateMaps(map, etaData)

nTraj = size(etaData,1);
nCols = size(etaData,2);
etaRec = cell(nTraj,nCols);
if nCols == 2
    for iTraj = 1:nTraj
        nPoints = length(etaData{iTraj,1});
        etaRec{iTraj,1} = etaData{iTraj,1};
        etaRec{iTraj,2} = iterateMap(map, nPoints, etaData{iTraj,2}(:,1));
    end
else
    for iTraj = 1:nTraj
        nPoints = length(etaData{iTraj,1});
        etaRec{iTraj,1} = etaData{iTraj,1}; etaRec{iTraj,3} = etaData{iTraj,3};
        etaRec{iTraj,2} = iterateMap(@(x) map(x,etaData{iTraj,3}),...
            nPoints, etaData{iTraj,2}(:,1));
    end
end
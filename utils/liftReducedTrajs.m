function yLifted = liftReducedTrajs(etaData, SSMFunction)

nTraj = size(etaData,1);
yLifted = cell(nTraj,2);

for iTraj = 1:nTraj
    yLifted{iTraj,1} = etaData{iTraj,1};
    yLifted{iTraj,2} = SSMFunction(etaData{iTraj,2});
end
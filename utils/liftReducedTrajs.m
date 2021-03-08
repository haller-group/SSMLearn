function xLifted = liftReducedTrajs(yData, SSMFunction)

nTraj = size(yData,1);
xLifted = cell(nTraj,2);

for iTraj = 1:nTraj
    xLifted{iTraj,1} = yData{iTraj,1};
    xLifted{iTraj,2} = SSMFunction(yData{iTraj,2});
end
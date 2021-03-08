clearvars
close all

nTraj = 6;
indTest = [1 2];
indTrain = setdiff(1:nTraj, indTest);
SSMDim = 2;

[F, IC] = oscillator(3, nTraj, SSMDim);
% [F, IC] = parabolicSyst(nTraj, 0.8, -0.01, 1, -0.13, [0,0,0], @(t,x) -[0;10*x(1,:).^3]);

observable = @(x) x;
tEnd = 500;
nSamp = 10000;

xSim = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC);

overEmbed = 0;
SSMOrder = 3;

xData = coordinates_embedding(xSim, SSMDim, 'OverEmbedding', overEmbed);

[V, SSMFunction, mfdInfo] = IMparametrization(xData(indTrain,:), SSMDim, SSMOrder, 'c1', 100, 'c2', 0.03);

plotReducedCoords(xData(indTest,2), V);

RRMS = getRMS(xData(indTest,2), SSMFunction, V)

plotReconstructedTrajectory(xData{indTest(1),1}, xData{indTest(1),2}, SSMFunction, V, 2)

plotSSMWithTrajectories(xData(indTrain,2), SSMFunction, [1,3,5], V, 50, 'SSMDimension', SSMDim)
% axis equal
view(50, 30)
%% Reduced dynamics
yData = cell(size(xData));
for iTraj = 1:nTraj
    yData{iTraj,1} = xData{iTraj,1};
    yData{iTraj,2} = V'*xData{iTraj,2};
end
[R,iT,N,T,Maps_info] = IMdynamics_map(yData(indTest,:), 'R_PolyOrd', SSMOrder, 'style', 'modal');

yRec = cell(nTraj,1); xRec = cell(nTraj,1);
reducedTrajDist = zeros(nTraj,1); fullTrajDist = zeros(nTraj,1); 
for iTraj = 1:nTraj
    nPoints = length(yData{iTraj,1});
    yRec{iTraj} = iterate_map(R, nPoints, yData{iTraj,2}(:,1));
    reducedTrajDist(iTraj) = norm(yRec{iTraj} - yData{iTraj,2}) / nPoints;
    xRec{iTraj} = SSMFunction(yRec{iTraj});
    fullTrajDist(iTraj) = norm(xRec{iTraj} - xData{iTraj,2}) / nPoints;
end
RMSE = mean(reducedTrajDist(indTest))
RRMSE = mean(fullTrajDist(indTest))

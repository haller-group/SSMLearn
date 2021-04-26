clearvars
close all

IC = [1e-3, 0.25];
nTraj = size(IC, 2);
indTrain = 2;
indTest = 1;
SSMDim = 2;

F = @(t,R) R-R.^3;

observable = @(x) x;
tEnd = 20;
dt = 0.1;
nSamp = tEnd/dt+1;

xData = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC);
%%
overEmbed = -3;
SSMOrder = 5;

yData = coordinates_embedding(xData, SSMDim, 'OverEmbedding', overEmbed, 'ShiftSteps', 1);

% [V, SSMFunction, mfdInfo] = IMparametrization(yData(indTrain,:), SSMDim, SSMOrder);
V = eye(2);
SSMFunction = @(q)V*q;
%%
% etaData = getProjectedTrajs(yData, V);
etaData = yData;
plotReducedCoords(etaData);

RRMS = getRMS(yData(indTest,:), SSMFunction, V)

% yLifted = liftReducedTrajs(etaData, SSMFunction);
% plotReconstructedTrajectory(yData(indTest(1),:), yLifted(indTest(1),:), 1)

% plotSSMWithTrajectories(yData(indTrain,:), SSMFunction, [1,2,3], V, 15000, 'SSMDimension', SSMDim)
% axis equal
% view(50, 30)
%% Reduced dynamics
ROMOrder = 3;
[R,~,~,~,MapsInfo] = IMdynamics_map(etaData(indTrain,:), 'R_PolyOrd', ROMOrder);

[etaRec, yRec] = iterateMaps(R, etaData, @(eta) SSMFunction(eta));

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(etaRec, yRec, etaData, yData);

RMSE = mean(reducedTrajDist(indTest))
RRMSE = mean(fullTrajDist(indTest))

plotReducedCoords(etaData(indTest(1),:), etaRec(indTest(1),:))
plotReconstructedTrajectory(yData(indTest(1),:), yRec(indTest(1),:), 2, 'm')

reconstructedEigenvalues = computeEigenvaluesMap(MapsInfo, dt)
% DSEigenvalues = lambda(1:SSMDim)


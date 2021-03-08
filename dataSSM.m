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
nSamp = 1000;
dt = tEnd/(nSamp-1);

xSim = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC);

overEmbed = 0;
SSMOrder = 3;

xData = coordinates_embedding(xSim, SSMDim, 'OverEmbedding', overEmbed);

[V, SSMFunction, mfdInfo] = IMparametrization(xData(indTrain,:), SSMDim, SSMOrder, 'c1', 100, 'c2', 0.03);
%%
yData = getProjectedTrajs(xData, V);
plotReducedCoords(yData(indTest,:));

RRMS = getRMS(xData(indTest,:), SSMFunction, V)

xLifted = liftReducedTrajs(yData, SSMFunction);
plotReconstructedTrajectory(xData(indTest(1),:), xLifted(indTest(1),:), 2)

plotSSMWithTrajectories(xData(indTrain,:), SSMFunction, [1,3,5], V, 50, 'SSMDimension', SSMDim)
% axis equal
view(50, 30)
%% Reduced dynamics
[R,iT,N,T,Maps_info] = IMdynamics_map(yData(indTest,:), 'R_PolyOrd', SSMOrder, 'style', 'modal');

[yRec, xRec] = iterateMaps(R, yData, SSMFunction);

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(yRec, xRec, yData, xData);

RMSE = mean(reducedTrajDist(indTest))
RRMSE = mean(fullTrajDist(indTest))

plotReducedCoords(yData(indTest(1),:), yRec(indTest(1),:))
plotReconstructedTrajectory(xData(indTest(1),:), xRec(indTest(1),:), 2)

computeEigenvaluesMap(Maps_info, dt)
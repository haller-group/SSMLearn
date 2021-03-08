clearvars
close all

nTraj = 6;
indTest = [1 2];
indTrain = setdiff(1:nTraj, indTest);

[F, IC] = parabolicSyst(nTraj, 0.8, -0.01, 1, -0.13, [0.5,0,0], @(t,x) -[0;10*x(1,:).^3]);

observable = @(x) x(3,:);
tEnd = 500;
nSamp = 3000;

xSim = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC);

SSMDim = 2;
overEmbed = 100;
SSMOrder = 3;

xData = coordinates_embedding(xSim, SSMDim, 'OverEmbedding', overEmbed);

[V, SSMFunction, mfdInfo] = IMparametrization(xData(indTrain,:), SSMDim, SSMOrder, 'c1', 100, 'c2', 0.03);

yData = getProjectedTrajs(xData, V);
plotReducedCoords(yData(indTest,:));

RRMS = getRMS(xData(indTest,:), SSMFunction, V)

xLifted = liftReducedTrajs(yData, SSMFunction);
plotReconstructedTrajectory(xData(indTest(1),:), xLifted(indTest(1),:), 2)

plotSSMWithTrajectories(xData(indTrain,:), SSMFunction, [1,3,5], V, 50, 'SSMDimension', SSMDim)
% axis equal
view(50, 30)
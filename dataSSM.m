clearvars
close all

nTraj = 5;
nTrajsOnMfd = 5;
indTest = [5];
indTrain = setdiff(1:nTraj, indTest);
SSMDim = 2;
ICRadius = 0.4;

[F, M, C, K, fnl, lambda] = oscillator(3);
ICOnMfd = getSSMIC(M, C, K, fnl, nTrajsOnMfd, ICRadius, SSMDim, 1);
% [F, ICOnMfd, lambda] = parabolicSyst(nTrajsOnMfd, ICRadius, -0.01, 1, -0.13, [0,0,0], @(t,x) -[0;10*x(1,:).^3]);
ICOffMfd = ICRadius * pickPointsOnHypersphere(nTraj-nTrajsOnMfd, 6, rand);
IC = [ICOnMfd, ICOffMfd];

observable = @(x) x(1,:);
tEnd = 500;
nSamp = 2000;
dt = tEnd/(nSamp-1);

xSim = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC);
%%
overEmbed = -2;
SSMOrder = 3;

xData = coordinates_embedding(xSim, SSMDim, 'OverEmbedding', overEmbed);
% xData = coordinates_embedding(xSim, SSMDim, 'ForceEmbedding', 1);

[V, SSMFunction, mfdInfo] = IMparametrization(xData(indTrain,:), SSMDim, SSMOrder, 'c1', 100, 'c2', 0.03);
%%
yData = getProjectedTrajs(xData, V);
plotReducedCoords(yData(indTest,:));

RRMS = getRMS(xData(indTest,:), SSMFunction, V)

xLifted = liftReducedTrajs(yData, SSMFunction);
plotReconstructedTrajectory(xData(indTest(1),:), xLifted(indTest(1),:), 2)

plotSSMWithTrajectories(xData(indTrain,:), SSMFunction, [1,2,3], V, 50, 'SSMDimension', SSMDim)
% axis equal
view(50, 30)
%% Reduced dynamics
[R_modal,iT_modal,N_modal,T_modal,Maps_info_modal] = IMdynamics_map(yData(indTrain,:), 'R_PolyOrd', 3, 'style', 'modal', 'c1', 100, 'c2', 0.03);

[yRecModal, xRecModal] = iterateMaps(R_modal, yData, SSMFunction);

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(yRecModal, xRecModal, yData, xData);

RMSE_modal = mean(reducedTrajDist(indTest))
RRMSE_modal = mean(fullTrajDist(indTest))

plotReducedCoords(yData(indTest(1),:), yRecModal(indTest(1),:))
plotReconstructedTrajectory(xData(indTest(1),:), xRecModal(indTest(1),:), 2, 'g')

reconstructedEigenvalues = computeEigenvaluesMap(Maps_info_modal, dt)
DSEigenvalues = lambda(1:SSMDim)
%% Normal form
[R,iT,N,T,Maps_info] = IMdynamics_map(yData(indTrain,:), 'R_PolyOrd', 7, 'style', 'normalform', 'c1', 0, 'c2', 0.03);

zData = transformToComplex(iT, yData);
[zRec, xRecNormal] = iterateMaps(N, zData, @(q) SSMFunction(T(q)));
yRecNormal = transformToComplex(T, zRec);

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(yRecNormal, xRecNormal, yData, xData);

RMSE_normal = mean(reducedTrajDist(indTest))
RRMSE_normal = mean(fullTrajDist(indTest))

plotReducedCoords(yData(indTest(1),:), yRecNormal(indTest(1),:))
plotReconstructedTrajectory(xData(indTest(1),:), xRecNormal(indTest(1),:), 2, 'c')

reconstructedEigenvalues = computeEigenvaluesMap(Maps_info, dt)
DSEigenvalues = lambda(1:SSMDim)
clearvars
close all

nTraj = 1;
indTest = [1];
indTrain = setdiff(1:nTraj, indTest);
ICRadius = 0.0001;
SSMDim = 2;

nElements = 2;
kappa = 10; % material damping modulus
E = 7000; % Young's modulus (Pa)
[M,C,K,fnl] = von_karman_model(nElements, E, kappa);
n = size(M,1); % mechanical dofs
[IC, mfd, DS, SSM] = getSSMIC(M, C, K, fnl, nTraj, ICRadius, SSMDim, 1)
% noise = 0.3;
% IC = IC .* (1 + noise * (1-2*rand(size(IC)))/2);
% IC = ICRadius * pickPointsOnHypersphere(nTraj, 2*n, 1)

F = @(t,x) DS.odefun(t,x);

% observable = @(x) x(3*nElements-1,:);
observable = @(x) x;
tEnd = 1;
nSamp = 5000;

tic
xSim = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC);
toc
% load vonkarmanfull

plotReconstructedTrajectory(xSim(1,:), xSim(1,:), n)

%%
overEmbed = 0;
SSMOrder = 3;

xData = coordinates_embedding(xSim, SSMDim, 'ForceEmbedding', 1);
% xData = coordinates_embedding(xSim, SSMDim, 'OverEmbedding', overEmbed);

[V, SSMFunction, mfdInfo] = IMparametrization(xData(indTrain,:), SSMDim, SSMOrder, 'c1', 1000, 'c2', 0.1);
%%
yData = getProjectedTrajs(xData, V);
plotReducedCoords(yData(indTest,:));

RRMS = getRMS(xData(indTest,:), SSMFunction, V)
%%
xLifted = liftReducedTrajs(yData, SSMFunction);
plotReconstructedTrajectory(xData(indTest(1),:), xLifted(indTest(1),:), n)
%%
plotSSMWithTrajectories(xData(indTest,:), SSMFunction, [1,3,5], V, 50, 'SSMDimension', SSMDim)
% axis equal
view(50, 30)

%% Reduced dynamics
[R,iT,N,T,Maps_info] = IMdynamics_map(yData(indTrain,:), 'R_PolyOrd', 3, 'style', 'modal', 'c1', 1000, 'c2', 0.1);

[yRec, xRec] = iterateMaps(R, yData, SSMFunction);

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(yRec, xRec, yData, xData);

RMSE = mean(reducedTrajDist(indTest))
RRMSE = mean(fullTrajDist(indTest))

plotReducedCoords(yData(indTest(1),:), yRec(indTest(1),:))
plotReconstructedTrajectory(xData(indTest(1),:), xRec(indTest(1),:), n, 'g')

computeEigenvaluesMap(Maps_info, yRec{1,1}(2)-yRec{1,1}(1))

clearvars
close all

nTraj = 4;
indTest = [1];
indTrain = setdiff(1:nTraj, indTest);
ICRadius = 0.15;
SSMDim = 2;

nElements = 5;
kappa = 4; % cubic spring
gamma = 0; % cubic damping
[M,C,K,fnl] = build_model(kappa, gamma, nElements);
[IC, mfd, DS, SSM] = getSSMIC(M, C, K, fnl, nTraj, ICRadius, SSMDim)

Minv = inv(M);
f = @(q,qdot) [zeros(DS.n-2,1); kappa*q(DS.n-1).^3; 0];

A = [zeros(DS.n), eye(DS.n);
    -Minv*K,     -Minv*C];
G = @(x) [zeros(DS.n,1);
         -Minv*f(x(1:DS.n),x(DS.n+1:2*DS.n))];
F = @(t,x) A*x + G(x);

% F = @(t,x) DS.odefun(t,x);

observable = @(x) x;
tEnd = 100;
nSamp = 15000;

% xSim = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC);
load bernoulli4dfull
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
plotReconstructedTrajectory(xData(indTest(1),:), xLifted(indTest(1),:), 9)
%%
plotSSMWithTrajectories(xData(indTest,:), SSMFunction, [5,9,19], V, 50, 'SSMDimension', SSMDim)
% axis equal
view(50, 30)

%% Reduced dynamics
[R,iT,N,T,Maps_info] = IMdynamics_map(yData(indTest,:), 'R_PolyOrd', 3, 'style', 'modal', 'c1', 1000, 'c2', 0.1);

[yRec, xRec] = iterateMaps(R, yData, SSMFunction);

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(yRec, xRec, yData, xData);

RMSE = mean(reducedTrajDist(indTest))
RRMSE = mean(fullTrajDist(indTest))

plotReducedCoords(yData(indTest(1),:), yRec(indTest(1),:))
plotReconstructedTrajectory(xData(indTest(1),:), xRec(indTest(1),:), 9)

computeEigenvaluesMap(Maps_info, yRec{1,1}(2)-yRec{1,1}(1))

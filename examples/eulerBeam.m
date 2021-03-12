clearvars
close all

nTraj = 4;
indTest = [1];
indTrain = setdiff(1:nTraj, indTest);
ICRadius = 0.4;
SSMDim = 2;

nElements = 5;
kappa = 4; % cubic spring
gamma = 0; % cubic damping
[M,C,K,fnl] = build_model(kappa, gamma, nElements);
n = size(M,1); % mechanical dofs
[IC, mfd, DS, SSM] = getSSMIC(M, C, K, fnl, nTraj, ICRadius, SSMDim, 1)
% noise = 0.3;
% IC = IC .* (1 + noise * (1-2*rand(size(IC)))/2);
% IC = ICRadius * pickPointsOnHypersphere(nTraj, 2*n, 1)

Minv = inv(M);
f = @(q,qdot) [zeros(n-2,1); kappa*q(n-1).^3; 0];

A = [zeros(n), eye(n);
    -Minv*K,     -Minv*C];
G = @(x) [zeros(n,1);
         -Minv*f(x(1:n),x(n+1:2*n))];

F = @(t,x) A*x + G(x);

% F = @(t,x) DS.odefun(t,x);

observable = @(x) x(9,:);
tEnd = 100;
nSamp = 15000;

tic
xSim = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC);
toc
% load bernoulli4dfull
%%
overEmbed = 0;
SSMOrder = 3;

% xData = coordinates_embedding(xSim, SSMDim, 'ForceEmbedding', 1);
xData = coordinates_embedding(xSim, SSMDim, 'OverEmbedding', overEmbed);

[V, SSMFunction, mfdInfo] = IMparametrization(xData(indTrain,:), SSMDim, SSMOrder, 'c1', 1000, 'c2', 0.1);
%%
yData = getProjectedTrajs(xData, V);
plotReducedCoords(yData(indTest,:));

RRMS = getRMS(xData(indTest,:), SSMFunction, V)
%%
xLifted = liftReducedTrajs(yData, SSMFunction);
plotReconstructedTrajectory(xData(indTest(1),:), xLifted(indTest(1),:), 1)
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
plotReconstructedTrajectory(xData(indTest(1),:), xRec(indTest(1),:), 1, 'g')

computeEigenvaluesMap(Maps_info, yRec{1,1}(2)-yRec{1,1}(1))

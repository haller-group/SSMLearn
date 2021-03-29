clearvars
close all

nTraj = 6;
nTrajsOnMfd = 1;
indTest = [1];
indTrain = setdiff(1:nTraj, indTest);
SSMDim = 2;
ICRadius = 0.15;
c1 = 0;
c2 = 0.4;

nElements = 5;
kappa = 4; % cubic spring
gamma = 0; % cubic damping
[M,C,K,fnl] = build_model(kappa, gamma, nElements);
n = size(M,1); % mechanical dofs
[ICOnMfd, mfd, DS, SSM] = getSSMIC(M, C, K, fnl, nTrajsOnMfd, ICRadius, SSMDim, 1);
ICOffMfd = ICRadius * pickPointsOnHypersphere(nTraj-nTrajsOnMfd, 2*n, 1);
IC = [ICOnMfd, ICOffMfd];
lambda = DS.spectrum.Lambda;

Minv = inv(M);
f = @(q,qdot) [zeros(n-2,1); kappa*q(n-1).^3; 0];

A = [zeros(n), eye(n);
    -Minv*K,     -Minv*C];
G = @(x) [zeros(n,1);
         -Minv*f(x(1:n),x(n+1:2*n))];

F = @(t,x) A*x + G(x);

observable = @(x) x(9,:);
tEnd = 100;
nSamp = 15000;
dt = tEnd/(nSamp-1);

tic
xSim = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC);
toc
% load bernoulli4dfull
%%
overEmbed = 0;
SSMOrder = 3;

% xData = coordinates_embedding(xSim, SSMDim, 'ForceEmbedding', 1);
xData = coordinates_embedding(xSim, SSMDim, 'OverEmbedding', overEmbed);

[V, SSMFunction, mfdInfo] = IMparametrization(xData(indTrain,:), SSMDim, SSMOrder, 'c1', c1, 'c2', c2);
% V = [1/sqrt(5)*ones(5,1), 1/sqrt(10)*[-2;-1;0;1;2]];
% SSMFunction = @(q)V*q;
%%
yData = getProjectedTrajs(xData, V);
plotReducedCoords(yData(indTest,:));

RRMS = getRMS(xData(indTest,:), SSMFunction, V)
%%
xLifted = liftReducedTrajs(yData, SSMFunction);
plotReconstructedTrajectory(xData(indTrain(1),:), xLifted(indTrain(1),:), 1)
%%
plotSSMWithTrajectories(xData(indTrain,:), SSMFunction, [1,3,5], V, 50, 'SSMDimension', SSMDim)
% axis equal
view(50, 30)

%% Reduced dynamics
[R,iT,N,T,Maps_info] = IMdynamics_map(sliceTrajectories(yData(indTrain,:), [0.5*tEnd,Inf]), 'R_PolyOrd', 3, 'style', 'modal', 'c1', c1, 'c2', c2);

[yRec, xRec] = iterateMaps(R, yData, SSMFunction);

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(yRec, xRec, yData, xData);

RMSE = mean(reducedTrajDist(indTest))
RRMSE = mean(fullTrajDist(indTest))

plotReducedCoords(yData(indTest(1),:), yRec(indTest(1),:))
plotReconstructedTrajectory(xData(indTest(1),:), xRec(indTest(1),:), 1, 'g')

reconstructedEigenvalues = computeEigenvaluesMap(Maps_info, yRec{1,1}(2)-yRec{1,1}(1))
DSEigenvalues = lambda(1:SSMDim)

%% Normal form
[~,iT,N,T,NormalFormInfo] = IMdynamics_map(sliceTrajectories(yData(indTrain,:), [0.5*tEnd,Inf]), 'R_PolyOrd', 7, 'style', 'normalform', 'c1', c1, 'c2', c2);

zData = transformToComplex(iT, yData);
[zRec, xRecNormal] = iterateMaps(N, zData, @(q) SSMFunction(T(q)));
yRecNormal = transformToComplex(T, zRec);

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(yRecNormal, xRecNormal, yData, xData);

RMSE_normal = mean(reducedTrajDist(indTest))
RRMSE_normal = mean(fullTrajDist(indTest))

plotReducedCoords(yData(indTest(1),:), yRecNormal(indTest(1),:))
plotReconstructedTrajectory(xData(indTest(1),:), xRecNormal(indTest(1),:), 2, 'c')

normalFormEigenvalues = computeEigenvaluesMap(NormalFormInfo, dt)
DSEigenvalues = lambda(1:SSMDim)

%% Backbone curves
N_info = NormalFormInfo.N;
[damp,freq] = nonres_normalform(N_info.coeff,N_info.exponents,dt);
figure
maxRho = abs(zData{indTest(1),2}(1,1));
backbonecurves(damp, freq, SSMFunction, T, 1, maxRho, 'norm');

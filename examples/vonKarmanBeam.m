clearvars
close all

nTraj = 6;
nTrajsOnMfd = 1;
indTest = [1 2];
indTrain = setdiff(1:nTraj, indTest);
ICRadius = 0.001;
SSMDim = 2;

nElements = 5;
E       = 70e9;  % 70e9 % 200e9 % Young's modulus
rho     = 2700; % 2700 % 7850 % density
nu      = 0.3;    % nu
kappa   = 1e8; % material damping modulus 1e8

[M,C,K,fnl] = von_karman_model(nElements, E, rho, nu, kappa);
n = size(M,1); % mechanical dofs (axial def, transverse def, angle)
[ICOnMfd, mfd, DS, SSM] = getSSMIC(M, C, K, fnl, nTrajsOnMfd, ICRadius, SSMDim, 1);
ICOffMfd = ICRadius * pickPointsOnHypersphere(nTraj-nTrajsOnMfd, 2*n, 1);
IC = [ICOnMfd, ICOffMfd];
lambda = DS.spectrum.Lambda;

F = @(t,x) DS.odefun(t,x);

% observable = @(x) x(3*nElements-1,:);
observable = @(x) x;
tEnd = 100;
nSamp = 50000;
sliceInt = [0.4*tEnd, Inf];

tic
% xSim = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC);
load vk5el6traj1stonmfddata.mat
toc
dt = xSim{1,1}(2)-xSim{1,1}(1);

% for iTraj = 1:nTraj
%     xSim{iTraj,2} = xSim{iTraj,2}(n-1,:)
% end

% plotReconstructedTrajectory(xSim(1,:), xSim(1,:), 1)

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
plotReconstructedTrajectory(xData(indTest(2),:), xLifted(indTest(2),:), n-1)
%%
plotSSMWithTrajectories(xData(indTest,:), SSMFunction, [n-4,n-1,2*n-1], V, 50, 'SSMDimension', SSMDim)
% axis equal
view(50, 30)

%% Reduced dynamics
[R,iT,N,T,Maps_info] = IMdynamics_map(sliceTrajectories(yData(indTrain,:), sliceInt), 'R_PolyOrd', 3, 'style', 'modal', 'c1', 0, 'c2', 0.3);

[yRec, xRec] = iterateMaps(R, yData, SSMFunction);

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(yRec, xRec, yData, xData);

RMSE = mean(reducedTrajDist(indTest))
RRMSE = mean(fullTrajDist(indTest))

plotReducedCoords(yData(indTest(1),:), yRec(indTest(1),:))
plotReconstructedTrajectory(xData(indTest(1),:), xRec(indTest(1),:), n-1, 'g')

reconstructedEigenvalues = computeEigenvaluesMap(Maps_info, dt)
DSEigenvalues = lambda(1:SSMDim)

%% Normal form
[~,iT,N,T,NormalFormInfo] = IMdynamics_map(sliceTrajectories(yData(indTrain,:), sliceInt), 'R_PolyOrd', 3, 'style', 'normalform', 'c1', 0, 'c2', 0.03);

zData = transformToComplex(iT, yData);
[zRec, xRecNormal] = iterateMaps(N, zData, @(q) SSMFunction(T(q)));
yRecNormal = transformToComplex(T, zRec);

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(yRecNormal, xRecNormal, yData, xData);

RMSE_normal = mean(reducedTrajDist(indTest))
RRMSE_normal = mean(fullTrajDist(indTest))

plotReducedCoords(yData(indTest(1),:), yRecNormal(indTest(1),:))
plotReconstructedTrajectory(xData(indTest(1),:), xRecNormal(indTest(1),:), n-1, 'c')

normalFormEigenvalues = computeEigenvaluesMap(NormalFormInfo, dt)
DSEigenvalues = lambda(1:SSMDim)

%% Backbone curves
N_info = NormalFormInfo.N;
[damp,freq] = nonres_normalform(N_info.coeff,N_info.exponents,dt);
figure
maxRho = abs(zData{indTest(1),2}(1,1));
backbonecurves(damp, freq, SSMFunction, T, 1, maxRho, 'norm');

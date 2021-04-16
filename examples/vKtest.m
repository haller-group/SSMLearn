clearvars
close all

nTraj = 1;
indTest = [1];
% indTrain = setdiff(1:nTraj, indTest);
indTrain = [1];
ICRadius = 1e-5;
SSMDim = 4;

nElements = 2;
E       = 70e9;   % Young's modulus
rho     = 2700;   % density
nu      = 0.3;    % nu
kappa   = 0e5;    % material damping modulus
l       = 1;
h       = 20e-3;
b       = 50e-3;

[M,C,K,fnl] = von_karman_model(nElements, E, rho, nu, kappa, l, h, b);
n = size(M,1);    % mechanical dofs (axial def, transverse def, angle)
DS = DynamicalSystem();
set(DS, 'M', M, 'C', C, 'K', K, 'fnl', fnl);
F = @(t,x) DS.odefun(t,x);

IC = ICRadius * pickPointsOnHypersphere(nTraj, 2*n, 1);

observable = @(x) x;
tEnd = 15;
nSamp = 50000;
dt = tEnd/(nSamp-1);
tic
xSim = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC);
toc

%%
SSMOrder = 3;
sliceInt = [4, tEnd];

xData = coordinates_embedding(xSim, SSMDim, 'ForceEmbedding', 1);

[V, SSMFunction, mfdInfo] = IMparametrization(sliceTrajectories(xData(indTrain,:), sliceInt), SSMDim, SSMOrder);
%%
yData = getProjectedTrajs(xData, V);
plotReducedCoords(yData(indTest,:));

RRMS = getRMS(xData(indTest,:), SSMFunction, V)
%%
xLifted = liftReducedTrajs(yData, SSMFunction);
plotReconstructedTrajectory(xData(indTest(1),:), xLifted(indTest(1),:), n-1)
%%
plotSSMWithTrajectories(xData(indTest,:), SSMFunction, n-1, V, 50, 'SSMDimension', SSMDim)
% axis equal
view(50, 30)

%% Normal form
xDataTrunc = sliceTrajectories(xData, sliceInt);
yDataTrunc = sliceTrajectories(yData, sliceInt);
[~,iT,N,T,NormalFormInfo] = IMdynamics_map(yDataTrunc, 'R_PolyOrd', 3, 'style', 'normalform');

zData = transformComplex(iT, yDataTrunc);
[zRec, xRecNormal] = iterateMaps(N, zData, @(q) SSMFunction(T(q)));
yRecNormal = transformComplex(T, zRec);

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(yRecNormal, xRecNormal, yDataTrunc, xDataTrunc);

RMSE_normal = mean(reducedTrajDist(indTest))
RRMSE_normal = mean(fullTrajDist(indTest))

plotReducedCoords(yData(indTest(1),:), yRecNormal(indTest(1),:))
plotReconstructedTrajectory(xData(indTest(1),:), xRecNormal(indTest(1),:), n-1, 'c')

SSM(DS).choose_E(1:SSMDim);
lambda = DS.spectrum.Lambda;
DSEigenvalues = lambda(1:SSMDim)
normalFormEigenvalues = computeEigenvaluesMap(NormalFormInfo, dt)

%% Backbone curves
N_info = NormalFormInfo.N;
[damp,freq] = nonres_normalform(N_info.coeff,N_info.exponents,dt);
figure
maxRho = abs(zData{indTest(1),2}(1,1));
backbonecurves(damp, freq, SSMFunction, T, n-1, maxRho, 'norm');

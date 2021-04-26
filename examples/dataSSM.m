clearvars
close all

nTraj = 4;
nTrajsOnMfd = 1;
indTest = [1];
indTrain = setdiff(1:nTraj, indTest);
SSMDim = 2;
ICRadius = 0.4;

[F, M, C, K, fnl, lambda] = oscillator(3, 1, 1, 0.006);
[ICOnMfd, mfd, DS, SSM] = getSSMIC(M, C, K, fnl, nTrajsOnMfd, ICRadius, SSMDim, 60);
% [F, ICOnMfd, lambda] = parabolicSyst(nTrajsOnMfd, ICRadius, -0.01, 1, -0.13, [0,0,0], @(t,x) -[0;10*x(1,:).^3]);
ICOffMfd = ICRadius * pickPointsOnHypersphere(nTraj-nTrajsOnMfd, 6, rand);
IC = [ICOnMfd, ICOffMfd];

observable = @(x) x;
tEnd = 1500;
nSamp = 5e3;
dt = tEnd/(nSamp-1);

xData = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC);
%%
overEmbed = 0;
SSMOrder = 5;

% yData = coordinates_embedding(xSim, SSMDim, 'OverEmbedding', overEmbed, 'ShiftSteps', 1);
yData = coordinates_embedding(xData, SSMDim, 'ForceEmbedding', 1);

sliceInt = [400, tEnd];
yDataTrunc = sliceTrajectories(yData, sliceInt);

[V, SSMFunction, mfdInfo] = IMparametrization(yDataTrunc(indTrain,:), SSMDim, SSMOrder);
%%
etaData = getProjectedTrajs(yData, V);
etaDataTrunc = sliceTrajectories(etaData, sliceInt);
plotReducedCoords(etaData);

RRMS = getRMS(yData(indTest,:), SSMFunction, V)

yLifted = liftReducedTrajs(etaData, SSMFunction);
plotReconstructedTrajectory(yData(indTrain(1),:), yLifted(indTrain(1),:), 2)

plotSSMWithTrajectories(yDataTrunc(indTrain,:), SSMFunction, 1, V, 50, 'SSMDimension', SSMDim)
% axis equal
view(50, 30)
%% Normal form
[~,iT,N,T,NormalFormInfo] = IMdynamics_flow(etaDataTrunc(indTrain,:), 'R_PolyOrd', 3, 'style', 'normalform');

zData = transformComplex(iT, etaDataTrunc);
[zRec, yRec] = integrateFlows(N, zData, @(z) SSMFunction(T(z)));
etaRec = transformComplex(T, zRec);

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(etaRec, yRec, etaDataTrunc, yDataTrunc);

RMSE = mean(reducedTrajDist(indTest))
RRMSE = mean(fullTrajDist(indTest))

plotReducedCoords(etaData(indTest(1),:), etaRec(indTest(1),:))
plotReconstructedTrajectory(yData(indTest(1),:), yRec(indTest(1),:), 2, 'm')

normalFormEigenvalues = computeEigenvaluesMap(NormalFormInfo, dt)
DSEigenvalues = lambda(1:SSMDim)

%% Backbone curves
N_info = NormalFormInfo.N;
[damp,freq] = nonres_normalform(N_info.coeff, N_info.exponents);
figure
maxRho = abs(zData{indTest(1),2}(1,1));
backbonecurves(damp, freq, SSMFunction, T, 1, maxRho, 'norm');

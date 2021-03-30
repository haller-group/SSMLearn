clearvars
close all

load SSM2D
xSim = X_traj;
indTest = [2];
indTrain = [1];
SSMDim = 2;

dt = mean(diff(xSim{1,1}(1,:))); % dt same in all experiments!
tEnd = max(xSim{1,1}(1,:)); % Same in both measurements
%%
overEmbed = -2;
SSMOrder = 1;
shiftSteps = 20;

xData = coordinates_embedding(xSim, SSMDim, 'OverEmbedding', overEmbed, 'ShiftSteps', shiftSteps);

[V, SSMFunction, mfdInfo] = IMparametrization(xData(indTrain,:), SSMDim, SSMOrder, 'c1', 100, 'c2', 0.03);
%%
yData = getProjectedTrajs(xData, V);
plotReducedCoords(yData(indTest,:));

RRMS = getRMS(xData(indTest,:), SSMFunction, V)

xLifted = liftReducedTrajs(yData, SSMFunction);
% plotReconstructedTrajectory(xData(indTest(1),:), xLifted(indTest(1),:), 2)

plotSSMWithTrajectories(xData(indTrain,:), SSMFunction, [1,2,3], V, 50, 'SSMDimension', SSMDim)
% axis equal
view(50, 30)
%% Reduced dynamics
[R,~,~,~,Maps_info] = IMdynamics_flow(yData(indTrain,:), 'R_PolyOrd', 5, 'style', 'modal');

[yRecModal, xRecModal] = integrateFlows(R, yData, SSMFunction);

plotReducedCoords(yData(indTest(1),:), yRecModal(indTest(1),:))
plotReconstructedTrajectory(xData(indTest(1),:), xRecModal(indTest(1),:), 2, 'g')

reconstructedEigenvalues = computeEigenvaluesFlow(Maps_info)
%% Normal form
[~,iT,N,T,NormalFormInfo] = IMdynamics_flow(yData(indTrain,:), 'R_PolyOrd', 3, 'style', 'normalform');

zData = transformToComplex(iT, yData);
[zRec, xRecNormal] = integrateFlows(N, zData, @(q) SSMFunction(T(q)));
yRecNormal = transformToComplex(T, zRec);
% [yRecNormal, xRecNormal] = integrateFlows(R, yData, SSMFunction);

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(yRecNormal, xRecNormal, yData, xData);

RMSE_normal = mean(reducedTrajDist(indTest))
RRMSE_normal = mean(fullTrajDist(indTest))

plotReducedCoords(yData(indTest(1),:), yRecNormal(indTest(1),:))
plotReconstructedTrajectory(xData(indTest(1),:), xRecNormal(indTest(1),:), 2, 'c')

normalFormEigenvalues = computeEigenvaluesFlow(NormalFormInfo)

%% Backbone curves
N_info = NormalFormInfo.N;
[damp,freq] = nonres_normalform(N_info.coeff,N_info.exponents,dt);
figure
maxRho = abs(zData{indTest(1),2}(1,1));
backbonecurves(damp, freq, SSMFunction, T, 1, maxRho, 'norm');

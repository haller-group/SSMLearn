clearvars
close all

IC = [1e-3, -1e-3];
nTraj = size(IC, 2);
indTrain = 1;
indTest = 1;
SSMDim = 1;

mu = 1;
F = @(t,R) mu*R-R.^3;

observable = @(x) x;
tEnd = 20;
dt = 0.1;
nSamp = tEnd/dt+1;

xData = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC);
%%
overEmbed = 0;
SSMOrder = 3;

yData = coordinates_embedding(xData, SSMDim, 'OverEmbedding', overEmbed, 'ShiftSteps', 1);
% yData = xData;

[V, SSMFunction, mfdInfo] = IMparametrization(yData(indTrain,:), SSMDim, SSMOrder);
% V = eye(1);
% SSMFunction = @(q)V*q;
%%
etaData = getProjectedTrajs(yData, V);
% etaData = yData;
% plotReducedCoords(etaData);

RRMS = getRMS(yData(indTest,:), SSMFunction, V)

yLifted = liftReducedTrajs(etaData, SSMFunction);
plotReconstructedTrajectory(yData(indTest(1),:), yLifted(indTest(1),:), 1)

plotSSMWithTrajectories(yData, SSMFunction, [1,2,3], V, 10, 'SSMDimension', SSMDim)
view(50, 30)
%% Reduced dynamics
ROMOrder = 3;
[R,~,~,~,MapsInfo] = IMdynamics_flow(etaData(indTrain,:), 'R_PolyOrd', ROMOrder);
MapsInfo.R.coeff

%%
figure
x = linspace(-1,1.5)';
plot(x, x-x.^3, 'DisplayName', 'True', 'LineWidth', 2)
hold on
plot(x,sum(MapsInfo.R.coeff.*x.^(1:ROMOrder),2), ':', 'DisplayName', ['ROM order ' num2str(ROMOrder)], 'LineWidth', 2)
xlabel('$R$', 'Interpreter', 'latex')
ylabel('$\dot{R}$', 'Interpreter', 'latex')
legend()
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)

figure
plot(xData{1,1}, xData{1,2}, 'LineWidth', 2)
xlabel('$t$ [s]', 'Interpreter', 'latex')
ylabel('$R$', 'Interpreter', 'latex')
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)

%%
[etaRec, yRec] = integrateFlows(R, etaData, @(eta) SSMFunction(eta));

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(etaRec, yRec, etaData, yData);

RMSE = mean(reducedTrajDist(indTest))
RRMSE = mean(fullTrajDist(indTest))

plotReconstructedTrajectory(yData(indTest(1),:), yRec(indTest(1),:), 1, 'm')

reconstructedEigenvalues = computeEigenvaluesFlow(MapsInfo)
DSEigenvalues = mu


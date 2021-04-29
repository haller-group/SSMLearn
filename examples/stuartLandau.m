clearvars
close all

IC = [1e-3, 1e-3];
nTraj = size(IC, 2);
indTrain = 1;
indTest = 2;
SSMDim = 0;

mu = 1;
F = @(t,R) mu*R-R.^3;

observable = @(x) x;
tEnd = 20;
dt = 0.01;
nSamp = tEnd/dt+1;

xData = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC);
%%
overEmbed = 10;
SSMOrder = 15;

if SSMDim > 0
    yData = coordinates_embedding(xData, SSMDim, 'OverEmbedding', overEmbed, 'ShiftSteps', 1);
    nRegVals = 30;
    [V, SSMFunction, mfdInfo] = IMparametrization(yData(indTrain,:), SSMDim, SSMOrder, 'n_folds', 2, 'l_vals', [0, logspace(-6, 6, nRegVals)]);
    
    etaData = getProjectedTrajs(yData, V);
    if(SSMDim>1); plotReducedCoords(etaData); end
    RRMS = getRMS(yData(indTest,:), SSMFunction, V)
    
    yLifted = liftReducedTrajs(etaData, SSMFunction);
    plotReconstructedTrajectory(yData(indTest(1),:), yLifted(indTest(1),:), 1)
    
    plotSSMWithTrajectories(yData, SSMFunction, [1,2,3], V, 10, 'SSMDimension', SSMDim)
    view(50, 30)
else
    yData = xData;
    V = eye(1);
    SSMFunction = @(q)V*q;
    etaData = yData;
end

%% Reduced dynamics
ROMOrder = 13;
nRegVals = 30;
[R,~,~,~,MapsInfo] = IMdynamics_flow(etaData(indTrain,:), 'R_PolyOrd', ROMOrder, 'n_folds', 2, 'l_vals', [0, logspace(-6, 6, nRegVals)]);
disp(MapsInfo.R.coeff)

if SSMDim == 0
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
end

%%
[etaRec, yRec] = integrateFlows(R, etaData, @(eta) SSMFunction(eta));

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(etaRec, yRec, etaData, yData);

RMSE = mean(reducedTrajDist(indTest))
RRMSE = mean(fullTrajDist(indTest))

plotReconstructedTrajectory(yData(indTest(1),:), yRec(indTest(1),:), 1, 'm')

DSEigenvalues = mu
reconstructedEigenvalues = computeEigenvaluesFlow(MapsInfo)

%% DMD
[Phi,omega,lambda,b,Xdmd,t] = DMD(yData{1,2}(:,1:end-1), yData{1,2}(:,2:end), 3, dt);
figure
plot(yData{1,1}, yData{1,2}(1,:), 'DisplayName', 'true', 'LineWidth', 2)
hold on
plot(t, Xdmd(1,:), ':', 'DisplayName', 'DMD', 'LineWidth', 2)
xlabel('$t$ [s]', 'Interpreter', 'latex')
ylabel('$R$', 'Interpreter', 'latex')
legend
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
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
n = size(M,1); % mechanical dofs
[ICOnMfd, mfd, DS, SSM] = getSSMIC(M, C, K, fnl, nTrajsOnMfd, ICRadius, SSMDim, 1)
ICOffMfd = ICRadius * pickPointsOnHypersphere(nTraj-nTrajsOnMfd, 2*n, 1);
IC = [ICOnMfd, ICOffMfd];
lambda = DS.spectrum.Lambda;

F = @(t,x) DS.odefun(t,x);

% observable = @(x) x(3*nElements-1,:);
observable = @(x) x;
tEnd = 100;
nSamp = 50000;

tic
% xSim = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC);
load vk5el6traj1stonmfddata.mat
toc

% for iTraj = 1:nTraj
%     xSim{iTraj,2} = xSim{iTraj,2}(15,:)
% end

% plotReconstructedTrajectory(xSim(1,:), xSim(1,:), n)

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
plotReconstructedTrajectory(xData(indTest(2),:), xLifted(indTest(2),:), n)
%%
plotSSMWithTrajectories(xData(indTest,:), SSMFunction, [n-2,n,2*n], V, 50, 'SSMDimension', SSMDim)
% axis equal
view(50, 30)

%% Reduced dynamics
cutoff = 20000;
for iTraj = indTrain
    yData{iTraj,1} = yData{iTraj,1}(:,cutoff:end);
    xData{iTraj,1} = xData{iTraj,1}(:,cutoff:end);
    yData{iTraj,2} = yData{iTraj,2}(:,cutoff:end);
    xData{iTraj,2} = xData{iTraj,2}(:,cutoff:end);
end
[R,iT,N,T,Maps_info] = IMdynamics_map(yData(indTrain,:), 'R_PolyOrd', 3, 'style', 'modal', 'c1', 0, 'c2', 0.3);

[yRec, xRec] = iterateMaps(R, yData, SSMFunction);

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(yRec, xRec, yData, xData);

RMSE = mean(reducedTrajDist(indTest))
RRMSE = mean(fullTrajDist(indTest))

plotReducedCoords(yData(indTest(1),:), yRec(indTest(1),:))
plotReconstructedTrajectory(xData(indTest(1),:), xRec(indTest(1),:), n, 'g')

reconstructedEigenvalues = computeEigenvaluesMap(Maps_info, yRec{1,1}(2)-yRec{1,1}(1))
DSEigenvalues = lambda(1:SSMDim)
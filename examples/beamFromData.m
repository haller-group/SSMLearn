clearvars
close all

load rspa20160759_si_002.mat
xSim = {NNM1a(1,:)-NNM1a(1,1), NNM1a(2,:);
        NNM1b(1,:)-NNM1b(1,1), NNM1b(2,:);
        NNM2(1,:)-NNM2(1,1), NNM2(2,:);
        NNM3(1,:)-NNM3(1,1), NNM3(2,:)};
indTest = [3];
indTrain = [3];
SSMDim = 2;
c1 = 1000; c2 = 6;

dt = mean(diff(NNM1a(1,:))); % dt same in all experiments!
omega = 45.8*2*pi; zeta = 0.38*0.01;
lambda = -omega*(zeta + [1i; -1i]*sqrt(1-zeta^2));
%%
overEmbed = -2;
SSMOrder = 5;

xData = coordinates_embedding(xSim, SSMDim, 'OverEmbedding', overEmbed);
% xData = coordinates_embedding(xSim, SSMDim, 'ForceEmbedding', 1);

[V, SSMFunction, mfdInfo] = IMparametrization(xData(indTrain,:), SSMDim, SSMOrder, 'c1', c1, 'c2', c2);
% V = [1/sqrt(3)*ones(3,1), 1/sqrt(2)*[-1;0;1]];
% SSMFunction = @(q)V*q;
%%
yData = getProjectedTrajs(xData, V);
plotReducedCoords(yData(indTest,:));

RRMS = getRMS(xData(indTest,:), SSMFunction, V)
%%
xLifted = liftReducedTrajs(yData, SSMFunction);
plotReconstructedTrajectory(xData(indTest,:), xLifted(indTest,:), 1)
%%
plotSSMWithTrajectories(xData(indTest,:), SSMFunction, [1,2,3], V, 50, 'SSMDimension', SSMDim)
% axis equal
view(50, 30)
%% Reduced dynamics
[R,iT,N,T,Maps_info] = IMdynamics_map(yData(indTrain,:), 'R_PolyOrd', 3, 'style', 'modal', 'c1', c1, 'c2', c2);
% [R,iT,N,T,Maps_info] = IMdynamics_flow(yData(indTrain,:), 'R_PolyOrd', 3, 'style', 'modal', 'c1', c1, 'c2', c2);

[yRec, xRec] = iterateMaps(R, yData, SSMFunction);
% yRec = cell(nTraj,2);
% for iTraj = 1:nTraj
% [tt,yy] = ode45(@(t,x) R(x),yData{iTraj,1}(:),yData{iTraj,2}(:,1)');
% yRec{iTraj,1} = tt'; yRec{iTraj,2} = yy';
% end
% xRec = liftReducedTrajs(yRec, SSMFunction);

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(yRec, xRec, yData, xData);

RMSE = mean(reducedTrajDist(indTest))
RRMSE = mean(fullTrajDist(indTest))

plotReducedCoords(yData(indTest(1),:), yRec(indTest(1),:))
plotReconstructedTrajectory(xData(indTest(1),:), xRec(indTest(1),:), 1, 'g')

DSEigenvalues = lambda(1:SSMDim)
reconstructedEigenvalues = computeEigenvaluesMap(Maps_info, dt)
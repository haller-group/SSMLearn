clearvars
close all

nTraj = 5;
indTest = [1];
indTrain = setdiff(1:nTraj, indTest);
ICRadius = 0.25;

nElements = 5;
kappa = 4; % cubic spring
gamma = 0; % cubic damping
[M,C,K,fnl] = build_model(kappa, gamma, nElements);
[IC, mfd, DS, SSM] = getSSMIC(M, C, K, fnl, nTraj, ICRadius)

Minv = inv(M);
f = @(q,qdot) [zeros(DS.n-2,1); kappa*q(DS.n-1).^3; 0];
           
A = [zeros(DS.n), eye(DS.n);
    -Minv*K,     -Minv*C];
G = @(x) [zeros(DS.n,1);
         -Minv*f(x(1:DS.n),x(DS.n+1:2*DS.n))];
F = @(t,x) A*x + G(x);

% F = @(t,x) DS.odefun(t,x);

observable = @(x) x;
tEnd = 500;
nSamp = 5000;

% xSim = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC);
load bernoullidata

SSMDim = 2;
overEmbed = 0;
SSMOrder = 3;

xData = coordinates_embedding(xSim, SSMDim, 'ForceEmbedding', 1);

[V, SSMFunction, mfdInfo] = IMparametrization(xData(indTrain,:), SSMDim, SSMOrder, 'c1', 0, 'c2', 0.03);

plotReducedCoords(xData(indTest,2), V);

RRMS = getRMS(xData(indTest,2), SSMFunction, V)

plotReconstructedTrajectory(xData{indTest(1),1}, xData{indTest(1),2}, SSMFunction, V, 2)

plotSSMWithTrajectories(xData(indTest,2), SSMFunction, [1,17,19], ICRadius)
% axis equal
view(50, 30)
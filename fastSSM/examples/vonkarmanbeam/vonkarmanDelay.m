clearvars
close all
addpath tools

nicered = [0.7,0.1,0.1];
niceblue = [0.1,0.1,0.7];
nicegreen = [0.1,0.9,0.1];
nicegray = [0.6,0.6,0.6];

%% Setup
nElements = 12;
[M, C, K, fnl, fExt, obsdof, PlotFieldonDefMesh] = buildModel(nElements);
obsdof = 0.25*nElements*3-1 % observe transverse displacement at 1/4 of the beam length

n = size(M,1);    % mechanical dofs (axial def, transverse def, angle)
[F, lambda, E] = functionFromTensors(M, C, K, fnl);
[~, lind] = sort(-real(lambda));
lambda = lambda(lind);
E = E(:,lind);
SSMDim = 6;
Ehat = real(E(1:n,1:SSMDim)./E(33,1:SSMDim));
Ehat = Ehat./vecnorm(Ehat);
IC = 0.1*Ehat(1:n,1:2:SSMDim)*[[1;0;0],[0;1;0],[0;0;1],[0.8;-0.8;0.8],[-0.1;0.8;0.8],[-0.6;-0.2;-0.8]];
IC = [IC;zeros(size(IC))];
nTraj = size(IC, 2);
indTrain = [4,5];
indTest = setdiff(1:nTraj, indTrain);

%% Load or simulate decaying data
fname = 'vonkarmandata/dataVKDecaynD.mat'
newMeasurement = false;
observable = @(x) x(obsdof,:);
slowTimeScale = 2*pi/abs(lambda(1));
numberPeriods = 25; numberPointsPerPeriod = 100;
if newMeasurement
    endTime = numberPeriods*slowTimeScale;
    nSamp = numberPeriods*numberPointsPerPeriod+1;
    dt = endTime/(nSamp-1);
    xDataFull = integrateTrajectories(F, endTime, IC, nSamp, 'odetol', 1e-5);
    xDataObs = xDataFull;
    for iTraj = 1:nTraj
        xDataObs{iTraj,2} = observable(xDataFull{iTraj,2});
    end
    DataInfo = struct('nElements', nElements);
    save(fname, 'DataInfo', 'xDataObs', 'xDataFull', 'dt', 'endTime', 'nSamp', 'lambda')
else
    load(fname)
    if nElements ~= DataInfo.nElements
       error('The loaded data comes from a model with a different number of elements.')
    end
end

%% System observation
xData = xDataObs;
embedDim = 50;
outdof = 1;
delaySteps = 1;
tau = dt*delaySteps;
yData = embedCoordinates(xData, embedDim, delaySteps);

sliceInt = [0, endTime];
yDataTrunc = sliceTrajectories(yData, sliceInt);
V = delayTangentSpace(tau, embedDim, lambda(1:SSMDim));
for iTraj = 1:nTraj
    xiData{iTraj,1} = yDataTrunc{iTraj,1};
    xiData{iTraj,2} = V\yDataTrunc{iTraj,2};
end
pt = 4; pm = 1; modeout = 2;
paperFigure('x',['$(V^\dagger y)_',num2str(pm*2-1),'$'],'y',['$(V^\dagger y)_',num2str(pm*2),'$'],'z',['$(V^\dagger y)_',num2str(modeout*2-1),'$']);
plot3(xiData{pt,2}(pm*2-1,:), xiData{pt,2}(pm*2,:), xiData{pt,2}(modeout*2-1,:), 'linewidth', 2, 'displayname', ['Projection of traj. ',num2str(pt), ' onto $V$'], 'color', nicered)
view(3); xticklabels([]); yticklabels([]); zticklabels([]); xticks([]); yticks([]); zticks([])

%% fastSSM modeling

SSMOrder = 3;
ROMOrder = 3;
NFOrder = 7;

[Mmap, iMmap, Tmap, iTmap, Nflow] = fastSSMplus(yDataTrunc(indTrain,:), SSMDim, SSMOrder, ROMOrder, NFOrder, 'mini', V);

clear yRec
for iTraj = 1:nTraj
    [t,zReci] = ode15s(Nflow, yDataTrunc{iTraj,1}, iTmap(iMmap(yDataTrunc{iTraj,2}(:,1))));
    yRec{iTraj,1} = t.';
    yRec{iTraj,2} = Mmap(Tmap(zReci.'));
    zRec{iTraj,1} = t.';
    zRec{iTraj,2} = zReci.';
end
NMTE = computeTrajectoryErrors(yRec, yDataTrunc)

%% Plot prediction
pt = 5;
paperFigure();
plot(yData{pt,1}(1,:), yData{pt,2}(outdof,:), 'Color', [0,0,0], 'LineWidth', 2.6, 'DisplayName', 'Training data');
plot(yRec{pt,1}(1,:), yRec{pt,2}(outdof,:), '--', 'Color', nicegreen, 'LineWidth', 2.6, 'DisplayName', 'Reconstruction');
xlabel('time [s]', 'Interpreter', 'latex');
ylabel('$w$ [m]', 'Interpreter', 'latex');
xlim([0,0.2])
legend

%% Plot trajectories in observable space
pdelays = [0,8,16];
scale = [1,1,2];
pind = 1+pdelays;
colors = colororder;
paperFigure('x','$w(t)$','y',['$w(t+',num2str(pdelays(2)),'\Delta t)$'],'z',['$w(t+',num2str(pdelays(3)),'\Delta t)$'],'legendcols',1);
% plot3(scale(1)*yDataTrunc{1,2}(pind(1),:), scale(1)*yDataTrunc{1,2}(pind(2),:), scale(1)*yDataTrunc{1,2}(pind(3),:), 'linewidth', 2, 'displayname','Mode 1','color',colors(2,:))
% plot3(scale(2)*yDataTrunc{2,2}(pind(1),:), scale(2)*yDataTrunc{2,2}(pind(2),:), scale(2)*yDataTrunc{2,2}(pind(3),:), 'linewidth', 2, 'displayname','Mode 2','color',colors(5,:))
% plot3(scale(3)*yDataTrunc{3,2}(pind(1),:), scale(3)*yDataTrunc{3,2}(pind(2),:), scale(3)*yDataTrunc{3,2}(pind(3),:), 'linewidth', 2, 'displayname',['Mode 3 ($\times$',num2str(scale(3)),')'],'color',colors(1,:))
plot3(yDataTrunc{4,2}(pind(1),:), yDataTrunc{4,2}(pind(2),:), yDataTrunc{4,2}(pind(3),:), 'linewidth', 2, 'displayname', 'Traj. 4', 'color', nicegray)
vscale = 1.5;
v1 = vscale*max(vecnorm(yData{1,2}(pind,:)))*V(pind,1)/norm(V(pind,1));
v2 = vscale*max(vecnorm(yData{1,2}(pind,:)))*V(pind,2)/norm(V(pind,2));
v3 = vscale*max(vecnorm(yData{1,2}(pind,:)))*V(pind,3)/norm(V(pind,3));
v4 = vscale*max(vecnorm(yData{1,2}(pind,:)))*V(pind,4)/norm(V(pind,4));
v5 = vscale*max(vecnorm(yData{1,2}(pind,:)))*V(pind,5)/norm(V(pind,5));
v6 = vscale*max(vecnorm(yData{1,2}(pind,:)))*V(pind,6)/norm(V(pind,6));
quiver3(0,0,0,v1(1),v1(2),v1(3), 'linewidth', 5, 'MaxHeadSize', 1, 'autoscalefactor', 1, 'displayname', '$\mathrm{Re}\,V_1$', 'color', colors(2,:))
% quiver3(0,0,0,v2(1),v2(2),v2(3), 'linewidth', 3, 'MaxHeadSize', 1, 'autoscalefactor', 1, 'displayname', '$\mathrm{Im}\,V_1$')
quiver3(0,0,0,v3(1),v3(2),v3(3), 'linewidth', 5, 'MaxHeadSize', 1, 'autoscalefactor', 1, 'displayname', '$\mathrm{Re}\,V_3$', 'color', colors(5,:))
% quiver3(0,0,0,v4(1),v4(2),v4(3), 'linewidth', 3, 'MaxHeadSize', 1, 'autoscalefactor', 1, 'displayname', '$\mathrm{Im}\,V_3$')
quiver3(0,0,0,v5(1),v5(2),v5(3), 'linewidth', 5, 'MaxHeadSize', 1, 'autoscalefactor', 1, 'displayname', '$\mathrm{Re}\,V_5$', 'color', colors(1,:))
% quiver3(0,0,0,v6(1),v6(2),v6(3), 'linewidth', 3, 'MaxHeadSize', 1, 'autoscalefactor', 1, 'displayname', '$\mathrm{Im}\,V_5$')
xticklabels([]); yticklabels([]); zticklabels([])
view(-134,16)

%% Transform data to normal form for visualization
for iTraj = 1:nTraj
    zData{iTraj,1} = yData{iTraj,1};
    disp(['transforming trajectory ', num2str(iTraj), ' out of ', num2str(nTraj), ' to normal form for visualization...'])
    for iN = 1:size(zRec{iTraj,1},2)
        zData{iTraj,2}(:,iN) = iTmap(iMmap(yData{iTraj,2}(:,iN)));
    end
end

%% Plot phase portrait in normal form
paperFigure('x','$\rho_1$','y','$\rho_2$','z','$\rho_3$','legendcols',3);
for iTraj = 1:nTraj
    rho = abs(zData{iTraj,2});
    plot3(rho(1,:),rho(3,:),rho(5,:), 'linewidth', 4, 'DisplayName', ['Traj. ', num2str(iTraj)])
end
view(45,30)

%% Visualize normal form
pmodes = [1,2];
modeout = 1;
% paperFigure('x',['$\rho_',num2str(pmodes(1)),'$'],'y',['$\rho_',num2str(pmodes(2)),'$'],'z',['$\dot{\theta}_',num2str(modeout),'$'],'legendcols',0);
paperFigure('x',['$\rho_',num2str(pmodes(1)),'$'],'y',['$\rho_',num2str(pmodes(2)),'$'],'z',['$\dot{\rho}_',num2str(modeout),'/\rho_',num2str(modeout),'$'],'legendcols',0);
for iTraj = 1:nTraj
    z = zData{iTraj,2};
    rho = abs(z);
    zDot = Nflow(zData{iTraj,1}, z);
%     dot = imag(zDot./z);
    dot = real(zDot./z);
    plot3(rho(pmodes(1)*2-1,:),rho(pmodes(2)*2-1,:),dot(2*modeout-1,:), 'linewidth', 3, 'DisplayName', ['Trajectory ', num2str(iTraj)])
end
[X,Y] = meshgrid(1e-6:0.02:1,1e-6:0.02:1);
XY = 1e-6*ones(SSMDim,length(X(:)));
XY((pmodes(1)*2-1):(pmodes(1)*2),:) = [X(:),X(:)]';
XY((pmodes(2)*2-1):(pmodes(2)*2),:) = [Y(:),Y(:)]';
% Nsurf = imag(Nflow(0,XY)./XY);
Nsurf = real(Nflow(0,XY)./XY);
Nsurf = reshape(Nsurf(2*modeout-1,:)',size(X));
surf(X, Y, Nsurf, 'facecolor', [0,0,0], 'EdgeColor', [0,0,0], 'facealpha', 0.1, 'edgealpha', 0.3)
view(3)
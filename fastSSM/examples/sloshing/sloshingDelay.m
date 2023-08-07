%% Finding an nD SSM from delay-embedded sloshing data
% 

clearvars
close all
clc

%% Example setup
load highampdata.mat % note: to save repository space, this data contains only every 10th point along the surface. The changes to the resulting model are negligible. The full dataset is available upon request.
nTraj = size(xData,1);
width = 500; % tank width in mm
modes = [1, 2, 4];
dt = 0.01;
nPoints = size(xData{1,2}, 1);
nicered = [0.7,0.1,0.1];
niceblue = [0.1,0.1,0.7];
nicegreen = [0.1,0.9,0.1];
nicegray = [0.6,0.6,0.6];

%% Delay embedding
shiftsteps = 5;
embedDim = 47;
outdof = 1; % leftmost point

yData = embedCoordinates(xData, embedDim, shiftsteps);

c = [-0.05,-0.07,-0.08,-0.09,-0.1];
omeg = sqrt(9.81*pi/(width/1000)*[1:5].*tanh(pi*(1:5)*400/width))
modesAnalysis = 1:5;
c = c(modesAnalysis); omeg = omeg(modesAnalysis);
lambda = reshape([c;c], 1, []) + 1i*reshape([omeg;-omeg], 1, []);

sliceInt = [1.2, Inf];
yDataTrunc = sliceTrajectories(yData, sliceInt);
clear yData

%% Tangent space
SSMDim = 2*length(modes);
SSMOrder = 4;
ROMOrder = 3;
NFOrder = 3;

W = cos((1/2+linspace(-1/2, 1/2, nPoints)')*modesAnalysis*pi);

V = delayTangentSpace(shiftsteps*dt, embedDim, lambda)

clear T
for mode = modesAnalysis
    T(:,2*mode-1) = reshape(W(:,mode)*V(:,2*mode-1).', [], 1);
    T(:,2*mode) = reshape(W(:,mode)*V(:,2*mode).', [], 1);
end
T = T./vecnorm(T);

%% Modal coordinates
xiData = transformTrajectories(@(y)pinv(T)*y, yDataTrunc);

%% Plot
pt = 2
pm = 1
modeout = 3
paperFigure('x',['$(T^\dagger y)_',num2str(pm*2-1),'$'],'y',['$(T^\dagger y)_',num2str(pm*2),'$'],'z',['$(T^\dagger y)_',num2str(modeout*2-1),'$']);
plot3(xiData{pt,2}(pm*2-1,:), xiData{pt,2}(pm*2,:), xiData{pt,2}(modeout*2-1,:), 'linewidth', 2, 'displayname', ['Projection of traj. ',num2str(pt), ' onto $V$'], 'color', nicered)
view(3);

%% Plot evolution
paperFigure('x', 'time [s]', 'y', 'Projected amplitude');
for ii = modesAnalysis
    plot(xiData{pt,1}, vecnorm(xiData{pt,2}(2*ii-1:2*ii,:)), 'linewidth', 2.6, 'displayname', ['Mode ',num2str(ii)])
end

%% fastSSM
[Mmap, iMmap, Tmap, iTmap, Nflow] = fastSSMplus(yDataTrunc, SSMDim, SSMOrder, ROMOrder, NFOrder, 'poly', T(:,modes*2+[-1;0]));

%% Transformation to normal form
clear zDataTrunc
for iTraj = 1:nTraj
    fprintf('transforming traj. %d out of %d to normal form...\n\n', iTraj, nTraj)
    zDataTrunc{iTraj,1} = yDataTrunc{iTraj,1};
    for iN = 1:size(yDataTrunc{iTraj,1},2)
        if mod(iN-1, 250) == 0
            fprintf('computed %5d points out of %5d... %5.1f %% done\n', iN-1, size(yDataTrunc{iTraj,1},2), 100*((iN-1)/size(yDataTrunc{iTraj,1},2)/nTraj+(iTraj-1)/nTraj))
        end
        zDataTrunc{iTraj,2}(:,iN) = iTmap(iMmap(yDataTrunc{iTraj,2}(:,iN)));
    end
    fprintf('computed %5d points out of %5d... %5.1f %% done\n\n', iN, size(yDataTrunc{iTraj,1},2), 100*iTraj/nTraj)
end

%% Reconstruction
zRec = integrateFlows(@(z)Nflow(0,z), zDataTrunc);
yRec = transformTrajectories(@(z)Mmap(Tmap(z)), zRec);
NMTE = computeTrajectoryErrors(yRec, yDataTrunc)

%% Plot decay prediction
pt = 1;
paperFigure('x', 'time', 'y', '$h_{x=0}$ [mm]');
plot(yDataTrunc{pt,1}, width/100*yDataTrunc{pt,2}(outdof,:), 'Color', [0,0,0], 'LineWidth', 2.6, 'DisplayName', 'Original');
plot(yRec{pt,1}, width/100*yRec{pt,2}(outdof,:), '--', 'Color', [0.2,0.9,0.2], 'LineWidth', 2.6, 'DisplayName', 'Reconstruction');
xlim([0,sliceInt(2)])

%% Plot normal form coordinates
paperFigure('x','$\rho_1$','y','$\rho_2$','z','$\rho_4$','legendcols',3);
for iTraj = 1:nTraj
    rho = abs(zDataTrunc{iTraj,2});
    if length(modes) > 2
    plot3(rho(1,:),rho(3,:),rho(5,:), 'linewidth', 4, 'DisplayName', ['Traj. ', num2str(iTraj)])
    else
    plot(rho(1,:),rho(3,:), 'linewidth', 4, 'DisplayName', ['Traj. ', num2str(iTraj)])
    end
end

%% Visualize normal form
pmodes = [1,2];
modeout = 1;
imretype = 'i';
if imretype == 'r'
    zl = ['$\dot{\rho}_',num2str(modeout),'/\rho_',num2str(modeout),'$'];
    imre = @real;
else
    zl = ['$\dot{\theta}_',num2str(modeout),'$'];
    imre = @imag;
end
paperFigure('x',['$\rho_',num2str(pmodes(1)),'$'],'y',['$\rho_',num2str(pmodes(2)),'$'],'z',zl,'legendcols',0);
for iTraj = 1:nTraj
    z = zDataTrunc{iTraj,2};
    rho = abs(z);
    zDot = Nflow(zDataTrunc{iTraj,1}, z);
    dot = imre(zDot./z);
    plot3(rho(pmodes(1)*2-1,:),rho(pmodes(2)*2-1,:),dot(2*modeout-1,:), 'linewidth', 3, 'DisplayName', ['Trajectory ', num2str(iTraj)])
end
[X,Y] = meshgrid(1e-6:0.02:1,1e-6:0.02:1);
XY = 1e-6*ones(SSMDim,length(X(:)));
XY((pmodes(1)*2-1):(pmodes(1)*2),:) = [X(:),X(:)]';
XY((pmodes(2)*2-1):(pmodes(2)*2),:) = [Y(:),Y(:)]';
Nsurf = imre(Nflow(0,XY)./XY);
Nsurf = reshape(Nsurf(2*modeout-1,:)',size(X));
surf(X, Y, Nsurf, 'facecolor', [0,0,0], 'EdgeColor', [0,0,0], 'facealpha', 0.1, 'edgealpha', 0.3)
view(3)

%% Display surface profile
iFreq = 2;
xCoords = linspace(0,width,nPoints);
paperFigure('x','$x$ [mm]','y','elevation [mm]');
frames = [];
xlim([0,width])
ylim([-130,250])
counter = 0;
for iTime = 1:2:length(yRec{1,1})/10
    cla
    area(xCoords, ...
        yDataTrunc{iFreq,2}(1:nPoints,iTime)*width/100, ...
        -130, 'FaceColor', [0.5,0.8,1], 'DisplayName', 'Experiment');
    plot(xCoords, yRec{iFreq,2}(1:nPoints,iTime)*width/100, ...
        ':', 'color', nicered, 'LineWidth', 8, 'DisplayName', 'Simulation');
    counter = counter+1;
    xlabel('$x$ [mm]', 'interpreter', 'latex'); xticklabels([0,100,200,300,400,500]);
    ylabel('Elevation $h$ [mm]', 'interpreter', 'latex'); yticklabels('auto');
    legend('location', 'NW');
    text(300, 220, ['$t = ', num2str(round(yRec{1,1}(iTime)*1000)/1000), '$ s'], 'fontsize', 28, 'interpreter', 'latex')
    frames = [frames, getframe];
end

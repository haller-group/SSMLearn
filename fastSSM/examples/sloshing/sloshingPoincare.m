%% Finding a 2D SSM from transient forced sloshing data
% See J. Axas, B. Bauerlein, K. Avila, and G. Haller. Data-driven modeling 
% of subharmonic forced response due to nonlinear resonance. Scientific 
% Reports, 2024.

clearvars
close all
clc
nicered = [0.7,0.1,0.1];
niceblue = [0.1,0.1,0.7];
nicegreen = [0.1,0.9,0.1];
nicegray = [0.6,0.6,0.6];

load sloshingForcedData.mat
% Note: Data dimensionality has been scaled down 20 times to save repo 
% space. To retrieve the exact results from Axas et al (2024), download the
% full dataset from the ETH library (url pending)

nTraj = size(xData, 1);
nSurfacePoints = size(xData{1,2}, 1);

embedDim = 22;
shiftsteps = 7;
yData = embedCoordinates(xData, embedDim, shiftsteps);
for iTraj = 1:nTraj
    yData{iTraj,3} = xData{iTraj,3};
end

%% Preprocessing
sliceInt = [40,Inf; 55,Inf; 5,Inf; 25,Inf];
yDataTrunc = sliceTrajectories(yData, sliceInt);
pt = 4;
cap = 5400;
paperFigure('x',['$y_{',num2str(1),'}$'],'y',['$y_{',num2str(nSurfacePoints*ceil(embedDim/2)+1),'}$'],'z',['$y_{',num2str(nSurfacePoints*(embedDim-1)+1),'}$'],'legendcols',1);
plot3(yDataTrunc{pt,2}(1,cap:end), yDataTrunc{pt,2}(nSurfacePoints*ceil(embedDim/2)+1,cap:end), yDataTrunc{pt,2}(nSurfacePoints*(embedDim-1)+1,cap:end), 'color', nicegreen,'DisplayName','Embedded 3-periodic orbit')
view(53, 27);
paperFigure('x','time [s]','y','$h_{x=0}$','legendcols',0);
plot(yData{pt,1}, yData{pt,2}(1,:), 'linewidth', 2, 'color', nicegreen)

%% Sample Poincare map
for iTraj = 1:nTraj
    [~,locs] = findpeaks(yDataTrunc{iTraj,3}');
    locs = locs+75; % shift sampling to the approximate max of the response
    locs = locs(1:end-1);
    yMap{iTraj,1} = yDataTrunc{iTraj,1}(:,locs);
    yMap{iTraj,2} = yDataTrunc{iTraj,2}(:,locs);
end

paperFigure('x','time','y','Tank displacement [mm]');
plot(xData{iTraj,1}, xData{iTraj,3}, 'r', 'linewidth', 2)
plot(yDataTrunc{iTraj,1}(locs), yDataTrunc{iTraj,3}(locs), 'ko', 'linewidth', 2)
legend({'Tank motion', 'Map sampling point'});
paperFigure('x','time','y','$h_0$ [mm]');
plot(yDataTrunc{iTraj,1}, yDataTrunc{iTraj,2}(2,:), 'b', 'linewidth', 2)
plot(yDataTrunc{iTraj,1}(locs), yDataTrunc{iTraj,2}(2,locs), 'ko', 'linewidth', 2)
legend({'Surface height', 'Map sampling point'});

%% Center data
ySteady = sliceTrajectories(yMap, [80,Inf]);
X = [];
for iTraj = 1:nTraj
    X = [X, ySteady{iTraj,2}];
end
yMean = mean(X, 2);

%% Identify SSM in Poincare map
SSMOrder = 2;
ROMOrder = 3;
SSMDim = 2;
indTrain = [4];
yMapCentered = yMap;
for iTraj = 1:nTraj
    yMapCentered{iTraj,2} = yMap{iTraj,2} - yMean;
end
[iMmap, Rflow, ~, MmapCentered] = fastSSMMap(yMapCentered(indTrain,:), SSMDim, SSMOrder, ROMOrder);
Mmap = @(eta) MmapCentered(eta) + yMean;
for iTraj = 1:nTraj
    etaMap(iTraj,:) = {yMapCentered{iTraj,1},iMmap(yMapCentered{iTraj,2})};
end

etaRec = iterateMaps(@(z)Rflow(0,z), etaMap);

yGeometry = transformTrajectories(Mmap, etaMap);
yRec = transformTrajectories(Mmap, etaRec);
NMTEGeometry = computeTrajectoryErrors(yGeometry, yMap, 1:nSurfacePoints)
NMTE = computeTrajectoryErrors(yRec, yMap, 1:nSurfacePoints)
%% Reduced coordinates
pt = 2;
paperFigure;
plot(etaMap{pt,2}(1,:), etaMap{pt,2}(2,:), 'o', 'linewidth', 0.5, 'markersize', 10, 'markeredgecolor', 'k', 'markerfacecolor', nicered, 'Color', nicered)
xlabel(['$\xi_',num2str(1),'$'], 'interpreter', 'latex')
ylabel(['$\xi_',num2str(2),'$'], 'interpreter', 'latex')
if NMTE(pt) < 2
    plot(etaRec{pt,2}(1,:), etaRec{pt,2}(2,:), '-s', 'linewidth', 1, 'markersize', 20, 'markeredgecolor', 'k', 'color', nicegreen)
end
legend({'Sampled data','Prediction'})
%% Plot SSM in observable space
etaDataTrunc = transformTrajectories(@(y)iMmap(y-yMean), yDataTrunc(:,1:2));
paperFigure('x',['$\xi_', num2str(1), '$'], 'y',['$\xi_', num2str(2), '$'], 'z','$h_{x=0}$', 'legendcols', 1);
plot3(etaDataTrunc{pt,2}(1,:), etaDataTrunc{pt,2}(2,:), yDataTrunc{pt,2}(1,:), 'color', 'k', 'linewidth', 0.3)
plot2DSSM(1, etaMap{pt,2}, Mmap, 50, 50, nicered);
plot3(etaMap{pt,2}(1,:), etaMap{pt,2}(2,:), yMap{pt,2}(1,:), 'o', 'markersize', 8, 'markeredgecolor', 'k', 'markerfacecolor', nicered, 'linewidth', 0.8, 'color', nicered)
legend({'Continuous trajectory', 'SSM', 'Sampling'})
view(45, 20)

%% surface profile prediction
times = 61:63;
xPlot = linspace(0, 500, nSurfacePoints);
for iPlot = 1:length(times)
    if iPlot == 1
        paperFigure('x', '$x$ [mm]', 'y', 'Elevation $h$ [mm]', 'legendcols', 1);
        yticklabels({-200:100:200})
    elseif iPlot == length(times)
        paperFigure('x', '$x$ [mm]', 'y', '', 'legendcols', 0);
        yticklabels({})
    else
        paperFigure('x', '$x$ [mm]', 'y', '', 'legendcols', 0);
        yticklabels({})
    end
    yticks(-200:100:200);
    grid off
    axis equal
    xlim([0,500])
    ylim([-200,300])
    area(xPlot, yMap{pt,2}(1:nSurfacePoints,times(iPlot)), -200, 'FaceColor', [0.5,0.8,1], 'DisplayName', 'Experiment')
    plot(xPlot, yRec{pt,2}(1:nSurfacePoints,times(iPlot)), ':', 'color', nicered, 'LineWidth', 8, 'DisplayName', 'SSM model')
    if iPlot == 1; legend('location', 'north'); end
end

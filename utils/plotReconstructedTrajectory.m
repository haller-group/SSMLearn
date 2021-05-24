function plotReconstructedTrajectory(yData, yRec, plotCoord, varargin)
% Plot reconstructed coordinate with index plotCoord over time in tData. 
% Plots the given coordinate from a full trajectory yData and as reconstructed
% with SSMFunction and V.
%
% INPUT
% yData        cell (nTraj x 2)   full trajectories
% yRec         cell (nTraj x 2)   reconstructed trajectories
% plotCoord    int                index in dim to plot

fullColor = 'r';
plotTitles = 0;
if ~isempty(varargin)
    fullColor = varargin{1};
end
if length(varargin) > 1
    titles = varargin{2};
    plotTitles = 1;
end
    
nTraj = size(yData,1);
fig = figure;
nCol = ceil(sqrt(nTraj));
nRow = nCol - (nCol * nCol - nTraj > nCol - 1);
tiledlayout(nRow,nCol, 'TileSpacing', 'compact')
maxAmp = max(max(horzcat(yData{:,2})));
minTime = min(min(horzcat(yData{:,1})));
maxTime = max(max(horzcat(yData{:,1})));

for iTraj = 1:nTraj
    nexttile
    hold on
    plot(yData{iTraj,1}, yData{iTraj,2}(plotCoord,:), fullColor, 'Linewidth', 0.3, 'DisplayName', 'Full')
    hh = plot(yRec{iTraj,1}, yRec{iTraj,2}(plotCoord,:), 'k:', 'Linewidth', 2, 'DisplayName', 'Reconstructed');
    hh.Color(4) = .4;
    xlim([-minTime,maxTime]);
    ylim([-maxAmp,maxAmp]);
    if nTraj > 9
        set(gca,'YTick',[])
        set(gca,'XTick', [])
    elseif iTraj == 1
        legend
    end
    if plotTitles
        title(titles{iTraj}, 'Interpreter', 'latex')
    end
    set(gca, 'fontname', 'times')
    set(gca, 'fontsize', 10)
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han, ['$y_{', num2str(plotCoord), '}$'], 'Interpreter', 'latex');
xlabel(han, 'time', 'Interpreter', 'latex');
title(han, 'Trajectory projection onto manifold');
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
hold off
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
if ~isempty(varargin)
    fullColor = varargin{1};
end
    
figure
hold on

nTraj = size(yData,1);
for iTraj = 1:nTraj
    plot(yData{iTraj,1}, yData{iTraj,2}(plotCoord,:), fullColor, 'Linewidth', 1, 'DisplayName', 'Full')
    hh = plot(yRec{iTraj,1}, yRec{iTraj,2}(plotCoord,:), 'k:', 'Linewidth', 2, 'DisplayName', 'Reconstructed');
    hh.Color(4) = .4;
end
xlabel('time', 'Interpreter', 'latex')
ylabel(['$y_{', num2str(plotCoord), '}$'], 'Interpreter', 'latex')
title('Trajectory projection onto manifold')
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
legend
hold off
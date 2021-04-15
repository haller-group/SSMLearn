function plotReconstructedTrajectory(xData, xRec, plotCoord, varargin)
% Plot reconstructed coordinate with index plotCoord over time in tData. 
% Plots the given coordinate from a full trajectory xData and as reconstructed
% with SSMFunction and V.
%
% INPUT
% xData        cell (nTraj x 2)   full trajectories
% xRec         cell (nTraj x 2)   reconstructed trajectories
% plotCoord    int                index in dim to plot

fullColor = 'r';
if ~isempty(varargin)
    fullColor = varargin{1};
end
    
figure
hold on

nTraj = size(xData,1);
for iTraj = 1:nTraj
    plot(xData{iTraj,1}, xData{iTraj,2}(plotCoord,:), fullColor, 'Linewidth', 1, 'DisplayName', 'Full')
    hh = plot(xRec{iTraj,1}, xRec{iTraj,2}(plotCoord,:), 'k:', 'Linewidth', 2, 'DisplayName', 'Reconstructed');
    hh.Color(4) = .6;
end
xlabel('time', 'Interpreter', 'latex')
ylabel(['$q_{', num2str(plotCoord), '}$'], 'Interpreter', 'latex')
title('Trajectory projection onto manifold')
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
legend
hold off
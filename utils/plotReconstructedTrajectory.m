function plotReconstructedTrajectory(xData, xRec, plotCoord)
% Plot reconstructed coordinate with index plotCoord over time in tData. 
% Plots the given coordinate from a full trajectory xData and as reconstructed
% with SSMFunction and V.
%
% INPUT
% xData        cell (nTraj x 2)   full trajectories
% xRec         cell (nTraj x 2)   reconstructed trajectories
% plotCoord     int         index in dim to plot

figure
hold on

nTraj = size(xData,1);
for iTraj = 1:nTraj
    plot(xData{iTraj,1}, xData{iTraj,2}(plotCoord,:),'k','Linewidth',2,'DisplayName','Full')
    plot(xRec{iTraj,1}, xRec{iTraj,2}(plotCoord,:),'r:','Linewidth',2,'DisplayName','Reconstructed')
end
xlabel('time', 'Interpreter', 'latex')
ylabel(['$q_{', num2str(plotCoord), '}$'], 'Interpreter', 'latex')
title('Trajectory projection onto manifold')
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
legend
hold off
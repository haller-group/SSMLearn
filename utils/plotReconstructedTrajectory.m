function plotReconstructedTrajectory(tData, xData, SSMFunction, V, plotCoord)
% Plot reconstructed coordinate with index plotCoord over time in tData. 
% Plots the given coordinate from a full trajectory xData and as reconstructed
% with SSMFunction and V.
%
% INPUT
% tData         (1 x n)     time
% xData         (dim x n)   trajectory
% SSMFunction   anonymous function for the SSM parametrisation
% V             (dim x 2)   2D subspace tangent to the SSM at the fixed point
% plotCoord     int         index in dim to plot

figure
hold on

plot(tData, xData(plotCoord,:),'k','Linewidth',2,'DisplayName','Full Trajectory')
xReconstructed = SSMFunction(V'*xData);
plot(tData, xReconstructed(plotCoord,:),'r:','Linewidth',2,'DisplayName','Reconstructed Trajectory')
xlabel('time', 'Interpreter', 'latex')
ylabel(['$q_', num2str(plotCoord), '$'], 'Interpreter', 'latex')
title('Trajectory projection onto manifold')
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
legend
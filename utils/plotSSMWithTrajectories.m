function plotSSMWithTrajectories(xData, SSMFunction, plotInds, plotRadius)
% Draws the shape of an SSM as well as trajectories in cell array xData, 
% in the space of components plotInds.

figure
hold on

plotSSM(plotRadius, plotInds, SSMFunction);

for iTraj = 1:length(xData)
    plot3(xData{iTraj}(plotInds(1),:), xData{iTraj}(plotInds(2),:), xData{iTraj}(plotInds(3),:))
end

xlabel(['$q_' num2str(plotInds(1)) '$'], 'Interpreter', 'latex')
ylabel(['$q_' num2str(plotInds(2)) '$'], 'Interpreter', 'latex')
zlabel(['$q_' num2str(plotInds(3)) '$'], 'Interpreter', 'latex')
title('Computed manifold with test set trajectories')
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
grid on
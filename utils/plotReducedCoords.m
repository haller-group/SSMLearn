function plotReducedCoords(xData, V)
% Plot the (first) two components of V'*xData for a cell array of trajectories
% xData against each other. 

figure
hold on
for i = 1:length(xData)
    y = V'*xData{i};
    plot(y(1,:), y(2,:))
end

xlabel('$\eta_1$', 'Interpreter', 'latex')
ylabel('$\eta_2$', 'Interpreter', 'latex')
title({'Test set trajectories in', 'reduced model coordinates'})
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
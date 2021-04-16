function plotReducedCoords(yData, varargin)
% Plot the (first) two components of yData for a cell array of trajectories
% yData against each other. 

if ~isempty(varargin)
    yRec = varargin{1};
end

figure
hold on
nTraj = size(yData,1);
for iTraj = 1:nTraj
    plot(yData{iTraj,2}(1,:), yData{iTraj,2}(2,:), 'DisplayName', ['Measurement ', num2str(iTraj)])
    if ~isempty(varargin)
        plot(yRec{iTraj,2}(1,:), yRec{iTraj,2}(2,:), ':', 'LineWidth', 2, 'DisplayName', ['Prediction ', num2str(iTraj)])
    end
end

legend();
xlabel('$y_1$', 'Interpreter', 'latex')
ylabel('$y_2$', 'Interpreter', 'latex')
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
hold off
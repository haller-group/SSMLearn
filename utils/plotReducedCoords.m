function plotReducedCoords(etaData, varargin)
% Plot the (first) two components of etaData for a cell array of trajectories
% etaData against each other. 

if ~isempty(varargin)
    etaRec = varargin{1};
end

figure
hold on
nTraj = size(etaData,1);
for iTraj = 1:nTraj
    plot(etaData{iTraj,2}(1,:), etaData{iTraj,2}(2,:), 'DisplayName', ['Measurement ', num2str(iTraj)])
    if ~isempty(varargin)
        plot(etaRec{iTraj,2}(1,:), etaRec{iTraj,2}(2,:), ':', 'LineWidth', 2, 'DisplayName', ['Prediction ', num2str(iTraj)])
    end
end

legend();
xlabel('$\eta_1$', 'Interpreter', 'latex')
ylabel('$\eta_2$', 'Interpreter', 'latex')
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
hold off
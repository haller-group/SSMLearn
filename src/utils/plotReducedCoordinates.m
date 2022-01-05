function plotReducedCoordinates(etaData, varargin)
% plotReducedCoordinates(etaData)
% Plot the (first) two or three components of etaData against each other.
% plotReducedCoordinates(etaData, etaRec)
% Plot the (first) two or three components of etaData and of etaRec. Useful
% to compare a reconstruction from reduced dynamics to the original data.
% plotReducedCoordinates(etaData, etaRec, [3,4,1])
% Plot the components eta_3, eta_4 and eta_1 of etaData and of etaRec.
%
% INPUT
% etaData   {nTraj x 2}   Cell array of trajectories to plot. First column
%                         contains time, second column contains state.
% etaRec    {nTraj x 2}   Optional (reconstructed) trajectories to plot in
%                         dashed style on top of etaData.
% plotinds  [3 x 1]       Optional component indices to plot.

p = inputParser;
validCell = @(x) iscell(x);
addOptional(p, 'etaRec', {}, validCell);
addOptional(p, 'plotinds', [1,2,3]);
parse(p, varargin{:});
etaRec = p.Results.etaRec; plotinds = p.Results.plotinds;

customFigure();
nTraj = size(etaData,1);
if size(etaData{1,2},1)==2
    for iTraj = 1:nTraj
        plot(etaData{iTraj,2}(1,:), etaData{iTraj,2}(2,:), ...
            'DisplayName', ['Measurement ', num2str(iTraj)])
        if ~isempty(etaRec)
            plot(etaRec{iTraj,2}(1,:), etaRec{iTraj,2}(2,:), ':', ...
                'LineWidth', 2, 'DisplayName', ['Prediction ', ...
                                                           num2str(iTraj)])
        end
    end
    legend();
xlabel('$\eta_1$', 'Interpreter', 'latex')
ylabel('$\eta_2$', 'Interpreter', 'latex')
else
    for iTraj = 1:nTraj
        plot3(etaData{iTraj,2}(plotinds(1),:), etaData{iTraj,2}(plotinds(2),:),  ...
            etaData{iTraj,2}(plotinds(3),:),'DisplayName', ['Measurement ', ...
            num2str(iTraj)])
        if ~isempty(etaRec)
            plot3(etaRec{iTraj,2}(plotinds(1),:), etaRec{iTraj,2}(plotinds(2),:), ...
                etaRec{iTraj,2}(plotinds(3),:), ':', 'LineWidth', 2, ...
                'DisplayName', ['Prediction ', num2str(iTraj)])
        end
    end

legend();
xlabel(['$\eta_{',num2str(plotinds(1)),'}$'], 'Interpreter', 'latex')
ylabel(['$\eta_{',num2str(plotinds(2)),'}$'], 'Interpreter', 'latex')
zlabel(['$\eta_{',num2str(plotinds(3)),'}$'], 'Interpreter', 'latex')
end
function plotSSMWithTrajectories(IMInfo, plotInds, yData, varargin)
% plotSSMWithTrajectories(IMInfo, plotInds, yData, 'Margin', m)
%
% Draws the shape of a 1D or 2D SSM along with the trajectories in yData, 
% in the space of components plotInds.
%   
% INPUT
% IMInfo     struct         Manifold from IMGeometry function. 
% plotInds   int or (3x1)   Either a vector of three components of the
%                           observable space to plot, or a single component
%                           to plot on the z axis. If a scalar value is
%                           passed the x and y axes are aligned with the
%                           reduced coordinate directions.
% yData      {nTraj x 2}    Cell array of trajectories. First column 
%                           contains time, second column contains state.
% Margin     real, percent  Optional. Default 20. Controls the margin
%                           between the manifold edge and the maximum
%                           extent of trajectories
% Colors     (nAmp x 3)     Optional RGB colors, one row per curve, or pass
%                           a single row vector for all curves.
% ColorSurf     (1 x 3)     Optional RGB color for the SSM surface, ...
%                           default [-1 -1 -1]='cool'

% EXAMPLES
% Plot a 1D or 2D manifold in the first three components of the observable
% space:
%    plotSSMWithTrajectories(IMInfo, [1,2,3], yData)
% Plot a 2D manifold with y_1 on the z axis and the reduced coordinates
% eta_1 and eta_2 on the x and y axes, with a margin of 50 %:
%    plotSSMWithTrajectories(IMInfo, 1, yData, 'Margin', 50)

customFigure(); colors = colororder;
p = inputParser;
addOptional(p, 'Margin', 20, @(x) isnumeric(x) && isscalar(x));
addOptional(p, 'Colors', colors);
addOptional(p, 'ColorSurf', [-1 -1 -1]);
parse(p, varargin{:});
colororder(p.Results.Colors)


SSMDim = IMInfo.parametrization.dimension;
SSMChart = IMInfo.chart.map;
SSMFunction = IMInfo.parametrization.map;

if SSMDim == 2
    etaData = SSMChart(cat(2,yData{:,2}));
    plot2DSSM(plotInds, etaData, SSMFunction, p.Results.Margin, 50, p.Results.ColorSurf);
elseif SSMDim == 1
    etaData = SSMChart(cat(2,yData{:,2}));
    etaMin = min(etaData) * (100+p.Results.Margin)/100;
    etaMax = max(etaData) * (100+p.Results.Margin)/100;
    parametrization = linspace(etaMin, etaMax);
    dVals = SSMFunction(parametrization);
    mfdplot = plot3(dVals(plotInds(1),:), dVals(plotInds(2),:), dVals(plotInds(3),:), 'b', 'LineWidth', 7);
    mfdplot.Color(4) = 0.1;
else
    disp("SSM plotting only available for 1D and 2D manifolds")
end

for iTraj = 1:size(yData,1)
    if length(plotInds) == 3
        plot3(yData{iTraj,2}(plotInds(1),:), yData{iTraj,2}(plotInds(2),:), yData{iTraj,2}(plotInds(3),:), 'LineWidth', 2)
        xlabel(['$y_{' num2str(plotInds(1)) '}$'], 'Interpreter', 'latex')
        ylabel(['$y_{' num2str(plotInds(2)) '}$'], 'Interpreter', 'latex')
        zlabel(['$y_{' num2str(plotInds(3)) '}$'], 'Interpreter', 'latex')
    elseif length(plotInds) == 1
        etaData = transformTrajectories(SSMChart, yData);
        plot3(etaData{iTraj,2}(1,:), etaData{iTraj,2}(2,:), yData{iTraj,2}(plotInds,:), 'LineWidth', 2)
        xlabel(['$\eta_1$'], 'Interpreter', 'latex')
        ylabel(['$\eta_2$'], 'Interpreter', 'latex')
        zlabel(['$y_{' num2str(plotInds) '}$'], 'Interpreter', 'latex')
    end
end
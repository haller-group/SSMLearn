function plotSSMWithTrajectories(yData, SSMFunction, plotInds, V, radiusIncr, varargin)
% Draws the shape of an SSM as well as trajectories in cell array yData, 
% in the space of components plotInds.

opts = struct('SSMDimension',2);
if rem(length(varargin),2) > 0 && length(varargin) > 1
    error('Error on input arguments. Missing or extra arguments.')
end
% Custom options
if nargin > 3
    for ii = 1:length(varargin)/2
        opts = setfield(opts,varargin{2*ii-1},...
            varargin{2*ii});
    end
end

figure
hold on

if opts.SSMDimension == 2
    etaData = transpose(V) * cat(2,yData{:,2});
    plot_2dSSM_surf(plotInds, etaData, SSMFunction, radiusIncr, 50, 0);
elseif opts.SSMDimension == 1
    etaData = transpose(V) * cat(2,yData{:,2});
    etaMin = - max(abs(etaData)) * (100+radiusIncr)/100;
    etaMax = + max(abs(etaData)) * (100+radiusIncr)/100;
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
        etaData = getProjectedTrajs(yData, V);
        plot3(etaData{iTraj,2}(1,:), etaData{iTraj,2}(2,:), yData{iTraj,2}(plotInds,:), 'LineWidth', 2)
        xlabel(['$\eta_1$'], 'Interpreter', 'latex')
        ylabel(['$\eta_2$'], 'Interpreter', 'latex')
        zlabel(['$y_{' num2str(plotInds) '}$'], 'Interpreter', 'latex')
    end
end

title('Computed manifold with trajectories')
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
grid on
hold off
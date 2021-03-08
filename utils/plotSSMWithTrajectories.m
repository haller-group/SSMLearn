function plotSSMWithTrajectories(xData, SSMFunction, plotInds, V, radiusIncr, varargin)
% Draws the shape of an SSM as well as trajectories in cell array xData, 
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
    yData = V'* cat(2,xData{:,2});
    plot_2dSSM_surf(plotInds, yData, SSMFunction, radiusIncr, 50, 0);
else
    disp("SSM plotting only available for 2D manifolds")
end

for iTraj = 1:length(xData)
    plot3(xData{iTraj,2}(plotInds(1),:), xData{iTraj,2}(plotInds(2),:), xData{iTraj,2}(plotInds(3),:))
end

xlabel(['$q_{' num2str(plotInds(1)) '}$'], 'Interpreter', 'latex')
ylabel(['$q_{' num2str(plotInds(2)) '}$'], 'Interpreter', 'latex')
zlabel(['$q_{' num2str(plotInds(3)) '}$'], 'Interpreter', 'latex')
title('Computed manifold with test set trajectories')
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
grid on
hold off
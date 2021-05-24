function plotTrajs(xData, varargin)

p = inputParser;
addOptional(p, 'colors', {'b'});
addOptional(p, 'legendnames', {''});
addOptional(p, 'styles', {''});
parse(p, varargin{:});

minTime = min(min(horzcat(xData{:,1})));
maxTime = max(max(horzcat(xData{:,1})));

nTraj = size(xData, 1);
figure; hold on;
for iTraj = 1:nTraj
    plot(xData{iTraj,1}, xData{iTraj,2}, p.Results.styles{iTraj}, 'DisplayName', ...
        p.Results.legendnames{iTraj}, 'LineWidth', 0.8, 'Color', p.Results.colors{iTraj});
end
legend('Interpreter', 'latex')
xlim([minTime, maxTime])
xlabel('Time [s]','Interpreter','latex')
ylabel('$\hat{X}$ [\%]','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
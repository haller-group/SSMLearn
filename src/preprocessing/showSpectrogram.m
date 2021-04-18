function [powerdensity, frequencies, times] = showSpectrogram(xData, varargin)
% Plot a spectrogram for a scalar signal xData{1,2} at times xData{1,1}
% Examples of use:
% showSpectrogram(xData)        plots the first component in xData{1,2}
% showSpectrogram(xData, 10)    plots the 10th component in xData{1,2}

plotcoord = 1;
if ~isempty(varargin)
    plotcoord = varargin{1};
end

stfourier = 0;
for iTraj = 1:size(xData, 1)
    t = xData{iTraj,1}; x = xData{iTraj,2}(plotcoord, :);
    Nwin = round(length(t)/50);
    % [stfourier, frequencies, times] = spectrogram(x, Nwin, round(Nwin*0.5), [], 1./(t(2)-t(1)));
    [stf, frequencies, times] = spectrogram(x, [], round(Nwin*0.5), [], 1./(t(2)-t(1))); % TODO: If timeseries vary between trajectories
    stfourier = stfourier + stf;
end

% spectrogram(x, Nwin, round(Nwin*0.5), [], 1./(t(2)-t(1)), 'yaxis');
% spectrogram(x, 'yaxis');
powerdensity = abs(stfourier);

surf(times, frequencies, powerdensity)
set(gca,'ColorScale','log')
xlim([min(times), max(times)])
ylim([min(frequencies), max(frequencies)])
xlabel('$t \, [$s$]$', 'Interpreter', 'latex')
ylabel('$f \, [$Hz$]$', 'Interpreter', 'latex')
set(gca,'fontname', 'times')
set(gca,'fontsize', 18)
view(2)
shading interp
c = colorbar;
c.Label.String = 'power spectral density [1/Hz]';
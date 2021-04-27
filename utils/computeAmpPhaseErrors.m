function [redAmplitudeError, redPhaseError] = computeAmpPhaseErrors(etaRec, etaData, varargin)

plotErrors = 1;
if length(varargin) == 1
    plotErrors = varargin{1};
end

nTraj = size(etaData,1);

redAmplitudeError = zeros(nTraj,1);
redPhaseError = zeros(nTraj,1);

for iTraj = 1:nTraj
    rhoData{iTraj,:} = vecnorm(etaData{iTraj,2});
    rhoRec{iTraj,:} = vecnorm(etaRec{iTraj,2});
    
    thetaData{iTraj,:} = unwrap(atan2(etaData{iTraj,2}(1,:), etaData{iTraj,2}(2,:)));
    thetaRec{iTraj,:} = unwrap(atan2(etaRec{iTraj,2}(1,:), etaRec{iTraj,2}(2,:)));
    
    redAmplitudeError(iTraj) = mean(abs(rhoRec{iTraj,:} - rhoData{iTraj,:})) / mean(rhoData{iTraj,:});
    redPhaseError(iTraj) = mean(abs(thetaRec{iTraj,:} - thetaData{iTraj,:})) / mean(abs(thetaData{iTraj,:}));
end

if plotErrors
    for iTraj = 1:nTraj
        figure(97)
        plot(etaData{iTraj,1}, rhoData{iTraj,:}, 'DisplayName', ['Measurement ' num2str(iTraj)], 'LineWidth', 2)
        hold on
        plot(etaRec{iTraj,1}, rhoRec{iTraj,:}, ':', 'DisplayName', ['Prediction ' num2str(iTraj)], 'LineWidth', 2);
        legend
        title('Amplitude')
        xlabel('time [s]')
        ylabel('amplitude')
        figure(98)
        plot(etaData{iTraj,1}, thetaData{iTraj,:}, 'DisplayName', ['Measurement ' num2str(iTraj)], 'LineWidth', 2)
        hold on
        plot(etaRec{iTraj,1}, thetaRec{iTraj,:}, ':', 'DisplayName', ['Prediction ' num2str(iTraj)], 'LineWidth', 2);
        legend
        title('Phase')
        xlabel('time [s]')
        ylabel('$\theta$ [rad]', 'Interpreter', 'latex')
    end
end
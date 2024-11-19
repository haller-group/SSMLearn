function [startTime, indStartTime, SSMDim] = SSM_startTime(data,indplot)
    t = data{indplot,1};

    [powerdensity, frequencies, times] = showSpectrogram(data, indplot);
    epsilon = 0.1;
    dominantFreqs_tot = analyzeSpectr(times,frequencies,powerdensity,epsilon);
    
    n_vect = zeros(1,length(times));
    Dimension_SSM = zeros(1,length(times));
    for ii = 1:length(times)
        n = length(dominantFreqs_tot{1,ii});
        n_vect(ii) = n;

        if n == 1 || n == 0
            Dimension_SSM(ii) = 2;
        else
            Dimension_SSM(ii) = 2*n;
        end

    end
    
    fig = customFigure();
    xlabel('time [s]')
    ylabel('SSM dimension')
    plot(times, Dimension_SSM,'k','LineWidth',2)

    xlim([times(1) times(end)]);
    ylim([0 max(Dimension_SSM) + 1]);

    time_single_harm = times(n_vect == 1);
    startTime = time_single_harm(1);
    [~,indStartTime] = min(abs(t-startTime));
    SSMDim = Dimension_SSM(end);
    
end
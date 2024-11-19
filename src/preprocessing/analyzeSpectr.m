function dominantFreqs_tot = analyzeSpectr(times,frequencies,powerdensity,epsilon)

dominantFreqs_tot = cell(1,length(times));

    for ii = 1:length(times)
    
        slice = powerdensity(:,ii);
%         [pks, locs] = findpeaks(slice);
        [pks, locs] = findpeaks(slice,'MinPeakDistance',10);
        if ii == 1
            max_pks = max(pks);
            minimum_pks = epsilon * max_pks;
        end

        dominantPeaks = pks(pks>=minimum_pks);
        dominantLocs = locs(pks>=minimum_pks);
        dominantFreqs = frequencies(dominantLocs);
        dominantFreqs_tot{ii} = dominantFreqs;

%         figure; hold on; grid on; 
%         xlabel('frequency [rad/s]')
%         ylabel('Power spectral density [1/Hz]')
%         plot(frequencies, slice,'k.-','MarkerSize',10,'linewidth',1)
%         plot(frequencies(locs), slice(locs),'ro','MarkerSize',10,'linewidth',2)
%         plot(frequencies(dominantLocs), slice(dominantLocs),'bo','MarkerSize',10,'linewidth',2)

    end

end
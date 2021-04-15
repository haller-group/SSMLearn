function tstart = pickTruncationFromSpectrogram(powerdensity, times, tolerance)

[~, locs] = findpeaks(powerdensity(:,1));

relativePowers = powerdensity(locs(2:end),:)./powerdensity(locs(1),:);
startIndex = find(relativePowers < tolerance, 1);

tstart = times(startIndex);
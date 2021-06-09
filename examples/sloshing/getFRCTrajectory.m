function [omega, x] = getFRCTrajectory(xData)

nTraj = size(xData,1);
for iTraj = 1:nTraj
    [amp,freq,damp,time] = PFF(xData{iTraj,1}, xData{iTraj,2}(1,:));
    
    len = floor(min(length(freq), length(amp))/2);
    omega{iTraj} = freq(1:len)*2*pi;
    x{iTraj} = amp(1:len);
end
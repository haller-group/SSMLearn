function [amp,ph,freq]=DFT(X,dt)
% DFT of a Time series matrix n.series x n.time_inst
[n,N]  = size(X);
if n>N
    X=transpose(X);
    [n,N]  = size(X);
end

df = 1/(N*dt);

Y  = fft(X,[],2);

if mod(N,2)==0
    freq = 0:df:(N/2*df);
    Y_pos = zeros(n,length(freq));
    Y_pos(:,1) = Y(:,1)/N;
    Y_pos(:,2:N/2) = 2*Y(:,2:N/2)/N;
    Y_pos(:,N/2+1) = Y(:,N/2+1)/N;
else
    freq = 0:df:((N-1)/2*df);
    Y_pos = zeros(n,length(freq));
    Y_pos(:,1) = Y(:,1)/N;
    Y_pos(:,2:(N+1)/2) = 2*Y(:,2:(N+1)/2)/N;
end;

amp=abs(Y_pos);
ph=angle(Y_pos);
end
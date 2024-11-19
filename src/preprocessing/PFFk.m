function [amp,freq,damp,time,t_zeros] = PFFk(t,x,kmean,varargin)
% [amp,freq,damp,time] = PFF(t,x,varargin)
%
% Implementation of the PFF algorithm for extraction of instantaneous
% damping and frequency of the time signal x known at times t. If the
% optional argument is given, this should be a vector of two values which
% are trated as two frequency in Hz that are used to filter the signal
% before processing

% See the paper below for more
% M. Jin, W. Chen, M. R. W. Brake, and H. Song. Identification of 
% instantaneous frequency and damping from transient decay data. Journal of 
% Vibration and Acoustics, 142(5):051111.

% Eventual filtering
if isempty(varargin)==0
   cut_freq = varargin{:};
   x = bandpass(x,2*pi*cut_freq,1/(t(2)-t(1))); 
end

if size(x,1)>size(x,2); x = transpose(x); end 
if size(t,1)>size(t,2); t = transpose(t); end 

% Find zero crossing points to compute frequencies
xx_1 = x(1:end-1).*x(2:end);
loc_change = find(xx_1<0);
t_a = t(loc_change); x_a = x(loc_change);
t_b = t(loc_change+1); x_b = x(loc_change+1);
t_zeros = t_a-x_a.*(t_b-t_a)./(x_b-x_a);
perd = diff(t_zeros)*2;

% Assume the signal to start right before a zero and to terminate right
% after it

idx = loc_change(1):(loc_change(end)+1);
x = x(idx); t = t(idx);


% Find maximal points to compute amplitudes, damping and times
xa = abs(x);
% [~,loc] = findpeaks(xa);
[~,loc] = findpeaks(xa,'MinPeakDistance',mean(diff(loc_change))/2);
X = [xa(loc-1); xa(loc); xa(loc+1)];
T = [t(loc-1);   t(loc);  t(loc+1)];

% figure
% hold on
% % plot(t,x,'r')
% % plot(t,x,'k*')
% plot(t_zeros,zeros(size(t_zeros)),'co');
% % plot(t(loc),x(loc),'go');
% plot(t,xa,'k*')
% plot(T(1,:),X(1,:),'g*')
% plot(T(2,:),X(2,:),'r*')
% plot(T(3,:),X(3,:),'b*')

amp = zeros(1,length(loc)); time = amp;

% for ii = 1:length(loc)
%    A = [T(:,ii).^2 T(:,ii) ones(3,1)];
%    P = A\X(:,ii);
%    time(ii) = -P(2)/2/P(1);
%    amp(ii) = -P(2).^2/4/P(1)+P(3);
% end

for ii = 1:length(loc)
   TT = T(:,ii)-T(1,ii);
   A = [TT.^2 TT ones(3,1)];
   P = A\X(:,ii);
   time(ii) = -P(2)/2/P(1);
   time(ii) = time(ii) + T(1,ii);
   amp(ii) = -P(2).^2/4/P(1)+P(3);
end

damp = diff(amp)./diff(time); damp = [damp damp(end)]./amp;

% Cleaning with moving mean
%kmean = 4;
perd = movmean(perd,kmean);
freq = 1./perd;
damp = movmean(damp,kmean);
freq = interp1(t_zeros(2:end)-(t_zeros(2)-t_zeros(1))/2, freq, time);
end



clearvars
close all

nElements = 2;
E       = 70e9;   % Young's modulus
rho     = 2700;   % density
nu      = 0.3;    % nu
kappa   = 3e6;    % linear damping
l       = 1;      % beam length
h       = 20e-3;  % height
b       = 50e-3;  % width

[M,C,K,fnl,fext,outdof] = von_karman_model(nElements, E, rho, nu, kappa, l, h, b);
n = size(M,1);    % mechanical dofs (axial def, transverse def, angle)
omega = linspace(95,110,10);

w0 = -K\fext; % linear initial guess
IC = [w0; zeros(n,1)];

for ii = 1:length(omega)
[F, lambda] = functionFromTensors(M, C, K, fnl, 3.5*fext, omega(ii));
observable = @(x) x;
tEnd = 100*2*pi/omega(ii);
nSamp = fix(50 * tEnd * abs(imag(lambda(1))) / (2*pi));
dt = tEnd/(nSamp-1);
tic
disp(['run ' num2str(ii), '/' num2str(length(omega)) ', omega ' num2str(omega(ii))])
xData(1,:) = integrateTrajectories(F, observable, tEnd, nSamp, 1, IC);
toc
ampli(ii) = max(abs(xData{1,2}(n-1,end-100:end)))
IC = xData{1,2}(:,end);
end
plot(omega,ampli,'o')
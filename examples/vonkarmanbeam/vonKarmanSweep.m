clearvars
close all

nElements = 3;
E       = 70e9;   % Young's modulus
rho     = 2700;   % density
nu      = 0.3;    % nu
kappa   = 3e6;    % linear damping
l       = 1;      % beam length
h       = 20e-3;  % height
b       = 50e-3;  % width

[M,C,K,fnl,fext,outdof] = von_karman_model(nElements, E, rho, nu, kappa, l, h, b);
n = size(M,1);    % mechanical dofs (axial def, transverse def, angle)
w_span = [120,200];
omega = linspace(w_span(1),w_span(2),20);
load FRC_data
f_full = [128];

w0 = -K\(f_full*fext); % linear initial guess
IC = [w0; zeros(n,1)];
%%
for ii = 1:length(omega)
[F, lambda] = functionFromTensors(M, C, K, fnl, f_full*fext, omega(ii));
observable = @(x) x;
tEnd = 200*2*pi/omega(ii);
nSamp = fix(50 * tEnd * abs(imag(lambda(1))) / (2*pi));
dt = tEnd/(nSamp-1);
tic
disp(['run ' num2str(ii), '/' num2str(length(omega)) ', omega ' num2str(omega(ii))])
xData(1,:) = integrateTrajectories(F, observable, tEnd, nSamp, 1, IC);
toc
ampli(ii) = max(abs(xData{1,2}(n-1,end-100:end)))
IC = xData{1,2}(:,end);
end
%%
FRC_full = getFRC_full(M, C, K, fnl, fext, f_full, n-1, w_span, 7); close all

%%
plot(omega,ampli,'o', 'DisplayName', 'Numerical integration')
hold on
plot(FRC_full.F1.Freq, FRC_full.F1.Amp, 'DisplayName', 'SSMTool')
plot(FRC_data.F2.Freq, FRC_data.F2.Amp, 'DisplayName', 'SSMLearn')
xlabel('$\Omega$ [rad/s]', 'Interpreter', 'latex')
ylabel('$u$ [m]', 'Interpreter', 'latex')
legend
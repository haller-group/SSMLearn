%% Finding a 2D SSM for a von Kármán beam
% 
% 
% This is an example of how to reconstruct a slow 2D SSM of a mechanical system 
% using synthetic measurements of a scalar quantity. In this example, we measure 
% the end tip displacement of a von Kármán beam. [1]
% 
% 
% 
% [1] Jain, S., Tiso, P. & Haller, G. (2017). Exact Nonlinear Model Reduction 
% for a von Karman beam: Slow-Fast Decomposition and Spectral Submanifolds. Journal 
% of Sound and Vibration. 423. 10.1016/j.jsv.2018.01.049. 

clearvars
close all
%% Example setup
% 

nElementsList = [2 3 4 5];
E       = 70e9;   % Young's modulus
rho     = 2700;   % density
nu      = 0.3;    % nu
kappa   = 3e6;    % linear damping
l       = 1;      % beam length
h       = 20e-3;  % height
b       = 50e-3;  % width

for iEl = 1:length(nElementsList)
nElements = nElementsList(iEl)

[M,C,K,fnl,fext,outdof] = von_karman_model(nElements, E, rho, nu, kappa, l, h, b);
n = size(M,1);    % mechanical dofs (axial def, transverse def, angle)
[F, lambda] = functionFromTensors(M, C, K, fnl);


% [Fforced, lambda2] = functionFromTensors(M, C, K, fnl, 7.0*fext, 105);
%% Generation of Synthetic Data
% 

nTraj = 1;
%% 
% 

loadvector = zeros(n,nTraj);
loadvector(n-1,1) = 700; % point load, [N]
w0 = -K\loadvector(:,1); % linear initial guess
options = optimoptions('fsolve', 'MaxFunctionEvaluations', 1400*n, 'MaxIterations', 10000, 'Display', 'off');
IC = zeros(2*n,nTraj);
for iLoad = 1:nTraj
    f_eq = @(w)([zeros(n) M]*F(0,[w; zeros(n,1)]) + loadvector(:,iLoad));
    [w0, ~, exitflag, output] = fsolve(f_eq, w0, options);
    if exitflag <= 0
        error('Warning: No solution found for loading configuration')
    end
    IC(:,iLoad) = [w0; zeros(n,1)];
end
%% 
% 
observable = @(x) x(n-1,:);
tEnd = 10;
nSamp = fix(50 * tEnd * abs(imag(lambda(1))) / (2*pi));
dt = tEnd/(nSamp-1);
tic
% xData(1,:) = integrateTrajectories(Fforced, observable, tEnd, nSamp, 1, IC(:,1));
xData(1,:) = integrateTrajectories(F, observable, tEnd, nSamp, 1, IC, 'odetol', 1e-7);
intTimes(iEl) = toc

%% Delay embedding
% 

indTest = [1];
indTrain = [1];

SSMDim = 2;
overEmbed = 3;
yData = coordinates_embedding(xData, SSMDim, 'OverEmbedding', overEmbed);

sliceInt = [0.1, tEnd];
yDataTrunc = sliceTrajectories(yData, sliceInt);
%% Datadriven manifold fitting

SSMOrder = 3;
[V, SSMFunction, mfdInfo] = IMparametrization(yDataTrunc(indTrain,:), SSMDim, SSMOrder);
%% Plot and validation
etaData = getProjectedTrajs(yData, V);
etaDataTrunc = getProjectedTrajs(yDataTrunc, V);
%% Reduced order model
ROMOrder = 7;
[~,Tinv,N,T,NormalFormInfo] = IMdynamics_flow(etaDataTrunc(indTrain,:), 'R_PolyOrd', ROMOrder, 'style', 'normalform', 'fig_disp_nf', 0);
%% 
zData = transformComplex(Tinv, etaDataTrunc);
%% Backbone curves
% 
N_info = NormalFormInfo.N;
[damp, freq] = polarnormalform(N_info.coeff, N_info.exponents, N_info.phi);
N_info.coeff
%%
hold on
maxRho = abs(zData{indTrain(1),2}(1,1));
[dmp, frq, amp, rho_plot] = backbonecurves(damp, freq, SSMFunction, T, 1, maxRho);
subplot(121); ylabel('$u \, [$m$]$','Interpreter','latex')
subplot(122); ylabel('$u \, [$m$]$','Interpreter','latex')
end
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
xData(iEl,:) = integrateTrajectories(F, observable, tEnd, nSamp, 1, IC, 'odetol', 1e-7);
intTimes(iEl) = toc
end
%%
hBar=bar(nElementsList, intTimes, 'y');
set(gca, 'YScale', 'log')
xlabel('# elements')
ylabel('time [s]')
ylim([1,100])
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
grid on
return
%% Delay embedding
% 

nTraj = length(nElementsList);
indTest = [2];
indTrain = [2];

SSMDim = 2;
overEmbed = 3;
yData = coordinates_embedding(xData, SSMDim, 'OverEmbedding', overEmbed);
%% Data filtering
% 

showSpectrogram(yData(indTrain,:), 1);
ylim([0,150])
%% 
% We plot the tip displacement  over time for closer inspection. 

plot(xData{1,1}, xData{1,2}, xData{2,1}, xData{2,2}, ':');
xlabel('$t \, [$s$]$','Interpreter','latex'); ylabel('$u \, [$m$]$','Interpreter','latex'); 
legend({'trajectory 1', 'trajectory 2'})
set(gca,'fontname', 'times'); set(gca,'fontsize', 18); title('Generated data')
%% 
% 

sliceInt = [0.1, tEnd];
yDataTrunc = sliceTrajectories(yData, sliceInt);
%% Datadriven manifold fitting

SSMOrder = 3;
[V, SSMFunction, mfdInfo] = IMparametrization(yDataTrunc(indTrain,:), SSMDim, SSMOrder);
%% Plot and validation
etaData = getProjectedTrajs(yData, V);
etaDataTrunc = getProjectedTrajs(yDataTrunc, V);
%% 
plotReducedCoords(etaData);
legend({'Test set trajectory', 'Training set trajectory'})
%% 
plotSSMWithTrajectories(yData(indTest,:), SSMFunction, 1, V, 10, 'SSMDimension', SSMDim)
view(-100,20); zlabel('$u \, [$m$]$','Interpreter','latex')
%% Reduced order model
ROMOrder = 7;
[~,Tinv,N,T,NormalFormInfo] = IMdynamics_flow(etaDataTrunc(indTrain,:), 'R_PolyOrd', ROMOrder, 'style', 'normalform');
%% 
zData = transformComplex(Tinv, etaDataTrunc);
[zRec, yRecNormal] = integrateFlows(N, zData, @(q) SSMFunction(T(q)));
etaRecNormal = transformComplex(T, zRec);
%% Evaluation of reduced dynamics
% 
[reducedTrajDist, fullTrajDist] = computeRecDynErrors(etaRecNormal, yRecNormal, etaDataTrunc, yDataTrunc);
RRMSE_normal = mean(fullTrajDist(indTest))
%% 
% 
plotReducedCoords(etaDataTrunc(indTest(1),:), etaRecNormal(indTest(1),:))
legend({'Test set (truncated)', 'Prediction'})
%% 
plotReconstructedTrajectory(yData(indTest(1),:), yRecNormal(indTest(1),:), 1, 'm')
legend({'Test set', 'Prediction'}); ylabel('$u \, [$m$]$','Interpreter','latex')
%% 
% 
DSEigenvalues = lambda(1:SSMDim)
normalFormEigenvalues = computeEigenvaluesFlow(NormalFormInfo)
%% Backbone curves
% 
N_info = NormalFormInfo.N;
[damp, freq] = polarnormalform(N_info.coeff, N_info.exponents, N_info.phi);
figure
maxRho = abs(zData{indTrain(1),2}(1,1));
[dmp, frq, amp, rho_plot] = backbonecurves(damp, freq, SSMFunction, T, 1, maxRho);
subplot(121); ylabel('$u \, [$m$]$','Interpreter','latex')
subplot(122); ylabel('$u \, [$m$]$','Interpreter','latex')
%% 
% The data-driven model can now be used for forced response predictions. We 
% first compute forced response using the full model with SSMTool that are used, 
% after calibration, to validate the predictions of our data-driven reduced-order 
% model.

% Compute with SSMTool
f_full = [0.7000 128.0000];
f_vec = loadvector(:,1)/max(abs(loadvector(:,1)));
w_span = [90 200];
FRC_full = getFRC_full(M, C, K, fnl, f_vec, f_full, n-1, w_span, ROMOrder);

% Calibration based on the maximal amplitude response of a validation FRC
idx_f_full = 1;
amp_max = max(FRC_full.(['F' num2str(idx_f_full)]).Amp);
[~,pos] = min(abs(amp-amp_max));
rho_max = rho_plot(pos);
f_red = abs(damp(rho_max)*rho_max) * f_full/f_full(idx_f_full);
ratio_force_data_full = f_red/f_full(idx_f_full);

% Compute with data-driven model
ddROM = struct('Dim',SSMDim,'Param',SSMFunction,'CCtoNormal',Tinv,'ReducedDynNormal',N,'CCfromNormal',T);
FRC_data = getFRC_ddROM(ddROM,f_red,w_span,1);
save FRC_data FRC_data

%%
% Plot
figure; hold on; grid on; box on; colors = colororder;
plot(frq, amp,'k','DisplayName', 'Backbone Curve - SSMLearn')% $\mathcal{O}(7)$
for ii = 1:length(f_full)
   freq_i = FRC_full.(['F' num2str(ii)]).Freq;
   amp_i  = FRC_full.(['F' num2str(ii)]).Amp;  
   stab_i = FRC_full.(['F' num2str(ii)]).Stab;
   [~,pos] = find(abs(diff(stab_i))==1);
   if isempty(pos); pos=length(stab_i); end;
   h_i = plot(freq_i(1:pos(1)),amp_i(1:pos(1)),'Color',colors(1,:),'Linewidth',2,...
        'DisplayName', 'FRC - stable -  SSMTool');
   if length(pos)>1
       h_ii = plot(freq_i(pos(1)+1:pos(2)),amp_i(pos(1)+1:pos(2)),'--','Color',colors(1,:),'Linewidth',2,...
        'DisplayName', 'FRC - unstable -  SSMTool');
   h_iii = plot(freq_i(pos(2)+1:end),amp_i(pos(2)+1:end),'Color',colors(1,:),'Linewidth',2);
   h_iii.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
   
   else
   h_ii = plot(freq_i(pos(1)+1:end),amp_i(pos(1)+1:end),'--','Color',colors(1,:),'Linewidth',2,...
        'DisplayName', 'FRC - unstable -  SSMTool');
   end
   if ii~= 1
       if ~isempty(h_i) h_i.Annotation.LegendInformation.IconDisplayStyle = 'off'; end
       if ~isempty(h_ii) h_ii.Annotation.LegendInformation.IconDisplayStyle = 'off'; end
   end
end
for ii = 1:length(f_red)
   freq_i = FRC_data.(['F' num2str(ii)]).Freq;
   amp_i  = FRC_data.(['F' num2str(ii)]).Amp;  
   stab_i = FRC_data.(['F' num2str(ii)]).Stab;
   [~,pos] = find(abs(diff(stab_i))==1);
   if isempty(pos); pos=length(stab_i); end;
   h_i = plot(freq_i(1:pos(1)),amp_i(1:pos(1)),'Color',colors(2,:),'Linewidth',2,...
        'DisplayName', 'FRC - stable -  SSMLearn');
   if length(pos)>1
   h_ii = plot(freq_i(pos(1)+1:pos(2)),amp_i(pos(1)+1:pos(2)),'--','Color',colors(2,:),'Linewidth',2,...
        'DisplayName', 'FRC - unstable -  SSMLearn');
   h_iii = plot(freq_i(pos(2)+1:end),amp_i(pos(2)+1:end),'Color',colors(2,:),'Linewidth',2);
   h_iii.Annotation.LegendInformation.IconDisplayStyle = 'off';
   else
   h_ii = plot(freq_i(pos(1)+1:end),amp_i(pos(1)+1:end),'--','Color',colors(2,:),'Linewidth',2,...
        'DisplayName', 'FRC - unstable -  SSMLearn');
   end
   if ii~= 1
       if ~isempty(h_i) h_i.Annotation.LegendInformation.IconDisplayStyle = 'off'; end
       if ~isempty(h_ii) h_ii.Annotation.LegendInformation.IconDisplayStyle = 'off'; end
   end
end
xlabel('$\Omega \, [$rad/s$]$','Interpreter','latex')
ylabel('$u \, [$m$]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
legend('interpreter','latex','location','NW')
xlim(w_span)
% ylim([0 0.03])
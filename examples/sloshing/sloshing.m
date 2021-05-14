%% Finding a 2D SSM from sloshing data
% 

clearvars
close all
clc

%% Example setup

load Data_decay_Omega=0.999_2018-04-03_17-57-11.csv
rawData = Data_decay_Omega_0_999_2018_04_03_17_57_11; clear Data_decay_Omega_0_999_2018_04_03_17_57_11;
%%
width = 500;
xData{1,1} = rawData(:,1)';
xData{1,2} = 100*rawData(:,3)'/width;
xData{1,1} = xData{1,1} - xData{1,1}(1);
tEnd = xData{1,1}(end);
nSamp = length(xData{1,1});
dt = tEnd/(nSamp-1);
xData{1,1} = 0:dt:tEnd;

nTraj = 1;
indTest = [1];
indTrain = 1;


%% Delay embedding

SSMDim = 2;
overEmbed = 0;
yData = coordinates_embedding(xData, SSMDim, 'OverEmbedding', overEmbed, 'ShiftSteps', 1);
showSpectrogram(yData(indTrain,:), 1);
%% 

sliceInt = [3.5, 60];
yDataTrunc = sliceTrajectories(yData, sliceInt);
%% Datadriven manifold fitting

SSMOrder = 1;
[V, SSMFunction, mfdInfo] = IMparametrization(yDataTrunc(indTrain,:), SSMDim, SSMOrder);
%% Plot and validation

etaData = getProjectedTrajs(yData, V);
etaDataTrunc = getProjectedTrajs(yDataTrunc, V);
plotReducedCoords(etaData);
plotSSMWithTrajectories(yData(indTrain,:), SSMFunction, 1, V, 10, 'SSMDimension', SSMDim)
view(-100,20); zlabel('$x_c$','Interpreter','latex')
%% Reduced order model

ROMOrder = 3;
[~,Tinv,N,T,NormalFormInfo] = IMdynamics_flow(etaDataTrunc(indTrain,:),'R_PolyOrd',ROMOrder,'style', 'normalform');

zData = transformComplex(Tinv, etaData);
zDataTrunc = transformComplex(Tinv, etaDataTrunc);
[zRec, yRec] = integrateFlows(N, zDataTrunc, @(q) SSMFunction(T(q)));
etaRec = transformComplex(T, zRec);
%%
[reducedTrajDist, fullTrajDist] = computeRecDynErrors(etaRec, yRec, etaDataTrunc, yDataTrunc);
RRMSE_normal = mean(fullTrajDist(indTest))

plotReducedCoords(etaDataTrunc(indTest(1),:), etaRec(indTest(1),:))
legend({'Test set (truncated)', 'Prediction'})

plotReconstructedTrajectory(yData(indTest(1),:), yRec(indTest(1),:), 1, 'm')
legend({'Test set', 'Prediction'}); ylabel('$u \, [$m$]$','Interpreter','latex')

normalFormEigenvalues = computeEigenvaluesFlow(NormalFormInfo)
lambda = normalFormEigenvalues;
%% Backbone curves

N_info = NormalFormInfo.N;
[damp, freq] = polarnormalform(N_info.coeff, N_info.exponents, N_info.phi);
figure
rhoCal = abs(zData{indTest(1),2}(1,1));
[dmp, frq, amp, rho_plot] = backbonecurves(damp, freq, SSMFunction, T, 1, rhoCal,'norm');
subplot(121); ylabel('$u \, [$m$]$','Interpreter','latex')
subplot(122); ylabel('$u \, [$m$]$','Interpreter','latex')

%%
% Calibration based on the maximal amplitude response of a validation FRC
% idx_f_full = 1;
% amp_max = max(FRC_full.(['F' num2str(idx_f_full)]).Amp);
% [~,pos] = min(abs(amp-amp_max));
% rho_max = rho_plot(pos);
% f_red = abs(damp(rho_max)*rho_max) * f_full/f_full(idx_f_full);
% ratio_force_data_full = f_red/f_full(idx_f_full)

% calibrationLoads = [0.000001];
% calibrationLoads = f_full(end);
% calloadvector = calibrationLoads.*f_vec;
% ICCal = getStaticResponse(K, M, F, calloadvector, 0, 0);
% uCal = observable(ICCal);
w_span = [0.77, 1.1]*7.8;
Omega = 0.91*7.8;
amplitudes = [0.64];
uCal = 6.2;
zCal = Tinv(V.'*repmat(uCal,size(yData{1,2},1),1));
rhoCal = abs(zCal(1));

fpsi = fsolve(@(fpsi) [damp(rhoCal)*rhoCal + fpsi(1)*sin(fpsi(2)); ...
    freq(rhoCal) - Omega + fpsi(1)/rhoCal*cos(fpsi(2))],...
    [0; 0]);
ratio = 2*0.5*abs(fpsi(1))/amplitudes(1);
% 0.5*imag(lambda(1))*abs(NormalFormInfo.iT.lintransf*sum(V).')

f_red = ratio*amplitudes*0.43;


%% Compute with data-driven model
ddROM = struct('Dim', SSMDim, 'Param', SSMFunction, 'CCtoNormal', Tinv, ...
    'ReducedDynNormal', N, 'CCfromNormal', T);
FRC_data = getFRC_ddROM(ddROM,f_red,w_span,1);

FRC_data.F1.Freq = FRC_data.F1.Freq/7.8;
frq = frq/7.8;

% Plot
figure(100); hold on;
plot(frq, amp,'k','DisplayName', 'Backbone - SSMlearn')
plot_FRC(FRC_data, 'b', 'SSMLearn')

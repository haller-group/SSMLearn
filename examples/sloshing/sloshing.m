%% Finding a 2D SSM from sloshing data
% 

clearvars
close all
clc

%% Example setup
datadir = 'decaydata/';
expAmp{1} = load('FRC_experimental_Amplitude0.09.csv');
expAmp{2} = load('FRC_experimental_Amplitude0.17.csv');
expAmp{3} = load('FRC_experimental_Amplitude0.32.csv');
expAmp{4} = load('FRC_experimental_Amplitude0.64.csv');
expPhase{1} = load('FRC_experimental_Phase0.09.csv');
expPhase{2} = load('FRC_experimental_Phase0.17.csv');
expPhase{3} = load('FRC_experimental_Phase0.32.csv');
expPhase{4} = load('FRC_experimental_Phase0.64.csv');
amplitudes = [0.09 0.17, 0.32, 0.64];
csvfiles = dir([datadir,'*.csv']);
ii = 0;
for file = csvfiles'
    ii = ii + 1;
    rawData{ii} = load([datadir, file.name]);
end
%%
width = 500;
nTraj = numel(rawData);
rawColInds = 3;

for iTraj = 1:nTraj
    xData{iTraj,1} = rawData{iTraj}(:,1)';
    xData{iTraj,2} = sgolayfilt(100*rawData{iTraj}(:,rawColInds)/width,7,29)';
%     xData{iTraj,2} = 100*rawData{iTraj}(:,rawColInds)'/width;
    xData{iTraj,1} = xData{iTraj,1} - xData{iTraj,1}(1);
    tEnd = xData{iTraj,1}(end);
    nSamp = length(xData{iTraj,1});
    dt = tEnd/(nSamp-1);
    xData{iTraj,1} = 0:dt:tEnd;
end
% clear rawData

indTrain = [17 19];
% indTrain = 1:nTraj;
indTest = indTrain;


%% Delay embedding

% showSpectrogram(xData(indTrain,:));
% for iTraj = 1:nTraj
%     xDataPass{iTraj,1} = xData{iTraj,1};
%     xDataPass{iTraj,2} = bandpass(xData{iTraj,2}, [0.9,1.8], 1/dt);
% end
% showSpectrogram(xDataPass(indTrain,:));
SSMDim = 2;
overEmbed = 0;
yData = coordinates_embedding(xData, SSMDim, 'OverEmbedding', overEmbed, 'ShiftSteps', 1);
%% 

sliceInt = [1.9, Inf];
yDataTrunc = sliceTrajectories(yData, sliceInt);
% Datadriven manifold fitting

SSMOrder = 1;
[V, SSMFunction, mfdInfo] = IMparametrization(yDataTrunc(indTrain,:), SSMDim, SSMOrder);
% Plot and validation

etaData = getProjectedTrajs(yData, V);
etaDataTrunc = getProjectedTrajs(yDataTrunc, V);
plotReducedCoords(etaData);
plotSSMWithTrajectories(yDataTrunc(indTrain,:), SSMFunction, 1, V, 10, 'SSMDimension', SSMDim)
%% Reduced order model

ROMOrder = 3;
% [~,Tinv,N,T,NormalFormInfo] = IMdynamics_flow(etaDataTrunc(indTrain,:),...
%     'R_PolyOrd',ROMOrder,'style', 'normalform', 'l_vals', [0.1]);
[~,Tinv,N,T,NormalFormInfo] = IMdynamics_flow(etaDataTrunc(indTrain,:),...
    'R_PolyOrd',ROMOrder,'style', 'normalform');

zData = transformComplex(Tinv, etaData);
zDataTrunc = transformComplex(Tinv, etaDataTrunc);
[zRec, yRec] = integrateFlows(N, zDataTrunc, @(q) SSMFunction(T(q)));
etaRec = transformComplex(T, zRec);
%%
[reducedTrajDist, fullTrajDist, fullAmpError] = computeRecDynErrors(etaRec, yRec, etaDataTrunc, yDataTrunc);
RRMSE_normal = mean(fullTrajDist(indTest))
%
plotReducedCoords(etaDataTrunc(indTest(1),:))%, etaRec(indTest(1),:))
legend({'Test set (truncated)', 'Prediction'});

ppw = length(indTest);
indPlots = indTest;
for ii = 0:floor(length(indPlots)/ppw-1)
    inds = ppw*ii+1:min(length(indPlots),ppw*ii+ppw);
    plotReconstructedTrajectory(yData(indPlots(inds),:), yRec(indPlots(inds),:), 1, 'm')
    ylabel('$\hat{Y} \, [\%]$','Interpreter','latex'); title('');
end

normalFormEigenvalues = computeEigenvaluesFlow(NormalFormInfo)
lambda = normalFormEigenvalues;
%% Backbone curves

N_info = NormalFormInfo.N;
% N_info.coeff = [-0.055 + 7.81i  -0.09 - 2i]
[damp, freq] = polarnormalform(N_info.coeff, N_info.exponents, N_info.phi);
figure
rhoMax = 1.*max(abs(zDataTrunc{17,2}(1,:)));
[dmp, frq, amp, rho_plot] = backbonecurves(damp, freq, SSMFunction, T, 1, rhoMax);
subplot(121); ylabel('$\hat{X} \, [\%]$','Interpreter','latex')
subplot(122); ylabel('$\hat{X} \, [\%]$','Interpreter','latex')
frq = frq/7.8;

%%
w_span = [0.77, 1.06]*7.8;
for iAmp = 1:length(amplitudes)
%     [Omega, pos] = max(expAmp{iAmp}(:,1));
%     Omega = Omega*7.8;
%     uCal = expAmp{iAmp}(pos,2);
    [uCal, pos] = max(expAmp{iAmp}(:,2));
    if iAmp == 4; [uCal, pos] = min(expAmp{iAmp}(:,2)); end
    Omega = expAmp{iAmp}(pos,1)*7.8;
    yCal = uCal*cos(Omega*dt*((1:size(V,1))-ceil(0.5*size(V,1)))).';
    zCal = Tinv(V.'*yCal);
    rhoCal = abs(zCal(1));
    f_red(iAmp) = sqrt(rhoCal^2*(freq(rhoCal)-Omega)^2+rhoCal^2*(damp(rhoCal))^2);
end
% f_red = amplitudes*f_red(iAmp)/amplitudes(iAmp)

% Compute FRC analytically
options = optimoptions('fsolve', 'MaxFunctionEvaluations', 10000);
for iAmp = 1:length(amplitudes)
    rhoTip = abs(fsolve(@(rho) 1e5*(f_red(iAmp)-(rho*damp(rho))), rhoCal, options));
    rhoFRC = logspace(log10(rhoTip*0.03), log10(rhoTip),1000);
    rhoFRC = [rhoFRC, -fliplr(rhoFRC)];
    OmegaFRC(iAmp,:) = 1/7.8*real(freq(rhoFRC) + -1./rhoFRC.*sqrt(f_red(iAmp)^2-(rhoFRC.*damp(rhoFRC)).^2));
    uFRC(iAmp,:) = max(abs(SSMFunction(T([rhoFRC;rhoFRC]))));
    psiFRC(iAmp,:) = acos((7.8*OmegaFRC(iAmp,:) - freq(rhoFRC)).*rhoFRC.*sign(rhoFRC)./f_red(iAmp));
end

%% Plot
figure(100); hold on; colors = colororder;
plot(frq, amp,'k','DisplayName', 'Backbone - SSMlearn')
figure(101); hold on; colors = colororder;
% plot_FRC(FRC_data, 'b', 'SSMLearn')
for iAmp = 1:length(amplitudes)
    figure(100);
    plot(OmegaFRC(iAmp,:), uFRC(iAmp,:), 'LineWidth', 2, 'Color', colors(iAmp+1,:),...
        'DisplayName', ['SSMLearn A = ', num2str(amplitudes(iAmp)), ' %'])
    plot(expAmp{iAmp}(:,1), expAmp{iAmp}(:,2), '.', 'MarkerSize', 12, ...
        'Color', colors(iAmp+1,:), 'DisplayName', ['Exp. A = ', num2str(amplitudes(iAmp)), ' %'],...
        'MarkerSize', 14)
    figure(101);
    plot(OmegaFRC(iAmp,:), 180*(psiFRC(iAmp,:)-pi)/pi, 'LineWidth', 2, 'Color', colors(iAmp+1,:),...
        'DisplayName', ['SSMLearn A = ', num2str(amplitudes(iAmp)), ' %'])
    plot(expPhase{iAmp}(:,1), expPhase{iAmp}(:,2), '.', 'MarkerSize', 12, ...
        'Color', colors(iAmp+1,:), 'DisplayName', ['Exp. A = ', num2str(amplitudes(iAmp)), ' %'])
end
figure(100)
xlim(w_span/7.8)
ylim([0,10]); %max(vertcat(expAmp{1,:})*[0;1])*1.4]
xlabel('Excitation frequency $\Omega$ [normalized]', 'Interpreter', 'latex')
ylabel('Amplitude $\hat{X}$ (\%)', 'Interpreter', 'latex')
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
legend('location', 'SW')
p = 1;
if p ~=1
    plot(freq(linspace(0,rhoMax))*p/7.8, ...
        max(abs(SSMFunction(T([linspace(0,rhoMax);linspace(0,rhoMax)])))), ...
        'b', 'DisplayName', [num2str(p),'\Omega'])
    xlim([min(p*w_span(1),w_span(1)),max(p*w_span(2),w_span(2))]/7.8)
end
figure(101)
xlim(w_span/7.8)
ylim([-180,0]); %max(vertcat(expAmp{1,:})*[0;1])*1.4]
xlabel('Excitation frequency $\Omega$ [normalized]', 'Interpreter', 'latex')
ylabel('Phase difference $\phi$', 'Interpreter', 'latex')
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
legend('location', 'SW')
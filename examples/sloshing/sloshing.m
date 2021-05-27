%% Finding a 2D SSM from sloshing data
% 

clearvars
close all
clc

%% Example setup
datadir = 'decaydata/';
expAmpC{1} = load('FRC_experimental_Amplitude0.09.csv');
expAmpC{2} = load('FRC_experimental_Amplitude0.17.csv');
expAmpC{3} = load('FRC_experimental_Amplitude0.32.csv');
expAmpC{4} = load('FRC_experimental_Amplitude0.64.csv');
expPhase{1} = load('FRC_experimental_Phase0.09.csv');
expPhase{2} = load('FRC_experimental_Phase0.17.csv');
expPhase{3} = load('FRC_experimental_Phase0.32.csv');
expPhase{4} = load('FRC_experimental_Phase0.64.csv');
expWall{1} = load('FRC_experimental_Wall0.09.csv');
expWall{2} = load('FRC_experimental_Wall0.17.csv');
expWall{3} = load('FRC_experimental_Wall0.32.csv');
expWall{4} = load('FRC_experimental_Wall0.64.csv');
amplitudes = [0.09 0.17, 0.32, 0.64];
csvfiles = dir([datadir,'*.csv']);
ii = 0;
rawData = cell(length(csvfiles),1);
for file = csvfiles'
    ii = ii + 1;
    rawData{ii} = load([datadir, file.name]);
    decayFrequencies(ii) = str2num(file.name(18:22));
    titles{ii} = ['$\Omega = ', num2str(decayFrequencies(ii)), '$'];
end
%%
width = 500;
nTraj = numel(rawData);
% rawColInds = [3];
% rawColInds = [3 5:250:1500];
rawColInds = [5];
if rawColInds(1) == 5
    expAmp = expWall;
    for ii=1:4
        expAmp{ii}(:,2) = 100*expAmp{ii}(:,2);
    end
else
    expAmp = expAmpC;
end

for iTraj = 1:nTraj
    cutoffPoint = find(abs(diff(rawData{iTraj}(1:100,2)')) + ...
        [abs(diff(rawData{iTraj}(1:100,2)',2)),0] < 0.01, 1, 'first');
    xData{iTraj,1} = rawData{iTraj}(cutoffPoint:end,1)';
    xData{iTraj,2} = sgolayfilt(100*rawData{iTraj}(cutoffPoint:end,rawColInds)/width,7,29)';
%     xData{iTraj,2} = 100*rawData{iTraj}(cutoffPoint:end,rawColInds)'/width;
    xData{iTraj,1} = xData{iTraj,1} - xData{iTraj,1}(1);
    tEnd = xData{iTraj,1}(end);
    nSamp = length(xData{iTraj,1});
    dt = tEnd/(nSamp-1);
    xData{iTraj,1} = 0:dt:tEnd;
end
% clear rawData

indTrain = 10:27;
% indTrain = 1:nTraj;
indTest = indTrain;


%% Delay embedding

SSMDim = 2;
overEmbed = 0;
if length(rawColInds) > 1; overEmbed = 4; end
[yData, opts_embd] = coordinates_embedding(xData, SSMDim, 'OverEmbedding', overEmbed, 'ShiftSteps', 1);
embedDim = size(yData{1,2},1)/length(rawColInds);
outdof = floor(embedDim/2)*length(rawColInds)+1;
%% 

sliceInt = [0.0, Inf];
yDataTrunc = sliceTrajectories(yData, sliceInt);

SSMOrder = 1;
[V, SSMFunction, mfdInfo] = IMparametrization(yDataTrunc(indTrain,:), SSMDim, SSMOrder);
%% Plot and validation

etaData = getProjectedTrajs(yData, V);
etaDataTrunc = getProjectedTrajs(yDataTrunc, V);
plotReducedCoords(etaDataTrunc);
plotSSMWithTrajectories(yDataTrunc(indTrain,:), SSMFunction, 1, V, 10, 'SSMDimension', SSMDim)
%% Reduced order model

% ROMOrder = 7;
% [~,Tinv,N,T,NormalFormInfo] = IMdynamics_flow(etaDataTrunc(indTrain,:),...
%     'R_PolyOrd',ROMOrder,'style', 'normalform', 'l_vals', [0.1]);
ROMOrder = 3;
[~,Tinv,N,T,NormalFormInfo] = IMdynamics_flow(etaDataTrunc(indTrain,:), ...
    'R_PolyOrd', ROMOrder, 'style', 'normalform');
% [~,Tinv,N,T,NormalFormInfo,ROMOrder,ROMerrs] = optimizeDynamicsFlow(etaDataTrunc(indTrain,:));

zData = transformComplex(Tinv, etaData);
zDataTrunc = transformComplex(Tinv, etaDataTrunc);
[zRec, yRec] = integrateFlows(N, zDataTrunc, @(q) SSMFunction(T(q)));
etaRec = transformComplex(T, zRec);
%%
[reducedTrajDist, fullTrajDist, fullAmpError] = computeRecDynErrors(etaRec, yRec, etaDataTrunc, yDataTrunc);
RRMSE_normal = mean(fullTrajDist(indTest))
%
plotReducedCoords(etaDataTrunc(indTest,:))%, etaRec(indTest(1),:))
legend({'Test set (truncated)', 'Prediction'});

indPlots = indTrain;
ppw = 9;%length(indTrain);
for ii = 0:floor(length(indPlots)/ppw-1)
    inds = ppw*ii+1:min(length(indPlots),ppw*ii+ppw);
    plotReconstructedTrajectory(yData(indPlots(inds),:), yRec(indPlots(inds),:), outdof, 'm', ...
        titles(indPlots(inds)))
    ylabel('$\hat{X} \, [\%]$','Interpreter','latex'); title('');
end
s=findobj('type','legend');
delete(s)
% savepic(['figs/',num2str(indTrain)])
%%
% plotTrajs(xData(indTrain,:), 'styles', {'-','-'}, 'colors', {'r','k'}, 'legendnames', titles(indTrain))
% savepic(['figs/',num2str(indTrain)])
%%
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
clear yCal
w_span = [0.77, 1.06]*7.8;
yObservable = @(y) abs(y(outdof,:));
for iAmp = 1:length(amplitudes)
    [uCal, pos] = min(expAmp{iAmp}(:,2));
    if iAmp == 4; [uCal, pos] = min(expAmp{iAmp}(:,2)); end
    Omega(iAmp) = expAmp{iAmp}(pos,1)*7.8;
    yCal(:,iAmp) = uCal*V(:,1)./V(outdof,1);
end
f_red = calibrateFRC(yCal, Omega, V, Tinv, damp, freq);

FRC_data = computeFRC(f_red, damp, freq, SSMFunction, T, yObservable);
for iAmp = 1:length(amplitudes)
    FRC_data.(['F',num2str(iAmp)]).Freq = FRC_data.(['F',num2str(iAmp)]).Freq/7.8;
    FRC_data.(['F',num2str(iAmp)]).Nf_Phs = -180+180/pi*FRC_data.(['F',num2str(iAmp)]).Nf_Phs;
end

%% Plot
figure(100); hold on; colors = colororder; colors = colors(2:end,:);
plot(frq, amp,'k','DisplayName', 'Backbone - SSMlearn')
figure(101); hold on;
for iAmp = 1:length(amplitudes)
    labels{iAmp} = ['A = ', num2str(amplitudes(iAmp)), ' %'];
    figure(100);
    plot(expAmp{iAmp}(:,1), expAmp{iAmp}(:,2), '.', 'MarkerSize', 12, ...
        'Color', colors(iAmp,:), 'DisplayName', ['Exp. A = ', num2str(amplitudes(iAmp)), ' %'],...
        'MarkerSize', 14)
    figure(101);
    plot(expPhase{iAmp}(:,1), expPhase{iAmp}(:,2), '.', 'MarkerSize', 12, ...
        'Color', colors(iAmp,:), 'DisplayName', ['Exp. A = ', num2str(amplitudes(iAmp)), ' %'])
end
figure(100)
plotFRC(FRC_data, colors, labels)
xlim(w_span/7.8)
ylim([0,50]);
xlabel('Excitation frequency $\Omega$ [normalized]', 'Interpreter', 'latex')
ylabel('Amplitude $\hat{X}$ (\%)', 'Interpreter', 'latex')
p = 1;
if p ~=1
    plot(freq(linspace(0,rhoMax))*p/7.8, ...
        max(abs(SSMFunction(T([linspace(0,rhoMax);linspace(0,rhoMax)])))), ...
        'b', 'DisplayName', [num2str(p),'\Omega'])
    xlim([min(p*w_span(1),w_span(1)),max(p*w_span(2),w_span(2))]/7.8)
end
figure(101)
plotFRC(FRC_data, colors, labels, 'y', 'Phase')
xlim(w_span/7.8)
ylim([-180,0]);
xlabel('Excitation frequency $\Omega$ [normalized]', 'Interpreter', 'latex')
ylabel('Phase difference $\phi$', 'Interpreter', 'latex')
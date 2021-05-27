%% Finding a 2D SSM from sloshing data
% 

clearvars
close all
clc

%% Example setup
datadir = 'decaydata/';
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
rawColInds = [3];
% rawColInds = [5:1:1500];
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

% F = @(t,x) [0,1;-7.8^2,-0.1]*x - [0;1]*2.*x(1).^3 - [0;1]*0.0003.*x(2).^3;
% xData(17,:) = integrateTrajectories(F, @(x)x(1,:), 80, 80/dt+1, 1, [4;0])

indTrain = [17];
% indTrain = 1:nTraj;
indTest = indTrain;


%% Delay embedding

SSMDim = 2;
overEmbed = 199;
if length(rawColInds) > 1; overEmbed = 4; end
[yData, opts_embd] = coordinates_embedding(xData(indTrain,:), SSMDim, 'OverEmbedding', overEmbed, 'ShiftSteps', 1);
embedDim = size(yData{1,2},1)/length(rawColInds);
outdof = floor(embedDim/2)*length(rawColInds)+1;
%% 
% plotTrajs(xData(indTrain,:), 'styles', {'-','-'}, 'colors', {'r','k'}, 'legendnames', titles(indTrain))

sliceInt = [0.0, Inf];
xDataTrunc = sliceTrajectories(xData, sliceInt);
yDataTrunc = sliceTrajectories(yData, sliceInt);
X1 = []; X2 = [];
indTrainY = 1:length(indTrain);
for iTraj = indTrainY
    X1 = [X1 yDataTrunc{iTraj,2}(:,1:end-1)];
    X2 = [X2 yDataTrunc{iTraj,2}(:,2:end)];
end

%%
figdir = 'figs/DMD/';
problemtype = 'ERRRR';
if isscalar(rawColInds) && rawColInds(1) == 3
    problemtype = 'xcenter'
elseif rawColInds(1) == 5
    problemtype = 'surface'
end
r = 8;
[Phi,omega,lambda,b,Xdmd,t] = DMD(X1,X2,r,dt);
figure(1)
plot(t(1:end/length(indTrain)),Xdmd(1,1:end/length(indTrain)), 'DisplayName', 'DMD reconstruction')
hold on
plot(xDataTrunc{indTrain(1),1}, xDataTrunc{indTrain(1),2}(1,:), 'DisplayName', 'Experiment')
legend
xlabel('$t$', 'Interpreter', 'latex')
ylabel('$\hat{X}$ [\%]', 'Interpreter', 'latex')
savepic([figdir,'DMDReconstruction_Delay',num2str(embedDim),'_', problemtype])
%
if strcmp(problemtype, 'xcenter')
modes = Phi\X1;
plotModes = 4;
dim = size(Phi,1);
yscale = max(max(abs(Phi)));
for iMode = 1:plotModes
    figure(2)
    subplot(plotModes,1,iMode)
    plot(linspace(0,1,dim),real(Phi(1:dim,2*iMode)), 'LineWidth', 2)
    ylabel(['Real($\Phi_',num2str(iMode),'$)'], 'Interpreter', 'latex')
    ylim([-yscale,yscale])
    set(gca,'fontsize',10)
    figure(3)
    subplot(plotModes,1,iMode)
    plot(linspace(0,1,dim),imag(Phi(1:dim,2*iMode)), 'LineWidth', 2)
    ylabel(['Imag($\Phi_',num2str(iMode),'$)'], 'Interpreter', 'latex')
    ylim([-yscale,yscale])
    set(gca,'fontsize',10)
end
figure(2)
xlabel('t/d', 'Interpreter', 'latex')
savepic([figdir,'realPhiDelay',num2str(embedDim),'_', problemtype])
figure(3)
xlabel('t/d', 'Interpreter', 'latex')
savepic([figdir,'imagPhiDelay',num2str(embedDim),'_', problemtype])
save([figdir,'omegaDelay',num2str(embedDim),'_', problemtype], 'omega')
elseif strcmp(problemtype, 'surface')
plotModes = 4;
dim = length(rawColInds);
yscale = max(max(abs(Phi)));
for iMode = 1:plotModes
    figure(4)
    subplot(plotModes,1,iMode)
    plot(linspace(0,1,dim),real(Phi(1:dim,2*iMode)), 'LineWidth', 2)
    ylabel(['Real($\Phi_',num2str(iMode),'$)'], 'Interpreter', 'latex')
    ylim([-yscale,yscale])
    set(gca,'fontsize',10)
    figure(5)
    subplot(plotModes,1,iMode)
    plot(linspace(0,1,dim),imag(Phi(1:dim,2*iMode)), 'LineWidth', 2)
    ylabel(['Imag($\Phi_',num2str(iMode),'$)'], 'Interpreter', 'latex')
    ylim([-yscale,yscale])
    set(gca,'fontsize',10)
end
figure(4)
xlabel('x/w', 'Interpreter', 'latex')
savepic([figdir,'realPhiDelay',num2str(embedDim),'_', problemtype])
figure(5)
xlabel('x/w', 'Interpreter', 'latex')
savepic([figdir,'imagPhiDelay',num2str(embedDim),'_', problemtype])
save([figdir,'omegaDelay',num2str(embedDim),'_', problemtype], 'omega')
end
%%
[u,s,v] = svd(X1,'econ');
figure, plot(diag(s)/sum(diag(s)),'o');
figure, plot(u(:,1:4))
figure, plot(v(:,5))
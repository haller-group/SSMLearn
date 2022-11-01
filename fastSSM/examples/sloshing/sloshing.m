%% Finding a 2D SSM from sloshing data
% 

clearvars
close all

%% Example setup
frcdir = 'frcdata/';
expAmp{1} = load([frcdir,'Measurements_A=0.09%.txt']);
expAmp{2} = load([frcdir,'Measurements_A=0.17%.txt']);
expAmp{3} = load([frcdir,'Measurements_A=0.32%.txt']);
amplitudes = [0.09, 0.17, 0.32];
load decaydata
decayFrequencies = [0.960,0.967];
width = 500;

%% Delay embedding
SSMDim = 2;
overEmbed = 9;
shiftSteps = 1;
dt = (xData{2,1}(end) - xData{2,1}(1))/(length(xData{2,1})-1);
if dt < 0.03; shiftSteps = 3; end
[yData, opts_embd] = coordinatesEmbedding(xData, SSMDim, 'OverEmbedding', overEmbed, 'ShiftSteps', shiftSteps);
embedDim = size(yData{1,2},1)/size(xData{1,2}, 1);
outdof = floor(embedDim/2)*size(xData{1,2}, 1)+1;

%% Train model
indTest = 1;
indTrain = 2;

SSMOrder = 5;

% [Mmap, iMmap, Tmap, iTmap, Nflow, yRecF, BBC] = fastSSMplus(yData(indTrain,:), SSMDim, SSMOrder, 5, 5);
[Mmap, iMmap, Tmap, iTmap, Nflow, yRecF, BBC] = fastSSM(yData(indTrain,:), SSMOrder);

zDataTrunc = transformTrajectories(@(y) iTmap(iMmap(y)), yData);
zRec = integrateFlows(@(z) Nflow(0,z), zDataTrunc);
yRec = transformTrajectories(@(z) Mmap(Tmap(z)), zRec);

NMTE = computeTrajectoryErrors(yRec, yData, [1,2:embedDim:1535])

%% Evaluate model
customFigure();
plot(yData{indTest(1),1}(1,:), yData{indTest(1),2}(outdof+1,:), 'Color', [0,0,0], 'LineWidth', 1.6, 'DisplayName', 'Original');
plot(yRec{indTest(1),1}(1,:), yRec{indTest(1),2}(outdof+1,:), '--', 'Color', [0.2,0.9,0.2], 'LineWidth', 1.6, 'DisplayName', 'Reconstruction');
xlabel('time [s]', 'Interpreter', 'latex');
ylabel('$h_{-w/2}$ [\%]', 'Interpreter', 'latex');
legend

IMInfo = struct();
IMInfo.parametrization.map = Mmap;
IMInfo.parametrization.dimension = SSMDim;
IMInfo.chart.map = iMmap;
plotSSMWithTrajectories(IMInfo, outdof+1, yData(indTrain,:), 'Colors', [0,0,0], 'ColorSurf', 'r')
xlabel('$\xi_1$', 'Interpreter', 'latex'); ylabel('$\xi_2$', 'Interpreter', 'latex');
zlabel('$h_{-w/2}$ [\%]', 'Interpreter', 'latex'); view(33,13)

%% FRC
clear yCal Omega
ampF = @(y) abs(y(outdof,:));
w_span = [0.87, 1.06]*7.8;
V = [Mmap(1e-5*[1;0]),Mmap(1e-5*[0;1])]/1e-5;
[~,indV] = max(ampF(V));
for iAmp = 1:length(amplitudes)
    [uCal, pos] = max(expAmp{iAmp}(:,2));
    if iAmp == 4; [uCal, pos] = min(expAmp{iAmp}(:,2)); end
    Omega(iAmp,1) = expAmp{iAmp}(pos,1)*7.8;
    yCal(:,iAmp) = uCal*V(:,indV)./V(outdof,indV);
end

[fRed, OmegaFm, OmegaFp, uF, psiF] = miniComputeFRC(yCal, Omega, BBC, Mmap, iMmap, Tmap, iTmap, ampF);
%
colors = get(0, 'DefaultAxesColorOrder'); colors = colors(2:end,:);
fig1 = customFigure();
fig2 = customFigure();
for iAmp = 1:length(amplitudes)
    rhoMax = fsolve(@(rho) abs(OmegaFm(rho,fRed(iAmp))-BBC.freq(rho)),0.01);
    rho = linspace(0,rhoMax,100);
    realinds = imag([OmegaFm(rho,fRed(iAmp)), OmegaFp(rho,fRed(iAmp))]) < 1e-3;
    pr = @(x) real(x(:,realinds));
    figure(fig1);
    plot(pr([OmegaFm(rho,fRed(iAmp)), OmegaFp(rho,fRed(iAmp))]), ...
        pr([uF(rho),uF(rho)]), 'Color', colors(iAmp,:), 'Linewidth', 2, ...
        'DisplayName', ['\texttt{fastSSM} $A = ', num2str(amplitudes(iAmp)), ' \%$']);
    plot(expAmp{iAmp}(:,1)*7.8, expAmp{iAmp}(:,2), '.', 'MarkerSize', 16, ...
        'Color', colors(iAmp,:), 'DisplayName', ['Exp. $A = ', num2str(amplitudes(iAmp)), ' \%$'],...
        'MarkerSize', 14)
    figure(fig2);
    plot(pr([OmegaFm(rho,fRed(iAmp)), OmegaFp(rho,fRed(iAmp))]), ...
        -180+180/pi*pr([psiF([rho,rho],fRed(iAmp),[OmegaFm(rho,fRed(iAmp)), OmegaFp(rho,fRed(iAmp))]),...
        psiF([rho,rho],fRed(iAmp),[OmegaFm(rho,fRed(iAmp)), OmegaFp(rho,fRed(iAmp))])]), ...
        'Color', colors(iAmp,:), 'Linewidth', 2, ...
        'DisplayName', ['\texttt{fastSSM} $A = ', num2str(amplitudes(iAmp)), ' \%$']);
    plot(expAmp{iAmp}(:,1)*7.8, expAmp{iAmp}(:,3), '.', 'MarkerSize', 16, ...
        'Color', colors(iAmp,:), 'DisplayName', ['Exp. $A = ', num2str(amplitudes(iAmp)), ' \%$'],...
        'MarkerSize', 14)
end
figure(fig1);
plot(BBC.frequency, BBC.amplitude, 'k', 'DisplayName', '\texttt{fastSSM} backbone curve')
xlim(w_span)
ylim([0,6]);
legend('Interpreter', 'latex', 'location', 'best')
xlabel('$\Omega$ [rad/s]', 'Interpreter', 'latex')
ylabel('$\hat{X}$ [\%]', 'Interpreter', 'latex')

figure(fig2);
xlim(w_span)
ylim([-180,50])
yticks([-180,-120,-60,0]);
legend('Interpreter', 'latex', 'location', 'north')
xlabel('$\Omega$ [rad/s]', 'Interpreter', 'latex')
ylabel('$\psi$ [$^\circ$]', 'Interpreter', 'latex')


%% Display surface profile
iFreq = indTest(1);
OmegaSim = 7.8*decayFrequencies(iFreq);
xCoords = linspace(0, width, size(xData{1,2}, 1)-1);
customFigure();
set(gca, 'fontsize', 28)
frames = [];
for iTime = 1:2:length(yRec{iFreq,1})/5
    cla
	set(gca,'Position',[0.2 0.2 0.7 0.7]);
    xlim([0,width])
    ylim([-130,250])
    area(xCoords, yData{iFreq,2}(2:size(xData{1,2}, 1),iTime)*width/100, ...
        -130, 'FaceColor', [0.5,0.8,1], 'DisplayName', 'Experiment');
    plot(xCoords, yRec{iFreq,2}(2:size(xData{1,2}, 1),iTime)*width/100, ...
        ':', 'LineWidth', 4, 'DisplayName', 'Prediction');
    xlabel('$x$ [mm]', 'interpreter', 'latex'); xticklabels([0,100,200,300,400]);
    ylabel('Elevation $h$ [mm]', 'interpreter', 'latex'); yticklabels('auto');
    legend('location', 'NW');
    text(300, 220, ['$t = ', num2str(round(yRec{1,1}(iTime)*1000)/1000), '$ s'], 'fontsize', 28, 'interpreter', 'latex')
    frames = [frames, getframe];
end
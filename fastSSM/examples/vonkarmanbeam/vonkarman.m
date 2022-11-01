%% Learning a 2D SSM with fastSSMplus from FE simulation data
% This is an example of how to reconstruct a slow 2D SSM of a mechanical system using synthetic measurements of a scalar quantity. 
% In this example, we measure the midpoint displacement of a clamped-clamped von Kármán beam. [1]
% We use fastSSMplus to predict backbone curves and forced response.
%
% [1] S. Jain, P. Tiso, and G. Haller, Exact nonlinear model reduction for a von Kármán beam: slow-fast decomposition and spectral 
% submanifolds, Journal of Sound and Vibration 423 (2018) 195-211. https://doi.org/10.1016/j.jsv.2018.01.049 
%

clearvars
close all
addpath tools

nElements = 12;
[M, C, K, fnl, fExt, obsdof, PlotFieldonDefMesh] = buildModel(nElements);

n = size(M,1);    % mechanical dofs (axial def, transverse def, angle)
[F, lambda] = functionFromTensors(M, C, K, fnl);
loads = -[12 14]*1e3;
nTraj = size(loads, 2);
indTest = 1;
indTrain = 2;
loadvector = loads.*fExt;
%%

newMeasurement = false;
observable = @(x) x(obsdof,:);
slowTimeScale = 2*pi/abs(lambda(1));
numberPeriods = 125; numberPointsPerPeriod = 100;
if newMeasurement
    IC = getStaticResponse(K, M, F, loadvector, false, PlotFieldonDefMesh);
    endTime = numberPeriods*slowTimeScale;
    nSamp = numberPeriods*numberPointsPerPeriod+1;
    dt = endTime/(nSamp-1);
    xData = integrateTrajectories(F, endTime, IC, nSamp, observable);
    DataInfo = struct('nElements', nElements, 'loadvector', loadvector);
    save('vonkarmandata/dataVKDecay.mat', 'DataInfo', 'xData', 'dt', 'endTime', 'nSamp', 'lambda')
else
    load vonkarmandata/dataVKDecay.mat
    if nElements ~= DataInfo.nElements
       error('The loaded data comes from a model with a different number of elements.')
    end
end

%%
SSMDim = 2;
overEmbed = 0;
yData = coordinatesEmbedding(xData, SSMDim, 'OverEmbedding', overEmbed);

slowTimeScale = 2*pi/abs(lambda(1));
sliceInt = [4*slowTimeScale, endTime];
yDataTrunc = sliceTrajectories(yData, sliceInt);
outdof = 3;

SSMOrder = 1;
ROMOrder = 5;
NFOrder = 11;


[Hmap, iHmap, Tmap, iTmap, Nflow, yRecF, BBC] = fastSSMplus(yDataTrunc(indTrain,:), 2, SSMOrder, ROMOrder, NFOrder);
yRec = {yDataTrunc{indTrain,1}, yRecF};
NMTE = computeTrajectoryErrors(yRec, yDataTrunc(indTrain,:))

%%
customFigure();
plot(yData{indTrain,1}(1,:), yData{indTrain,2}(outdof,:), 'Color', [0,0,0], 'LineWidth', 0.8, 'DisplayName', 'Training data');
plot(yRec{1,1}(1,:), yRec{1,2}(outdof,:), '--', 'Color', [0.1,0.9,0.1], 'LineWidth', 0.8, 'DisplayName', 'Reconstruction');
xlabel('time [s]', 'Interpreter', 'latex');
ylabel('$u$ [m]', 'Interpreter', 'latex');
xlim([0,sliceInt(2)])
legend

IMInfo = struct();
IMInfo.parametrization.map = Hmap;
IMInfo.parametrization.dimension = SSMDim;
IMInfo.chart.map = iHmap;
plotSSMWithTrajectories(IMInfo, outdof, yDataTrunc(indTrain,:), 'Colors', [0,0,0], 'ColorSurf', 'r')
xlabel('$\xi_1$', 'Interpreter', 'latex')
ylabel('$\xi_2$', 'Interpreter', 'latex')
zlabel('$u$', 'Interpreter', 'latex')
view(50, 35)

%% FRC
newSweep = false;
fFull = 75;
if newSweep
    fSweep = fFull(end)*(M\fExt);
    FRCSweep = integrateSweep(F, [650,800], n, fSweep, observable, ...
        'initRevolutions',150,'stepRevolutions', 150, 'step', 5);
    DataInfo = struct('nElements', nElements, 'fSweep', fSweep, ...
        'obsdof', obsdof);
    save('vonkarmandata/dataVKSweep.mat', 'DataInfo', 'FRCSweep')
else
    load('vonkarmandata/dataVKSweep.mat', 'DataInfo', 'FRCSweep')
end
%%
clear yCal Omega
omegaSpan = abs(imag(lambda(1)))*[0.9 1.2];
ampF = @(y) abs(y(outdof,:));
V = [Hmap(1e-10*[1;0]),Hmap(1e-10*[0;1])]/1e-10;
[~,indV] = max(ampF(V));
for iAmp = 1:1
    [uCal, pos] = max(FRCSweep.amp);
    Omega(iAmp,1) = FRCSweep.omega(pos);
    yCal(:,iAmp) = uCal*V(:,indV)./V(outdof,indV);
end

[fRed, OmegaFm, OmegaFp, uF, psiF] = miniComputeFRC(yCal, Omega, BBC, Hmap, iHmap, Tmap, iTmap, ampF);

colors = [0.1,0.7,0.1];
customFigure();
for iAmp = 1:1
    rhoMax = fsolve(@(rho) abs(OmegaFm(rho,fRed(iAmp))-BBC.freq(rho)),0.01);
    rho = linspace(0,rhoMax,100);
    realinds = imag([OmegaFm(rho,fRed(iAmp)), OmegaFp(rho,fRed(iAmp))]) < 1e-3;
    pr = @(x) x(:,realinds);
    plot(pr([OmegaFm(rho,fRed(iAmp)), OmegaFp(rho,fRed(iAmp))]), ...
        pr([uF(rho),uF(rho)]), 'Color', colors(iAmp,:), 'Linewidth', 2, ...
        'DisplayName', ['\texttt{fastSSM}$^+$']);
    plot(FRCSweep.omega, FRCSweep.amp, '.', 'MarkerSize', 16, ...
        'Color', [0,0,0], 'DisplayName', ['Simulation'],...
        'MarkerSize', 14)
end
plot(BBC.frequency, BBC.amplitude, 'k', 'DisplayName', '\texttt{fastSSM}$^+$ backbone curve')
xlim(omegaSpan)
legend('Interpreter', 'latex', 'location', 'best')
xlabel('$\Omega$ [rad/s]', 'Interpreter', 'latex')
ylabel('$u$ [m]', 'Interpreter', 'latex')
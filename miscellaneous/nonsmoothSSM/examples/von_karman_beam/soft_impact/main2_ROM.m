% Von Karman beam with soft impact: data-driven SSM-based model reduction.

clearvars; figure; classicColors = colororder; close all; format shortg; clc

load dataVKDecay2DFriction

%% Fit and validate SSMs
delta = DataInfo.delta;
%%%%%% Negative case
xDataCurr = xDataMinus;
n = size(M,1);
% Linear parts
[VMinusComplex,RMinusComplex] = eigSorted(AMinus);
WMinusComplex = inv(VMinusComplex);
VMinus = [];
for ii = 1:n
    VMinus =[VMinus real(VMinusComplex(:,ii)) imag(VMinusComplex(:,ii))];
end
WMinus = inv(VMinus);
V0 = VMinus(:,[1:2]); W0 = WMinus([1:2],:);
R0 = W0*AMinus*V0;
W0Minus = W0; V0Minus = V0;
% Chart map
chartMap = @(x) W0*x;
% Shift
shiftMinus = @(x) x-fixedpointMinus;
% Set reduced coordinates
xDataCurrShifted = funToCell(xDataCurr,shiftMinus);
etaData = funToCell(xDataCurrShifted,chartMap);
indTrain = 2; indTest = 1; SSMDim = 2; SSMOrder = 5; ROMOrder = SSMOrder;

% Parametrization
plotCoordSurf = [outdof outdof+n outdof+2];
IMInfoSMinus = IMGeometry(xDataCurrShifted(indTrain,:), SSMDim, SSMOrder, 'reducedCoordinates',etaData(indTrain,:),'Ve',V0);
IMInfoSMinus.chart.map = chartMap;
disp('Results of geometry fit:')
xRec = liftTrajectories(IMInfoSMinus, etaData);
normedTrajDist = computeTrajectoryErrors(xRec, xDataCurrShifted);
meanErrorGeoTrain = mean(normedTrajDist(indTrain))*100
meanErrorGeoTest = mean(normedTrajDist(indTest))*100
plotSSMWithTrajectories(IMInfoSMinus, plotCoordSurf, xDataCurrShifted,'Colors', [0.5; 0.1]*[1 1 1],'ColorSurf', classicColors(1,:));
xlabel('$q_{mid}^-$ [m]','interpreter','latex');
zlabel('$q_{mid,ax}^-$ [m]','interpreter','latex');
ylabel('$\dot{q}_{mid}$ [m/s]','interpreter','latex');
set(gca,'FontSize',20);
view(-10,11); legend('off')
plotSSMandTrajectories(IMInfoSMinus, @(x) WMinus(6,:)*x, xDataCurrShifted, etaData,'ColorTraj', [0.5; 0.1]*[1 1 1],'ColorSSM', classicColors(1,:));

% from shifted to unshifted
IMInfoMinus = IMInfoSMinus;
paraMap = IMInfoSMinus.parametrization.map;
IMInfoMinus.parametrization.map = @(eta) fixedpointMinus + paraMap(eta);
IMInfoMinus.chart.map = @(x) chartMap(shiftMinus(x));
plotSSMWithTrajectories(IMInfoMinus, plotCoordSurf, xDataCurr,'Colors', [0.5; 0.1]*[1 1 1],'ColorSurf', classicColors(2,:));
xlabel('$q_{mid}$ [m]','interpreter','latex');
zlabel('$q_{mid,ax}$ [m]','interpreter','latex');
ylabel('$\dot{q}_{mid}$ [m/s]','interpreter','latex');
set(gca,'FontSize',20);
view(21,18); legend('off')

% Reduced dynamics
RDInfoMinus = IMDynamicsFlow(etaData(indTrain,:), 'R_PolyOrd', ROMOrder, 'style', 'default','R_coeff',R0);
disp('Errors for the dynamics on reduced-coordinates:')
etaRec = advectRD(RDInfoMinus, etaData);
normedTrajDist = computeTrajectoryErrors(etaRec, etaData);
meanErrorDynTrain = mean(normedTrajDist(indTrain))*100
meanErrorDynTest = mean(normedTrajDist(indTest))*100
plotTrajectories(etaData(indTest,:), etaRec(indTest,:), 'm','PlotCoordinate', 1, 'DisplayName', {'Test set', 'Prediction'})
ylabel('$y_1^-$','interpreter','latex')
xlabel('$t$ [s]','interpreter','latex')
set(gca,'FontSize',20);
xlim(xDataCurrShifted{indTest,1}([1 end]))
disp('Errors for the reduced-order model (geometry & dynamics):')
xRec = liftTrajectories(IMInfoSMinus, etaRec);
normedTrajDist = computeTrajectoryErrors(xRec, xDataCurrShifted);
meanErrorROMTrain = mean(normedTrajDist(indTrain))*100
meanErrorROMTest = mean(normedTrajDist(indTest))*100
plotTrajectories(xRec(indTest,:), xDataCurrShifted(indTest,:), 'm','PlotCoordinate', outdof, 'DisplayName', {'Test set', 'Prediction'})
ylabel('$q_{mid}^-$','interpreter','latex')
xlabel('$t$ [s]','interpreter','latex')
set(gca,'FontSize',20);
xlim(xDataCurrShifted{indTest,1}([1 end]))

% Full model error
[xRec, etaRec, ~] = advect(IMInfoMinus, RDInfoMinus, xDataCurr);
normedTrajDist = computeTrajectoryErrors(xRec, xDataCurr);
meanErrorROMTrain = mean(normedTrajDist(indTrain))*100
meanErrorROMTest = mean(normedTrajDist(indTest))*100
plotTrajectories(xRec(indTest,:), xDataCurr(indTest,:), 'm','PlotCoordinate', outdof, 'DisplayName', {'Test set', 'Prediction'})
ylabel('$q_{mid}$ [m]','interpreter','latex');
xlabel('$t$ [s]','interpreter','latex')
set(gca,'FontSize',20);
xlim(xDataCurrShifted{indTest,1}([1 end]))

%%%%%% Positive case
xDataCurr = xDataPlus; tolZero = 5e-6;
% Linear parts
[VPlusComplex,RPlusComplex] = eigSorted(APlus);
WPlusComplex = inv(VPlusComplex);
VPlus = [];
for ii = 1:n
    VPlus =[VPlus real(VPlusComplex(:,ii)) imag(VPlusComplex(:,ii))];
end
WPlus = inv(VPlus);
V0 = VPlus(:,[1:2]); W0 = WPlus([1:2],:);
R0 = W0*APlus*V0;
W0Plus = W0; V0Plus = V0;
% Chart map
chartMap = @(x) W0*x;
% Shift
shiftPlus = @(x) x-fixedpointPlus;
% Set reduced coordinates
xDataCurrShifted = funToCell(xDataCurr,shiftPlus);
etaData = funToCell(xDataCurrShifted,chartMap);

% Parametrization
IMInfoSPlus = IMGeometry(xDataCurrShifted(indTrain,:), SSMDim, SSMOrder, 'reducedCoordinates',etaData(indTrain,:),'Ve',V0);
IMInfoSPlus.chart.map = chartMap;
disp('Results of geometry fit:')
xRec = liftTrajectories(IMInfoSPlus, etaData);
normedTrajDist = computeTrajectoryErrors(xRec, xDataCurrShifted);
meanErrorGeoTrain = mean(normedTrajDist(indTrain))*100;
meanErrorGeoTest = mean(normedTrajDist(indTest))*100;
plotSSMWithTrajectories(IMInfoSPlus, plotCoordSurf, xDataCurrShifted,'Colors', [0.5; 0.1]*[1 1 1],'ColorSurf', classicColors(2,:));
xlabel('$q_{mid}^+$ [m]','interpreter','latex');
zlabel('$q_{mid,ax}^+$ [m]','interpreter','latex');
ylabel('$\dot{q}_{mid}$ [m/s]','interpreter','latex');
set(gca,'FontSize',20);
view(-10,11); legend('off')

% From shifted to unshifted
IMInfoPlus = IMInfoSPlus;
paraMap = IMInfoSPlus.parametrization.map;
IMInfoPlus.parametrization.map = @(eta) fixedpointPlus + paraMap(eta);
IMInfoPlus.chart.map = @(x) chartMap(shiftPlus(x));
plotSSMWithTrajectories(IMInfoPlus, plotCoordSurf, xDataCurr,'Colors', [0.5; 0.1]*[1 1 1],'ColorSurf', classicColors(1,:));
xlabel('$q_{mid}^+$ [m]','interpreter','latex');
zlabel('$q_{mid,ax}^+$ [m]','interpreter','latex');
ylabel('$\dot{q}_{mid}$ [m/s]','interpreter','latex');
set(gca,'FontSize',20);
view(21,18);

% Reduced dynamics
RDInfoPlus = IMDynamicsFlow(etaData(indTrain,:), 'R_PolyOrd', ROMOrder, 'style', 'default','R_coeff',R0);
disp('Errors for the dynamics on reduced-coordinates:')
etaRec = advectRD(RDInfoPlus, etaData);
normedTrajDist = computeTrajectoryErrors(etaRec, etaData);
meanErrorDynTrain = mean(normedTrajDist(indTrain))*100
meanErrorDynTest = mean(normedTrajDist(indTest))*100
plotTrajectories(etaData(indTest,:), etaRec(indTest,:), 'm','PlotCoordinate', 1, 'DisplayName', {'Test set', 'Prediction'})
ylabel('$y_1^+$','interpreter','latex')
xlabel('$t$ [s]','interpreter','latex')
set(gca,'FontSize',20);
xlim(xDataCurrShifted{indTest,1}([1 end]))
disp('Errors for the reduced-order model (geometry & dynamics):')
xRec = liftTrajectories(IMInfoSPlus, etaRec);
normedTrajDist = computeTrajectoryErrors(xRec, xDataCurrShifted);
meanErrorROMTrain = mean(normedTrajDist(indTrain))*100
meanErrorROMTest = mean(normedTrajDist(indTest))*100
plotTrajectories(xRec(indTest,:), xDataCurrShifted(indTest,:), 'm','PlotCoordinate', outdof, 'DisplayName', {'Test set', 'Prediction'})
ylabel('$q_{mid}^+$ [m]','interpreter','latex')
xlabel('$t$ [s]','interpreter','latex')
set(gca,'FontSize',20);
xlim(xDataCurrShifted{indTest,1}([1 end]))

% Full model error
[xRec, etaRec, ~] = advect(IMInfoPlus, RDInfoPlus, xDataCurr);
normedTrajDist = computeTrajectoryErrors(xRec, xDataCurr);
meanErrorROMTrain = mean(normedTrajDist(indTrain))*100
meanErrorROMTest = mean(normedTrajDist(indTest))*100
plotTrajectories(xRec(indTest,:), xDataCurr(indTest,:), 'm','PlotCoordinate', outdof, 'DisplayName', {'Test set', 'Prediction'})
ylabel('$q_{mid}$ [m]','interpreter','latex');
xlabel('$t$ [s]','interpreter','latex')
set(gca,'FontSize',20);
xlim(xDataCurrShifted{indTest,1}([1 end]))

%% Get SSM intersections

sFunctionPhys = @(x) x(outdof);
% Definitions
wPlus = @(t,x) IMInfoPlus.chart.map(x);
vPlus = @(t,eta) IMInfoPlus.parametrization.map(eta);
rPlus = @(t,eta) RDInfoPlus.reducedDynamics.map(eta);
wMinus = @(t,x) IMInfoMinus.chart.map(x);
vMinus = @(t,eta) IMInfoMinus.parametrization.map(eta);
rMinus = @(t,eta) RDInfoMinus.reducedDynamics.map(eta);
wFunction = @(t,x,l) 0.5*(1+l)*wPlus(t,x)+0.5*(1-l)*wMinus(t,x);
rFunction = @(t,eta,l) 0.5*(1+l)*rPlus(t,eta)+0.5*(1-l)*rMinus(t,eta);
vFunction = @(t,eta,l) 0.5*(1+l)*vPlus(t,eta)+0.5*(1-l)*vMinus(t,eta);

sToEta = @(eta,l) [0 1]*eta;
eta2_Max = 35; 

l0 = -1; 
eta2_vec = linspace(0,eta2_Max,101);
eta1_vecPlus = zeros(1,length(eta2_vec));
eta1_0 = 0;
for ii = 1:length(eta2_vec)
    eta1_0 = fsolve(@(eta1) ...
               sFunctionPhys(vFunction(0,[eta1; eta2_vec(ii)],l0)),eta1_0,options); 
    eta1_vecPlus(ii) = eta1_0;
end
eta1_vecMinus = zeros(1,length(eta2_vec));
eta1_0 = 0;
for ii = 2:length(eta2_vec)
    eta1_0 = fsolve(@(eta1) ...
               sFunctionPhys(vFunction(0,[eta1; eta2_vec(ii)],l0)),eta1_0,options);  
    eta1_vecMinus(ii) = eta1_0;
end
eta2_vec= [-fliplr(eta2_vec(2:end)) eta2_vec(1:end)];
eta1_vec= [fliplr(eta1_vecMinus(2:end)) eta1_vecPlus(1:end)];
gammaMinus = @(s) [s; interp1(eta2_vec, eta1_vec, s)];

l0 = 1; 
eta2_vec = linspace(0,eta2_Max,1001);
eta1_vecPlus = zeros(1,length(eta2_vec));
eta1_0 = 0;
for ii = 1:length(eta2_vec)
    eta1_0 = fsolve(@(eta1) ...
               sFunctionPhys(vFunction(0,[eta1; -eta2_vec(ii)],l0)),eta1_0,options); 
    eta1_vecPlus(ii) = eta1_0;
end
eta1_vecMinus = zeros(1,length(eta2_vec));
eta1_0 = 0;
for ii = 2:length(eta2_vec)
    eta1_0 = fsolve(@(eta1) ...
               sFunctionPhys(vFunction(0,[eta1; -eta2_vec(ii)],l0)),eta1_0,options); 
    eta1_vecMinus(ii) = eta1_0;
end
eta2_vec= [-fliplr(eta2_vec(2:end)) eta2_vec(1:end)];
eta1_vec= [fliplr(eta1_vecMinus(2:end)) eta1_vecPlus(1:end)];
gammaPlus = @(s) [s; interp1(eta2_vec, eta1_vec, s)];

gammaFunction = @(s,l) 0.5*(1+l)*gammaPlus(s)+0.5*(1-l)*gammaMinus(s); 
WeightMatrix= eye(2*n); 

%% Integrate non-smooth model

idxIni = 1;
timeSpan = [xDataNonSmooth{1,1}(idxIni) xDataNonSmooth{1,1}(xDataNonSmooth{1,3}(end))];

% Integration
[t,xRec,indSwitch,~] = integratePWSROM(wFunction,rFunction,vFunction,sFunctionPhys,timeSpan, xDataNonSmooth{1,2}(:,idxIni),...
                       'RICmethod','fiber','WeightMatrix', WeightMatrix,'gammaStoEta',gammaFunction,'gammaEtatoS',sToEta);
xRec = transpose(xRec);

slowTimeScale = 2*pi/abs(lambda(1));

colorFOM = [0 0 0];
colorROM = classicColors(5,:);
customFigure('subPlot',[2 1]);tLimits = timeSpan([1 end]);
subplot(211); xlim(tLimits)
plot(xDataNonSmooth{1,1},xDataNonSmooth{1,2}(outdof,:),'Linewidth',2,'Color', colorFOM);
plot(t,xRec(outdof,:),'Linewidth',2,'Color', colorROM);
plot(xDataNonSmooth{1,1}(xDataNonSmooth{1,3}),xDataNonSmooth{1,2}(outdof,xDataNonSmooth{1,3}),'.','Color', colorFOM,'MarkerSize',14);
plot(t(indSwitch),xRec(outdof,indSwitch),'.','Color', colorROM,'MarkerSize',10);
legend('FOM','ROM','S-FOM','S-ROM','numColumns',2)
xlabel('$t$ [s]','interpreter','latex')
ylabel('$q_{mid}$ [m]','interpreter','latex')
set(gca,'FontSize',20)
subplot(212); xlim(tLimits)
plot(xDataNonSmooth{1,1},xDataNonSmooth{1,2}(outdof+n,:),'Linewidth',2,'Color', colorFOM);
plot(t,xRec(outdof+n,:),'Linewidth',2,'Color', colorROM);
plot(xDataNonSmooth{1,1}(xDataNonSmooth{1,3}),xDataNonSmooth{1,2}(outdof+n,xDataNonSmooth{1,3}),'.','Color', colorFOM,'MarkerSize',14);
plot(t(indSwitch),xRec(outdof+n,indSwitch),'.','Color', colorROM,'MarkerSize',10);
xlabel('$t$ [s]','interpreter','latex')
ylabel('$\dot{q}_{mid}$ [m/s]','interpreter','latex');
set(gca,'FontSize',20)

%%%%%%%%%
customFigure('subPlot',[2 2]); tLimitsZ = [0 2*slowTimeScale]+xDataNonSmooth{1,1}(idxIni);
subplot(221); xlim(tLimitsZ)
plot(xDataNonSmooth{1,1},xDataNonSmooth{1,2}(outdof,:),'Linewidth',2,'Color', colorFOM);
plot(t,xRec(outdof,:),'--','Linewidth',2,'Color', colorROM);
plot(xDataNonSmooth{1,1}(xDataNonSmooth{1,3}),xDataNonSmooth{1,2}(outdof,xDataNonSmooth{1,3}),'.','Color', colorFOM,'MarkerSize',14);
plot(t(indSwitch),xRec(outdof,indSwitch),'.','Color', colorROM,'MarkerSize',10);
xlabel('$t$ [s]','interpreter','latex')
ylabel('$q_{mid}$ [m]','interpreter','latex')
subplot(223); xlim(tLimitsZ)
plot(xDataNonSmooth{1,1},xDataNonSmooth{1,2}(outdof+n,:),'Linewidth',2,'Color', colorFOM);
plot(t,xRec(outdof+n,:),'--','Linewidth',2,'Color', colorROM);
plot(xDataNonSmooth{1,1}(xDataNonSmooth{1,3}),xDataNonSmooth{1,2}(outdof+n,xDataNonSmooth{1,3}),'.','Color', colorFOM,'MarkerSize',14);
plot(t(indSwitch),xRec(outdof+n,indSwitch),'.','Color', colorROM,'MarkerSize',10);
xlabel('$t$ [s]','interpreter','latex')
ylabel('$\dot{q}_{mid}$ [m/s]','interpreter','latex');
tLimitsZ = [-2*slowTimeScale 0]+xDataNonSmooth{1,1}(xDataNonSmooth{1,3}(end));
subplot(222); xlim(tLimitsZ)
plot(xDataNonSmooth{1,1},xDataNonSmooth{1,2}(outdof,:),'Linewidth',2,'Color', colorFOM);
plot(t,xRec(outdof,:),'--','Linewidth',2,'Color', colorROM);
plot(xDataNonSmooth{1,1}(xDataNonSmooth{1,3}),xDataNonSmooth{1,2}(outdof,xDataNonSmooth{1,3}),'.','Color', colorFOM,'MarkerSize',14);
plot(t(indSwitch),xRec(outdof,indSwitch),'.','Color', colorROM,'MarkerSize',10);
xlabel('$t$ [s]','interpreter','latex')
ylabel('$q_{mid}$ [m]','interpreter','latex')
subplot(224); xlim(tLimitsZ)
plot(xDataNonSmooth{1,1},xDataNonSmooth{1,2}(outdof+n,:),'Linewidth',2,'Color', colorFOM);
plot(t,xRec(outdof+n,:),'--','Linewidth',2,'Color', colorROM);
plot(xDataNonSmooth{1,1}(xDataNonSmooth{1,3}),xDataNonSmooth{1,2}(outdof+n,xDataNonSmooth{1,3}),'.','Color', colorFOM,'MarkerSize',14);
plot(t(indSwitch),xRec(outdof+n,indSwitch),'.','Color', colorROM,'MarkerSize',10);
xlabel('$t$ [s]','interpreter','latex')
ylabel('$\dot{q}_{mid}$ [m/s]','interpreter','latex');

%% Non-autonomous system

load frc_data

% Define forced models
[IMInfoFMinus,RDInfoFMinus] = forcedSSMROMgraph(IMInfoMinus,RDInfoMinus,...
    forcingVectors_na,WMinusComplex([SSMDim+1:end],:),...
    diag(RMinusComplex([SSMDim+1:end],[SSMDim+1:end])),...
    VMinusComplex(:,[SSMDim+1:end]),W0Minus,V0Minus);

[IMInfoFPlus,RDInfoFPlus] = forcedSSMROMgraph(IMInfoPlus,RDInfoPlus,...
    forcingVectors_na,WPlusComplex([SSMDim+1:end],:),...
    diag(RPlusComplex([SSMDim+1:end],[SSMDim+1:end])),...
    VPlusComplex(:,[SSMDim+1:end]),W0Plus,V0Plus);

amplitudeFunction = @(x) x(outdof,:);

ii = 0;
FRC_reduced = cell(2,length(forcingAmplitudes));

for forAmp = forcingAmplitudes
    ii = ii + 1;
   
    forFreq = freq_0;
     % Definitions
    wPlus = @(t,x) IMInfoFPlus.chart.map(x);
    vPlus = @(t,eta) IMInfoFPlus.parametrization.map(t,eta,forFreq,forAmp,0);
    rPlus = @(t,eta) RDInfoFPlus.reducedDynamics.map(t,eta,forFreq,forAmp,0);
    wMinus = @(t,x) IMInfoFMinus.chart.map(x);
    vMinus = @(t,eta) IMInfoFMinus.parametrization.map(t,eta,forFreq,forAmp,0);
    rMinus = @(t,eta) RDInfoFMinus.reducedDynamics.map(t,eta,forFreq,forAmp,0);
    wFunction = @(t,x,l) 0.5*(1+l)*wPlus(t,x)+0.5*(1-l)*wMinus(t,x);
    rFunction = @(t,eta,l) 0.5*(1+l)*rPlus(t,eta)+0.5*(1-l)*rMinus(t,eta);
    vFunction = @(t,eta,l) 0.5*(1+l)*vPlus(t,eta)+0.5*(1-l)*vMinus(t,eta);


    % Initial condition
    ampsLin = (K-eye(n)*forFreq^2)\forcingVectors_na(n+1:end)*forAmp;
    x0 = [ampsLin; zeros(n,1)];
    timeSpan = [0 nPers*2*pi/forFreq];
    % Integration
    [t,xRec,indSwitch,~] = integratePWSROM(wFunction,rFunction,vFunction,sFunctionPhys,timeSpan, x0,...
            'RICmethod','not_physical','WeightMatrix', WeightMatrix,'gammaStoEta',gammaFunction,'gammaEtatoS',sToEta);
    xRec = transpose(xRec);
    idxLastLoop = sum(t<(nPers-1)*2*pi/forFreq);
    amp = max(abs(amplitudeFunction(xRec(:,idxLastLoop:end))));
    IC = xRec(:,end);

    nstepsLow = 101;
    freqsLow = linspace(detunings(1),forFreq,nstepsLow);
    ICsLow = zeros(2*n,nstepsLow);
    ampsLow = zeros(1,nstepsLow);
    ICsLow(:,end) = IC; ampsLow(end) = amp;
    nstepsRef = 101;
    freqsHigh = linspace(forFreq,detunings(2),nstepsLow);
    nstepsHigh = length(freqsHigh);
    ICsHigh  = zeros(2*n,nstepsHigh);
    ampsHigh  = zeros(1,nstepsHigh);
    ICsHigh(:,1) = IC; ampsHigh(1) = amp;

    for iFreq = nstepsLow-1:-1:1
            forFreq = freqsLow(iFreq); x0 = ICsLow(:,iFreq+1);
            timeSpan = [0 nPers*2*pi/forFreq];
            % Definitions
            wPlus = @(t,x) IMInfoFPlus.chart.map(x);
            vPlus = @(t,eta) IMInfoFPlus.parametrization.map(t,eta,forFreq,forAmp,0);
            rPlus = @(t,eta) RDInfoFPlus.reducedDynamics.map(t,eta,forFreq,forAmp,0);
            wMinus = @(t,x) IMInfoFMinus.chart.map(x);
            vMinus = @(t,eta) IMInfoFMinus.parametrization.map(t,eta,forFreq,forAmp,0);
            rMinus = @(t,eta) RDInfoFMinus.reducedDynamics.map(t,eta,forFreq,forAmp,0);
            wFunction = @(t,x,l) 0.5*(1+l)*wPlus(t,x)+0.5*(1-l)*wMinus(t,x);
            rFunction = @(t,eta,l) 0.5*(1+l)*rPlus(t,eta)+0.5*(1-l)*rMinus(t,eta);
            vFunction = @(t,eta,l) 0.5*(1+l)*vPlus(t,eta)+0.5*(1-l)*vMinus(t,eta);
            % Integration
            [t,xRec,indSwitch,~] = integratePWSROM(wFunction,rFunction,vFunction,sFunctionPhys,timeSpan, x0,...
                'RICmethod','non_physical','WeightMatrix', WeightMatrix,'gammaStoEta',gammaFunction,'gammaEtatoS',sToEta);
            xRec = transpose(xRec);
            idxLastLoop = sum(t<(nPers-1)*2*pi/forFreq);
            ampsLow(iFreq) = max(abs(amplitudeFunction(xRec(:,idxLastLoop:end))));
            ICsLow(:,iFreq) = xRec(:,end);            
                      
    end
    for iFreq = 2:nstepsHigh
            forFreq = freqsHigh(iFreq); x0 = ICsHigh(:,iFreq-1);
            timeSpan = [0 nPers*2*pi/forFreq];
            % Definitions
            wPlus = @(t,x) IMInfoFPlus.chart.map(x);
            vPlus = @(t,eta) IMInfoFPlus.parametrization.map(t,eta,forFreq,forAmp,0);
            rPlus = @(t,eta) RDInfoFPlus.reducedDynamics.map(t,eta,forFreq,forAmp,0);
            wMinus = @(t,x) IMInfoFMinus.chart.map(x);
            vMinus = @(t,eta) IMInfoFMinus.parametrization.map(t,eta,forFreq,forAmp,0);
            rMinus = @(t,eta) RDInfoFMinus.reducedDynamics.map(t,eta,forFreq,forAmp,0);
            wFunction = @(t,x,l) 0.5*(1+l)*wPlus(t,x)+0.5*(1-l)*wMinus(t,x);
            rFunction = @(t,eta,l) 0.5*(1+l)*rPlus(t,eta)+0.5*(1-l)*rMinus(t,eta);
            vFunction = @(t,eta,l) 0.5*(1+l)*vPlus(t,eta)+0.5*(1-l)*vMinus(t,eta);
            % Integration
            [t,xRec,indSwitch,~] = integratePWSROM(wFunction,rFunction,vFunction,sFunctionPhys,timeSpan, x0);
            xRec = transpose(xRec);
            idxLastLoop = sum(t<(nPers-1)*2*pi/forFreq);
            ampsHigh(iFreq) = max(abs(amplitudeFunction(xRec(:,idxLastLoop:end))));
            ICsHigh(:,iFreq) = xRec(:,end);

    end
    FRC_reduced{1,ii} = [ampsLow ampsHigh(2:end)];
    FRC_reduced{2,ii} = [freqsLow freqsHigh(2:end)];

end

save('FRC_reduced','FRC_reduced');

FRC_reduced_reverse = cell(2,length(forcingAmplitudes));

detunings = [freq_0 * 0.97 freq_0 * 1.05];

ii = 0;
for forAmp = forcingAmplitudes
    ii = ii + 1;
    forFreq = detunings(2);
     % Definitions
    wPlus = @(t,x) IMInfoFPlus.chart.map(x);
    vPlus = @(t,eta) IMInfoFPlus.parametrization.map(t,eta,forFreq,forAmp,0);
    rPlus = @(t,eta) RDInfoFPlus.reducedDynamics.map(t,eta,forFreq,forAmp,0);
    wMinus = @(t,x) IMInfoFMinus.chart.map(x);
    vMinus = @(t,eta) IMInfoFMinus.parametrization.map(t,eta,forFreq,forAmp,0);
    rMinus = @(t,eta) RDInfoFMinus.reducedDynamics.map(t,eta,forFreq,forAmp,0);
    wFunction = @(t,x,l) 0.5*(1+l)*wPlus(t,x)+0.5*(1-l)*wMinus(t,x);
    rFunction = @(t,eta,l) 0.5*(1+l)*rPlus(t,eta)+0.5*(1-l)*rMinus(t,eta);
    vFunction = @(t,eta,l) 0.5*(1+l)*vPlus(t,eta)+0.5*(1-l)*vMinus(t,eta);


    % Initial condition
    ampsLin = (K-eye(n)*forFreq^2)\forcingVectors_na(n+1:end)*forAmp;
    x0 = [ampsLin; zeros(n,1)];
    timeSpan = [0 nPers*2*pi/forFreq];
    % Integration
    [t,xRec,indSwitch,~] = integratePWSROM(wFunction,rFunction,vFunction,sFunctionPhys,timeSpan, x0,...
            'RICmethod','not_physical','WeightMatrix', WeightMatrix,'gammaStoEta',gammaFunction,'gammaEtatoS',sToEta);
    xRec = transpose(xRec);
    idxLastLoop = sum(t<(nPers-1)*2*pi/forFreq);
    amp = max(abs(amplitudeFunction(xRec(:,idxLastLoop:end))));
    IC = xRec(:,end);

    nstepsLow = 51;
    freqsLow = linspace(detunings(1),forFreq,nstepsLow);
    ICsLow = zeros(2*n,nstepsLow);
    ampsLow = zeros(1,nstepsLow);
    ICsLow(:,end) = IC; ampsLow(end) = amp;

    for iFreq = nstepsLow-1:-1:1
            forFreq = freqsLow(iFreq); x0 = ICsLow(:,iFreq+1);
            timeSpan = [0 nPers*2*pi/forFreq];
            % Definitions
            wPlus = @(t,x) IMInfoFPlus.chart.map(x);
            vPlus = @(t,eta) IMInfoFPlus.parametrization.map(t,eta,forFreq,forAmp,0);
            rPlus = @(t,eta) RDInfoFPlus.reducedDynamics.map(t,eta,forFreq,forAmp,0);
            wMinus = @(t,x) IMInfoFMinus.chart.map(x);
            vMinus = @(t,eta) IMInfoFMinus.parametrization.map(t,eta,forFreq,forAmp,0);
            rMinus = @(t,eta) RDInfoFMinus.reducedDynamics.map(t,eta,forFreq,forAmp,0);
            wFunction = @(t,x,l) 0.5*(1+l)*wPlus(t,x)+0.5*(1-l)*wMinus(t,x);
            rFunction = @(t,eta,l) 0.5*(1+l)*rPlus(t,eta)+0.5*(1-l)*rMinus(t,eta);
            vFunction = @(t,eta,l) 0.5*(1+l)*vPlus(t,eta)+0.5*(1-l)*vMinus(t,eta);
            % Integration
            [t,xRec,indSwitch,~] = integratePWSROM(wFunction,rFunction,vFunction,sFunctionPhys,timeSpan, x0,...
                'RICmethod','non_physical','WeightMatrix', WeightMatrix,'gammaStoEta',gammaFunction,'gammaEtatoS',sToEta);
            xRec = transpose(xRec);
            idxLastLoop = sum(t<(nPers-1)*2*pi/forFreq);
            ampsLow(iFreq) = max(abs(amplitudeFunction(xRec(:,idxLastLoop:end))));
            ICsLow(:,iFreq) = xRec(:,end);            
                    
    end
    
    FRC_reduced_reverse{1,ii} = [ampsLow];
    FRC_reduced_reverse{2,ii} = [freqsLow];

end

save('FRC_reduced_reverse','FRC_reduced_reverse');

%% plot the curves separately
colors = colororder;

load frc_data 
load FRC_reduced
load frc_data_reverse
load FRC_reduced_reverse

% amplitude 25
customFigure();
plot(FRC_full{2,1},FRC_full{1,1},'k','Marker','o','linestyle','none','MarkerSize',12,'linewidth',3,'DisplayName','Full-order model');
plot(FRC_full_reverse{2,1},FRC_full_reverse{1,1},'k','Marker','o','linestyle','none','MarkerSize',12,'linewidth',3,'HandleVisibility','off');
plot(FRC_reduced{2,1},FRC_reduced{1,1},'.g','MarkerSize',12,'DisplayName','Reduced-order model');
plot(FRC_reduced_reverse{2,1},FRC_reduced_reverse{1,1},'.g','MarkerSize',12,'HandleVisibility','off');
xlabel('Forcing frequency [rad/s]','Interpreter','latex')
ylabel('Amplitude [m]','Interpreter','latex')
ax = gca; 
ax.FontSize = 20;
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';
ylim([0 8e-3])
legend('location','NE','Interpreter','latex')
ylim([0 8e-3])

% amplitude 35
customFigure();
plot(FRC_full{2,2},FRC_full{1,2},'k','Marker','o','linestyle','none','MarkerSize',12,'linewidth',3,'DisplayName','Full-order model');
plot(FRC_full_reverse{2,2},FRC_full_reverse{1,2},'k','Marker','o','linestyle','none','MarkerSize',12,'linewidth',3,'HandleVisibility','off');
plot(FRC_reduced{2,2},FRC_reduced{1,2},'.g','MarkerSize',12,'DisplayName','Reduced-order model');
plot(FRC_reduced_reverse{2,2},FRC_reduced_reverse{1,2},'.g','MarkerSize',12,'HandleVisibility','off');
xlabel('Forcing frequency [rad/s]','Interpreter','latex')
ylabel('Amplitude [m]','Interpreter','latex')
ax = gca; 
ax.FontSize = 20;
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';
ylim([0 8e-3])
legend('location','NE','Interpreter','latex')
ylim([0 8e-3])

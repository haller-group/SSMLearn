% Von Karman beam with with friction on a moving belt: data-driven SSM-based model reduction.

clearvars; figure; classicColors = colororder; close all; format shortg; clc

load dataVKDecay2DFriction
%% Fit and validate SSMs
delta = DataInfo.delta;

%%%%%% Negative case
xDataCurr = xDataMinus; tolZero = 2e-6;
n = size(M,1);

% Chart map
W0 = zeros(2,2*n); W0(1,outdof) = 1; W0(2,outdof+n) = 1;
chartMap = @(x) W0*x;
% Shift
shiftMinus = @(x) x-fixedpointMinus;

% Set reduced coordinates
xDataCurrShifted = funToCell(xDataCurr,shiftMinus);
etaData = funToCell(xDataCurrShifted,chartMap);
indTrain = 2; indTest = 1; SSMDim = 2; SSMOrder = 9; ROMOrder = SSMOrder;

% Linear parts
[VeMinus,ReMinus] = eigSorted(AMinus); WeMinus = inv(VeMinus);
VMinus = VeMinus; WMinus = WeMinus; RMinus = ReMinus;
VeMinus = VeMinus(:,[1 n+1]); WeMinus = WeMinus([1 n+1],:); ReMinus = ReMinus([1 n+1], [1 n+1]);
PassMat = W0*VeMinus;
V0 = real(VeMinus/PassMat); V0 = V0.*(abs(V0)>tolZero);
R0 = real((PassMat*ReMinus)/PassMat); R0 = R0.*(abs(R0)>tolZero);

V0Minus = V0;
W0Minus = W0;
R0Minus = R0;

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
set(gca,'FontSize',20)
view(-10,11); legend('off')

% from shifted to unshifted
IMInfoMinus = IMInfoSMinus;
paraMap = IMInfoSMinus.parametrization.map;
IMInfoMinus.parametrization.map = @(eta) fixedpointMinus + paraMap(eta);
IMInfoMinus.chart.map = @(x) chartMap(shiftMinus(x));
plotSSMWithTrajectories(IMInfoMinus, plotCoordSurf, xDataCurr,'Colors', [0.5; 0.1]*[1 1 1],'ColorSurf', classicColors(1,:));
xlabel('$q_{mid}$ [m]','interpreter','latex');
zlabel('$q_{mid,ax}$ [m]','interpreter','latex');
ylabel('$\dot{q}_{mid}$ [m/s]','interpreter','latex');
set(gca,'FontSize',20)
view(21,18); legend('off')

% Reduced dynamics

RDInfoMinus = IMDynamicsMech(etaData(indTrain,:), 'R_PolyOrd', ROMOrder, 'style', 'default','R_coeff',R0);
disp('Errors for the dynamics on reduced-coordinates:')
etaRec = advectRD(RDInfoMinus, etaData);
normedTrajDist = computeTrajectoryErrors(etaRec, etaData);
meanErrorDynTrain = mean(normedTrajDist(indTrain))*100
meanErrorDynTest = mean(normedTrajDist(indTest))*100
plotTrajectories(etaData(indTest,:), etaRec(indTest,:), 'm','PlotCoordinate', 1, 'DisplayName', {'Test set', 'Prediction'})
ylabel('$q_1$ [m]','interpreter','latex')
disp('Errors for the reduced-order model (geometry & dynamics):')
xRec = liftTrajectories(IMInfoSMinus, etaRec);
normedTrajDist = computeTrajectoryErrors(xRec, xDataCurrShifted);
meanErrorROMTrain = mean(normedTrajDist(indTrain))*100
meanErrorROMTest = mean(normedTrajDist(indTest))*100
plotTrajectories(xRec(indTest,:), xDataCurrShifted(indTest,:), 'm','PlotCoordinate', outdof, 'DisplayName', {'Test set', 'Prediction'})
ylabel('$q_{mid}^-$','interpreter','latex')
xlabel('$t [s]$','interpreter','latex')
set(gca,'FontSize',20)
xlim(xDataCurrShifted{indTest,1}([1 end]))


% Full model error
[xRec, etaRec, ~] = advect(IMInfoMinus, RDInfoMinus, xDataCurr);
normedTrajDist = computeTrajectoryErrors(xRec, xDataCurr);
meanErrorROMTrain = mean(normedTrajDist(indTrain))*100
meanErrorROMTest = mean(normedTrajDist(indTest))*100
plotTrajectories(xRec(indTest,:), xDataCurr(indTest,:), 'm','PlotCoordinate', outdof, 'DisplayName', {'Test set', 'Prediction'})
ylabel('$q_{mid}$ [m]','interpreter','latex');
xlabel('$t$ [s]','interpreter','latex')
set(gca,'FontSize',20)
xlim(xDataCurrShifted{indTest,1}([1 end]))

% positive case 
xDataCurr = xDataPlus; tolZero = 2e-6; tLimits = [0 600];
n = size(M,1);

% Chart map
W0 = zeros(2,2*n); W0(1,outdof) = 1; W0(2,outdof+n) = 1;
chartMap = @(x) W0*x;
% Shift
shiftPlus = @(x) x-fixedpointPlus;

% Set reduced coordinates
xDataCurrShifted = funToCell(xDataCurr,shiftPlus);
etaData = funToCell(xDataCurrShifted,chartMap);
indTrain = 2; indTest = 1; SSMDim = 2; SSMOrder = 9; ROMOrder = SSMOrder;

% Linear parts
[VePlus,RePlus] = eigSorted(APlus); WePlus = inv(VePlus);
VPlus = VePlus; WPlus = WePlus; RPlus = RePlus;
VePlus = VePlus(:,[1 n+1]); WePlus = WePlus([1 n+1],:); RePlus = RePlus([1 n+1], [1 n+1]);
PassMat = W0*VePlus;
V0 = real(VePlus/PassMat); V0 = V0.*(abs(V0)>tolZero);
R0 = real((PassMat*RePlus)/PassMat); R0 = R0.*(abs(R0)>tolZero);

V0Plus = V0;
W0Plus = W0;
R0Plus = R0;

% Parametrization
plotCoordSurf = [outdof outdof+n outdof+2];
IMInfoSPlus = IMGeometry(xDataCurrShifted(indTrain,:), SSMDim, SSMOrder, 'reducedCoordinates',etaData(indTrain,:),'Ve',V0);
IMInfoSPlus.chart.map = chartMap;
disp('Results of geometry fit:')
xRec = liftTrajectories(IMInfoSPlus, etaData);
normedTrajDist = computeTrajectoryErrors(xRec, xDataCurrShifted);
meanErrorGeoTrain = mean(normedTrajDist(indTrain))*100
meanErrorGeoTest = mean(normedTrajDist(indTest))*100
plotSSMWithTrajectories(IMInfoSPlus, plotCoordSurf, xDataCurrShifted,'Colors', [0.5; 0.1]*[1 1 1],'ColorSurf', classicColors(1,:));
xlabel('$q_{mid}^-$ [m]','interpreter','latex');
zlabel('$q_{mid,ax}^-$ [m]','interpreter','latex');
ylabel('$\dot{q}_{mid}$ [m/s]','interpreter','latex');
set(gca,'FontSize',20)
view(-10,11); legend('off')

% from shifted to unshifted
IMInfoPlus = IMInfoSPlus;
paraMap = IMInfoSPlus.parametrization.map;
IMInfoPlus.parametrization.map = @(eta) fixedpointPlus + paraMap(eta);
IMInfoPlus.chart.map = @(x) chartMap(shiftPlus(x));
plotSSMWithTrajectories(IMInfoPlus, plotCoordSurf, xDataCurr,'Colors', [0.5; 0.1]*[1 1 1],'ColorSurf', classicColors(1,:));
xlabel('$q_{mid}$ [m]','interpreter','latex');
zlabel('$q_{mid,ax}$ [m]','interpreter','latex');
ylabel('$\dot{q}_{mid}$ [m/s]','interpreter','latex');
set(gca,'FontSize',20)
view(21,18); legend('off')

% Reduced dynamics

RDInfoPlus = IMDynamicsMech(etaData(indTrain,:), 'R_PolyOrd', ROMOrder, 'style', 'default','R_coeff',R0);
disp('Errors for the dynamics on reduced-coordinates:')
etaRec = advectRD(RDInfoPlus, etaData);
normedTrajDist = computeTrajectoryErrors(etaRec, etaData);
meanErrorDynTrain = mean(normedTrajDist(indTrain))*100
meanErrorDynTest = mean(normedTrajDist(indTest))*100
plotTrajectories(etaData(indTest,:), etaRec(indTest,:), 'm','PlotCoordinate', 1, 'DisplayName', {'Test set', 'Prediction'})
ylabel('$q_1$ [m]','interpreter','latex')
disp('Errors for the reduced-order model (geometry & dynamics):')
xRec = liftTrajectories(IMInfoSPlus, etaRec);
normedTrajDist = computeTrajectoryErrors(xRec, xDataCurrShifted);
meanErrorROMTrain = mean(normedTrajDist(indTrain))*100
meanErrorROMTest = mean(normedTrajDist(indTest))*100
plotTrajectories(xRec(indTest,:), xDataCurrShifted(indTest,:), 'm','PlotCoordinate', outdof, 'DisplayName', {'Test set', 'Prediction'})
ylabel('$q_{mid}^-$','interpreter','latex')
xlabel('$t$ [s]','interpreter','latex')
set(gca,'FontSize',20)
xlim(xDataCurrShifted{indTest,1}([1 end]))

% Full model error
[xRec, etaRec, ~] = advect(IMInfoPlus, RDInfoPlus, xDataCurr);
normedTrajDist = computeTrajectoryErrors(xRec, xDataCurr);
meanErrorROMTrain = mean(normedTrajDist(indTrain))*100
meanErrorROMTest = mean(normedTrajDist(indTest))*100
plotTrajectories(xRec(indTest,:), xDataCurr(indTest,:), 'm','PlotCoordinate', outdof, 'DisplayName', {'Test set', 'Prediction'})
ylabel('$q_{mid}$ [m]','interpreter','latex');
xlabel('$t$ [s]','interpreter','latex')
set(gca,'FontSize',20)
xlim(xDataCurrShifted{indTest,1}([1 end]))

%% Integrate non-smooth model

ndof = n;
for iTrajNS = 1:1
    rdof = 1; velTransp = v_ground; 
    sFunction = @(x) x(rdof+1)-velTransp; velVec = zeros(1,2*rdof); velVec(rdof+1) = 1; 
    DxsFunction = @(x) velVec;
    redEqPlus = IMInfoSPlus.chart.map(fixedpointPlus);
    redEqMinus = IMInfoSMinus.chart.map(fixedpointMinus);
    fPlus = @(t,eta) RDInfoPlus.reducedDynamics.map(eta-redEqPlus);
    fMinus = @(t,eta) RDInfoMinus.reducedDynamics.map(eta-redEqMinus);
    fFunction = @(t,x,l) 0.5*(1+l)*fPlus(t,x)+0.5*(1-l)*fMinus(t,x);
    idxIni = 1;
    [t,eta,indSwitch,l0vect] = odenonsmooth_1switchsurf(fFunction, sFunction, DxsFunction, [xDataNonSmooth{iTrajNS,1}(idxIni) xDataNonSmooth{iTrajNS,1}(end)], xDataNonSmooth{iTrajNS,2}([outdof n+outdof],idxIni));
    eta = transpose(eta);
    timeSpan = [xDataNonSmooth{iTrajNS,1}(idxIni) xDataNonSmooth{iTrajNS,1}(end)];
    tLimits = timeSpan([1 end]);
    customFigure('subPlot',[2 1]);
    subplot(211); xlim(tLimits)
    plot(xDataNonSmooth{iTrajNS,1},xDataNonSmooth{iTrajNS,2}(outdof,:),'Linewidth',1,'Color', classicColors(1,:));
    plot(t,eta(1,:),'--','Linewidth',1,'Color', classicColors(2,:));
    plot(xDataNonSmooth{iTrajNS,1}(xDataNonSmooth{iTrajNS,3}),xDataNonSmooth{iTrajNS,2}(outdof,xDataNonSmooth{iTrajNS,3}),'.','Color', classicColors(1,:),'MarkerSize',14);
    plot(t(indSwitch),eta(1,indSwitch),'.','Color', classicColors(2,:),'MarkerSize',10);
    legend('FOM','ROM','S-FOM','S-ROM')
    xlabel('$t$ [s]','interpreter','latex')
    ylabel('$q_1$ [m]','interpreter','latex')
    set(gca,'FontSize',20)
    subplot(212); xlim(tLimits)
    plot(xDataNonSmooth{iTrajNS,1},xDataNonSmooth{iTrajNS,2}(outdof+n,:),'Linewidth',1,'Color', classicColors(1,:));
    plot(t,eta(rdof+1,:),'--','Linewidth',1,'Color', classicColors(2,:));
    plot(xDataNonSmooth{iTrajNS,1}(xDataNonSmooth{iTrajNS,3}),xDataNonSmooth{iTrajNS,2}(outdof+n,xDataNonSmooth{iTrajNS,3}),'.','Color', classicColors(1,:),'MarkerSize',14);
    plot(t(indSwitch),eta(rdof+1,indSwitch),'.','Color', classicColors(2,:),'MarkerSize',10);
    xlabel('$t$ [s]','interpreter','latex')
    ylabel('$\dot{q}_1$ [m]','interpreter','latex')
    set(gca,'FontSize',20)

end

%% Comparison in phase space 

figure
hold on
plot(xDataNonSmooth{iTrajNS,2}(outdof,:),xDataNonSmooth{iTrajNS,2}(outdof+n,:),'Linewidth',2,'Color','k','DisplayName','FOM');
plot(xDataNonSmooth{iTrajNS,2}(outdof,xDataNonSmooth{iTrajNS,3}(1)),xDataNonSmooth{iTrajNS,2}(ndof+outdof,xDataNonSmooth{iTrajNS,3}(1)),'.','Color', 'k','MarkerSize',14,'LineWidth',10,'DisplayName','Switch regime - FOM');
plot(xDataNonSmooth{iTrajNS,2}(outdof,xDataNonSmooth{iTrajNS,3}(2:end)),xDataNonSmooth{iTrajNS,2}(ndof+outdof,xDataNonSmooth{iTrajNS,3}(2:end)),'.','Color', 'k','MarkerSize',14,'LineWidth',10,'HandleVisibility','off');

plot(eta(1,:),eta(2,:),'Linewidth', 2, 'Color','g','DisplayName','ROM');
plot(eta(1,indSwitch(1)), eta(2,indSwitch(1)),'.','Color', 'r','MarkerSize',14,'LineWidth',10,'DisplayName','Switch regime - ROM');
plot(eta(1,indSwitch(2:end)), eta(2,indSwitch(2:end)),'.','Color', 'r','MarkerSize',14,'LineWidth',10,'HandleVisibility','off');
% plot(fixedpointPlus(outdof), fixedpointPlus(outdof+n), 'b*', 'MarkerSize',10,'LineWidth',2);
% plot(fixedpointMinus(outdof), fixedpointMinus(outdof+n), 'r*', 'MarkerSize',10,'LineWidth',2);

grid on
xlabel('$q_{mid}$ [m]','interpreter','latex');
ylabel('$\dot{q}_{mid}$ [m/s]','interpreter','latex')
set(gca,'Fontsize',20)
ylim([-1 1]);

legend('location','NW','Interpreter','latex')

%% Non-autonomous system 
load trajectory_full

[v_minus,r_minus] = NonAutoParam_RD(n,V0Minus,W0Minus,AMinus,forcingVectors_na);
[v_plus,r_plus] = NonAutoParam_RD(n,V0Plus,W0Plus,APlus,forcingVectors_na);

IMInfoFMinus = IMInfoSMinus;
autParam_Minus = IMInfoFMinus.parametrization.map;
nonAutParam_Minus = @(t,q,fFreqs,fAmpls) autParam_Minus(q) + v_minus(t,fFreqs,fAmpls);
IMInfoFMinus.parametrization.map = nonAutParam_Minus;
RDInfoFMinus = RDInfoMinus;
autRedDyn_Minus = RDInfoFMinus.reducedDynamics.map;
nonAutRedDyn_Minus = @(t,x,fFreqs,fAmpls) autRedDyn_Minus(x) + r_minus(t,fFreqs,fAmpls);
RDInfoFMinus.reducedDynamics.map = nonAutRedDyn_Minus;

IMInfoFPlus = IMInfoSPlus;
autParam_Plus = IMInfoFPlus.parametrization.map;
nonAutParam_Plus = @(t,q,fFreqs,fAmpls) autParam_Plus(q) + v_plus(t,fFreqs,fAmpls);
IMInfoFPlus.parametrization.map = nonAutParam_Plus;
RDInfoFPlus = RDInfoPlus;
autRedDyn_Plus = RDInfoFPlus.reducedDynamics.map;
nonAutRedDyn_Plus = @(t,x,fFreqs,fAmpls) autRedDyn_Plus(x) + r_plus(t,fFreqs,fAmpls);
RDInfoFPlus.reducedDynamics.map = nonAutRedDyn_Plus;


rdof = 1; velTransp = v_ground; 
sFunction = @(x) x(rdof+1)-velTransp; velVec = zeros(1,2*rdof); velVec(rdof+1) = 1; 
DxsFunction = @(x) velVec;
redEqPlus = IMInfoFPlus.chart.map(fixedpointPlus);
redEqMinus = IMInfoFMinus.chart.map(fixedpointMinus);
fPlus_forced = @(t,eta) RDInfoFPlus.reducedDynamics.map(t,eta-redEqPlus,forFreq,forAmp);
fMinus_forced = @(t,eta) RDInfoFMinus.reducedDynamics.map(t,eta-redEqMinus,forFreq,forAmp);
fFunction_forced = @(t,x,l) 0.5*(1+l)*fPlus_forced(t,x)+0.5*(1-l)*fMinus_forced(t,x);

timeSpanForced = [solution_full{1,1}(1) solution_full{1,1}(end)];
ic_forced = solution_full{1,2}([outdof n+outdof],1);
[t_forced,eta_forced,indSwitch_forced,l0vect_forced] = odenonsmooth_1switchsurf(fFunction_forced, sFunction, DxsFunction, timeSpanForced, ic_forced);
eta_forced = transpose(eta_forced);

figure
hold on
plot(solution_full{1,2}(outdof,:), solution_full{1,2}(outdof+n,:), 'Linewidth', 2, 'Color','k','DisplayName','Full order model')
plot(eta_forced(1,:),eta_forced(2,:),'Linewidth', 2, 'Color','g','LineStyle','-','DisplayName','Reduced order model');
grid on
xlabel('$q_{mid}$ [m]','interpreter','latex');
ylabel('$\dot{q}_{mid}$ [m/s]','interpreter','latex')
set(gca,'Fontsize',20)
ylim([-0.25 0.25]);
legend('location','SW','Interpreter','latex')

figure 
hold on
plot(solution_full{1,1}, solution_full{1,2}(outdof,:), 'Linewidth', 2, 'Color','k','DisplayName','Full order model')
plot(t_forced,eta_forced(1,:),'Linewidth', 2, 'Color','g','LineStyle','--','DisplayName','Reduced order model');
grid on
xlabel('$t$ [s]','interpreter','latex');
ylabel('$q_{mid}$ [m]','interpreter','latex');
set(gca,'Fontsize',20)
xlim([0 0.2])
ylim([-3.5*1e-4 3.5*1e-4])
legend('location','SW','Interpreter','latex')


% Shaw-Pierre example: equation-driven SSM-based model reduction.

clearvars; figure; classicColors = colororder; close all; format shortg; clc

load sim_data
optFRC = 1;

%% linear part quantities

ndof = size(M,1);
Mass = M;

%%% negative case
KMinus = K+KIII;
AMinus = [zeros(ndof) eye(ndof); -M\KMinus -M\C];
[VMinusComplex,RMinusComplex] = eigSorted(AMinus);
WMinusComplex = inv(VMinusComplex);
VMinus = [];
for ii = 1:ndof
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

%%% positive case
KPlus = K+KIII;
APlus = [zeros(ndof) eye(ndof); -M\KPlus -M\C];
[VPlusComplex,RPlusComplex] = eigSorted(APlus);
WPlusComplex = inv(VPlusComplex);
VPlus = [];
for ii = 1:ndof
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

SSMDim = 2;
M = 3; %polynomial order of the Taylor approximation of the SSM

%% parametrization and reduced dynamics
q0 = abs(fixedpointMinus(1));

[phi,Exp_mat] = multivariatePolynomial(SSMDim,1,M);
optsGeomtery.l = 0;
optsGeomtery.c1 = 0;
optsGeomtery.c2 = 0;
l_opt = 0;
Err = [];

%%%%%% negative case (lambda = -1)

[z1_minus, z2_minus, coeff_param_minus, coeff_dyn_minus] = parametrization_real(RMinusComplex, VMinus, q0, gamma, -1, R0);
IMInfoMinus.chart.map = @(x) chartMap(shiftMinus(x));


%parametrization
H1_minus = [eye(SSMDim), zeros(SSMDim,size(coeff_param_minus,2));
           zeros(SSMDim, SSMDim) coeff_param_minus];
H_minus = VMinus * H1_minus;
IMParam_minus_H1 = @(q) H1_minus * phi(q);
IMParam_minus = @(q) H_minus * phi(q);
paramInfo_minus = struct('map',IMParam_minus,'polynomialOrder',M,...
            'dimension', SSMDim, 'tangentSpaceAtOrigin',H_minus(:,1:SSMDim),...
            'nonlinearCoefficients',H_minus(:,SSMDim+1:end),'phi',phi,...
            'exponents',Exp_mat,'l',optsGeomtery.l,...
            'c1',optsGeomtery.c1,'c2',optsGeomtery.c2);
IMInfoMinus.parametrization = paramInfo_minus;
paraMap_minus = IMInfoMinus.parametrization.map;
IMInfoMinus.parametrization.map = @(eta) fixedpointMinus + paraMap_minus(eta);

%reduced dynamics
W_r_minus = coeff_dyn_minus;
fun_reduced_minus = @(x) W_r_minus*phi(x);
R_info_minus = struct('map',fun_reduced_minus,'coefficients',W_r_minus,'polynomialOrder',...
    M,'phi',phi,'exponents',Exp_mat,'l_opt',l_opt,...
    'CV_error',Err);
T_info_minus = struct('map',@(x) x,'coefficients',eye(SSMDim),'polynomialOrder',...
                    M,'phi',@(x) x,'exponents',eye(SSMDim));
iT_info_minus = T_info_minus;
N_info_minus = R_info_minus;
[V_minus,D_minus,d_minus] = eigSorted(W_r_minus(:,1:SSMDim));
RDInfoMinus = struct('reducedDynamics',R_info_minus,'inverseTransformation',...
             iT_info_minus,'conjugateDynamics',N_info_minus,'transformation',T_info_minus,...
             'conjugacyStyle','default','dynamicsType','flow',...
             'eigenvaluesLinPartFlow',d_minus,'eigenvectorsLinPart',V_minus);


%%%%%%  positive case (lambda = 1)

[z1_plus, z2_plus, coeff_param_plus, coeff_dyn_plus] = parametrization_real(RPlusComplex, VPlus, q0, gamma, 1, R0);

%parametrization
IMInfoPlus.chart.map = @(x) chartMap(shiftPlus(x));
H1_plus = [eye(SSMDim), zeros(SSMDim,size(coeff_param_plus,2));
           zeros(SSMDim, SSMDim) coeff_param_plus];
H_plus = VPlus * H1_plus;
IMParam_plus_H1 = @(q) H1_plus * phi(q);
IMParam_plus = @(q) H_plus * phi(q);
paramInfo_plus = struct('map',IMParam_plus,'polynomialOrder',M,...
            'dimension', SSMDim, 'tangentSpaceAtOrigin',H_plus(:,1:SSMDim),...
            'nonlinearCoefficients',H_plus(:,SSMDim+1:end),'phi',phi,...
            'exponents',Exp_mat,'l',optsGeomtery.l,...
            'c1',optsGeomtery.c1,'c2',optsGeomtery.c2);
IMInfoPlus.parametrization = paramInfo_plus;
paraMap_plus = IMInfoPlus.parametrization.map;
IMInfoPlus.parametrization.map = @(eta) fixedpointPlus + paraMap_plus(eta);

%reduced dynamics
W_r_plus = coeff_dyn_plus;
fun_reduced_plus = @(x) W_r_plus*phi(x);
R_info_plus = struct('map',fun_reduced_plus,'coefficients',W_r_plus,'polynomialOrder',...
    M,'phi',phi,'exponents',Exp_mat,'l_opt',l_opt,...
    'CV_error',Err);
T_info_plus = struct('map',@(x) x,'coefficients',eye(SSMDim),'polynomialOrder',...
                    M,'phi',@(x) x,'exponents',eye(SSMDim));
iT_info_plus = T_info_plus;
N_info_plus = R_info_plus;
[V_plus,D_plus,d_plus] = eigSorted(W_r_plus(:,1:SSMDim));
RDInfoPlus = struct('reducedDynamics',R_info_plus,'inverseTransformation',...
             iT_info_plus,'conjugateDynamics',N_info_plus,'transformation',T_info_plus,...
             'conjugacyStyle','default','dynamicsType','flow',...
             'eigenvaluesLinPartFlow',d_plus,'eigenvectorsLinPart',V_plus);

%% Get SSM intersections

sFunctionPhys = @(x) x(ndof+1);

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
eta2_Max = 3;
eta2_vec = linspace(-eta2_Max, eta2_Max, 41);

for l0 = [1 -1]
    eta2_vec_plus = eta2_vec(eta2_vec>0);
    eta1_vecPlus = zeros(1,length(eta2_vec_plus));
    eta1_0 = 0;
    options = optimset('Display','off');
    for ii = 1:length(eta2_vec_plus)
        eta1_0 = fsolve(@(eta1) ...
            sFunctionPhys(vFunction(0,[eta1; eta2_vec_plus(ii)],l0)),eta1_0, options);
        eta1_vecPlus(ii) = eta1_0;
    end
    eta2_vec_minus = sort(-eta2_vec_plus);
    eta1_vecMinus = zeros(1,length(eta2_vec_minus));
    eta1_0 = 0;
    for ii = length(eta2_vec_minus):-1:1
        eta1_0 = fsolve(@(eta1) ...
            sFunctionPhys(vFunction(0,[eta1; eta2_vec_minus(ii)],l0)),eta1_0, options);
        eta1_vecMinus(ii) = eta1_0;
    end
    if l0 == 1
        gammaPlus = @(s) [interp1([eta2_vec_minus(1:end-1) eta2_vec_plus], [eta1_vecMinus(1:end-1) eta1_vecPlus],s); s];
    else
        gammaMinus = @(s) [interp1([eta2_vec_minus(1:end-1) eta2_vec_plus], [eta1_vecMinus(1:end-1) eta1_vecPlus],s); s];
    end
end

eta2_vec = [eta2_vec_minus(1:end-1) eta2_vec_plus];
gammaFunction = @(s,l) 0.5*(1+l)*gammaPlus(s)+0.5*(1-l)*gammaMinus(s);
WeightMatrix= eye(2*ndof);

%% Plot the SSMs and switching surface (physical coordinates)
eta1_vals = linspace(-3,3,501);
eta2_vals = linspace(-3,3,501);
[ETA1, ETA2] = meshgrid(eta1_vals, eta2_vals);

Y1 = zeros(size(ETA1));
Y2 = zeros(size(ETA1));
Z1 = zeros(size(ETA1));
Z2 = zeros(size(ETA1));
X1 = zeros(size(ETA1));
X2 = zeros(size(ETA1));
X3 = zeros(size(ETA1));
X4 = zeros(size(ETA1));

customFigure;
plot3(fixedpointMinus(1),fixedpointMinus(2),fixedpointMinus(3),'*','Color',...
    'red','MarkerSize',10,'linewidth',3)
plot3(fixedpointPlus(1),fixedpointPlus(2),fixedpointPlus(3),'*','Color',...
    'blue','MarkerSize',10,'linewidth',3)

manifold_plus = cell(3,1);
manifold_minus = cell(3,1);

for l = [1 -1]
    for ii = 1:size(ETA1,1)
        for jj = 1:size(ETA1,2) 
            if l == 1
                Y1(ii,jj) = ETA1(ii,jj);
                Y2(ii,jj) = ETA2(ii,jj);
                Z1(ii,jj) = [0 0 1 0] * IMParam_plus_H1([Y1(ii,jj); Y2(ii,jj)]);
                Z2(ii,jj) = [0 0 0 1] * IMParam_plus_H1([Y1(ii,jj); Y2(ii,jj)]);
                X1(ii,jj) = [1 0 0 0] * IMInfoPlus.parametrization.map([Y1(ii,jj); Y2(ii,jj)]);
                X2(ii,jj) = [0 1 0 0] * IMInfoPlus.parametrization.map([Y1(ii,jj); Y2(ii,jj)]);
                X3(ii,jj) = [0 0 1 0] * IMInfoPlus.parametrization.map([Y1(ii,jj); Y2(ii,jj)]);
                X4(ii,jj) = [0 0 0 1] * IMInfoPlus.parametrization.map([Y1(ii,jj); Y2(ii,jj)]);
            else
                Y_negative =  p2n_reduced_coordinates([ETA1(ii,jj),ETA2(ii,jj)],q0, VMinus);
                Y1(ii,jj) = Y_negative(1);
                Y2(ii,jj) = Y_negative(2);
                Z1(ii,jj) = [0 0 1 0] * IMParam_minus_H1([Y1(ii,jj); Y2(ii,jj)]);
                Z2(ii,jj) = [0 0 0 1] * IMParam_minus_H1([Y1(ii,jj); Y2(ii,jj)]);
                X1(ii,jj) = [1 0 0 0] * IMInfoMinus.parametrization.map([Y1(ii,jj); Y2(ii,jj)]);
                X2(ii,jj) = [0 1 0 0] * IMInfoMinus.parametrization.map([Y1(ii,jj); Y2(ii,jj)]);
                X3(ii,jj) = [0 0 1 0] * IMInfoMinus.parametrization.map([Y1(ii,jj); Y2(ii,jj)]);
                X4(ii,jj) = [0 0 0 1] * IMInfoMinus.parametrization.map([Y1(ii,jj); Y2(ii,jj)]);
            end
        end
    end

    if l == 1
        X3(find(X3<0)) = nan; 
        h = surf(X1, X2, X3);
        h.FaceColor = classicColors(1,:);
        h.EdgeColor = 'none';
        h.FaceAlpha = .5;

        manifold_plus{1} = X1;
        manifold_plus{2} = X2;
        manifold_plus{3} = X3;

    else
        X3(find(X3>0)) = nan; 
        h = surf(X1, X2, X3);
        h.FaceColor = classicColors(2,:);
        h.EdgeColor = 'none';
        h.FaceAlpha = .5;

        manifold_minus{1} = X1;
        manifold_minus{2} = X2;
        manifold_minus{3} = X3;
    end
end

% Switching surface
[XX, YY] = meshgrid([-1 1]*max(max(abs(X1))),[-1 1]*max(max(abs(X2))));
ZZ = 0*XX;
h = surf(XX,YY,ZZ);
h.FaceColor = classicColors(3,:);
h.EdgeColor = 'none';
h.FaceAlpha = .5;

xlabel('$q_1$ [m]','interpreter','latex')
ylabel('$q_2$ [m]','interpreter','latex')
zlabel('$\dot{q}_1$ [m/s]','interpreter','latex');

view(-10,11)

eta_intersection_minus = gammaMinus(eta2_vals);
x_intersection_minus = vFunction(0,eta_intersection_minus,-1);

eta_intersection_plus = gammaPlus(eta2_vals);
x_intersection_plus = vFunction(0,eta_intersection_plus,1);

plot3(x_intersection_minus(1,:), x_intersection_minus(2,:), x_intersection_minus(3,:), 'r', 'linewidth',3);
plot3(x_intersection_plus(1,:), x_intersection_plus(2,:), x_intersection_plus(3,:), 'b', 'linewidth',3);

xlim([-1.5 1.5])
ylim([-1.5 1.5])
zlim([-1.5 1.5])

grid on

set(gca,'FontSize',18)
set(get(gca,'XLabel'),'Rotation',0)
set(get(gca,'YLabel'),'Rotation',0)
set(get(gca,'ZLabel'),'Rotation',0)

plot3(x_full_DataNonSmooth{1,2}(1,:),x_full_DataNonSmooth{1,2}(2,:),x_full_DataNonSmooth{1,2}(3,:),'k','LineWidth',2)


%% Integrate smooth models separately (negative and positive)
xDataCurr = xData(2,1:2);
xDataCurrShifted = funToCell(xDataCurr,shiftPlus);
etaData = funToCell(xDataCurrShifted,chartMap);
etaRec_positive = advectRD(RDInfoPlus, etaData);
normedTrajDist = computeTrajectoryErrors(etaRec_positive, etaData);

plotTrajectories(etaData(1,:), etaRec_positive(1,:), 'm','PlotCoordinate', 1, 'DisplayName', {'Test set', 'Prediction'})
ylabel('$y_1^+$','interpreter','latex')
xlabel('$t $[s]','interpreter','latex')
xlim([0 xDataCurr{1,1}(end)]);
set(gca,'FontSize',20)

%% Integrate non-smooth model
idxIni = 1;
idxEnd = x_full_DataNonSmooth{1,3}(end);
timeSpan = [x_full_DataNonSmooth{1,1}(idxIni) x_full_DataNonSmooth{1,1}(idxEnd)+10.1];

% Integration 
[t,xRec,indSwitch,~] = integratePWSROM(wFunction,rFunction,vFunction,sFunctionPhys,timeSpan, x_full_DataNonSmooth{1,2}(:,idxIni),...
    'RICmethod','non_physical','WeightMatrix', WeightMatrix,'gammaStoEta',gammaFunction,'gammaEtatoS',sToEta);
xRec = transpose(xRec);

colorFOM = [0 0 0];
colorROM = classicColors(5,:);
customFigure('subPlot',[2 1]);tLimits = timeSpan([1 end]);%[x_full_DataNonSmooth{1,1}(idxIni) timeSpan(end)]; %[0 700];
subplot(211); xlim(tLimits)
plot(x_full_DataNonSmooth{1,1},x_full_DataNonSmooth{1,2}(1,:),'Linewidth',2,'Color', colorFOM);
plot(t,xRec(1,:),'Linewidth',2,'Color', colorROM);
plot(x_full_DataNonSmooth{1,1}(x_full_DataNonSmooth{1,3}),x_full_DataNonSmooth{1,2}(1,x_full_DataNonSmooth{1,3}),'.','Color', colorFOM,'MarkerSize',14);
plot(t(indSwitch),xRec(1,indSwitch),'.','Color', colorROM,'MarkerSize',10);
legend('FOM','ROM','S-FOM','S-ROM','numColumns',2)
xlabel('$t$ [s]','interpreter','latex')
ylabel('$q_1$ [m]','interpreter','latex')
subplot(212); xlim(tLimits)
plot(x_full_DataNonSmooth{1,1},x_full_DataNonSmooth{1,2}(ndof+1,:),'Linewidth',2,'Color', colorFOM);
plot(t,xRec(ndof+1,:),'Linewidth',2,'Color', colorROM);
plot(x_full_DataNonSmooth{1,1}(x_full_DataNonSmooth{1,3}),x_full_DataNonSmooth{1,2}(ndof+1,x_full_DataNonSmooth{1,3}),'.','Color', colorFOM,'MarkerSize',14);
plot(t(indSwitch),xRec(ndof+1,indSwitch),'.','Color', colorROM,'MarkerSize',10);
xlabel('$t$ [s]','interpreter','latex')
ylabel('$\dot{q}_1$ [m/s]','interpreter','latex');

%% Non-autonomous system
if optFRC == 1
    load frc_data
    % Define forced models
    [IMInfoFMinus,RDInfoFMinus] = forcedSSMROMgraph(IMInfoMinus,RDInfoMinus,...
        forcingVectors,WMinusComplex([SSMDim+1:end],:),...
        diag(RMinusComplex([SSMDim+1:end],[SSMDim+1:end])),...
        VMinusComplex(:,[SSMDim+1:end]),W0Minus,V0Minus);

    [IMInfoFPlus,RDInfoFPlus] = forcedSSMROMgraph(IMInfoPlus,RDInfoPlus,...
        forcingVectors,WPlusComplex([SSMDim+1:end],:),...
        diag(RPlusComplex([SSMDim+1:end],[SSMDim+1:end])),...
        VPlusComplex(:,[SSMDim+1:end]),W0Plus,V0Plus);

    freq0_vals = abs(RMinusComplex(1,1));
    amplitudeFunction = @(x) x(1,:); nPers = 140;
    customFigure('figNumber',1000);
    plotFRC(FRC_FOM,'color', 'k','datalabel','');
    xlim(detunings)
    ylim([0 1])
    legend('FOM','location','NW')
    xlabel('forcing frequency $\Omega$','Interpreter','latex');
    ylabel('amplitude $\max|q_1|$','Interpreter','latex');
    set(gca,'FontSize',20)

    % Check sticking condition
    check_function = @(x1, x3, x4,omega,t,epsilon) abs(-2 * stiffness/mass * x1 + stiffness/mass * x3 + damping/mass * x4 - gamma * x1.^3 + epsilon * forcingVectors(3) * cos(omega * t));

    for forAmp = forcingAmplitudes
        forFreq = .95*freq0_vals(1);
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
        x0 = real((1i*forFreq*eye(2*ndof)-AMinus)\forcingVectors)*forAmp;
        timeSpan = [0 nPers*2*pi/forFreq];
        % Integration
        [t,xRec,indSwitch,~] = integratePWSROM(wFunction,rFunction,vFunction,sFunctionPhys,timeSpan, x0,...
            'RICmethod','not_physical','WeightMatrix', WeightMatrix,'gammaStoEta',gammaFunction,'gammaEtatoS',sToEta);
        xRec = transpose(xRec);
        idxLastLoop = sum(t<(nPers-1)*2*pi/forFreq);
        amp = max(abs(amplitudeFunction(xRec(:,idxLastLoop:end))));
        IC = xRec(:,end);

        nstepsLow = 41;
        freqsLow = linspace(detunings(1),forFreq,nstepsLow);
        ICsLow = zeros(2*ndof,nstepsLow);
        ampsLow = zeros(1,nstepsLow);
        ICsLow(:,end) = IC; ampsLow(end) = amp;
        nstepsRef = 41;
        freqsHigh = linspace(forFreq,detunings(2),nstepsLow);
        nstepsHigh = length(freqsHigh);
        ICsHigh  = zeros(2*ndof,nstepsHigh);
        ampsHigh  = zeros(1,nstepsHigh);
        ICsHigh(:,1) = IC; ampsHigh(1) = amp;

        % Loop
        for iFreq = nstepsLow-1:-1:1
            forFreq = freqsLow(iFreq); x0 = ICsLow(:,iFreq+1);
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

            timeSpan = [0 nPers*2*pi/forFreq];
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

            timeSpan = [0 nPers*2*pi/forFreq];
            % Integration
            [t,xRec,indSwitch,~] = integratePWSROM(wFunction,rFunction,vFunction,sFunctionPhys,timeSpan, x0);
            xRec = transpose(xRec);
            idxLastLoop = sum(t<(nPers-1)*2*pi/forFreq);
            ampsHigh(iFreq) = max(abs(amplitudeFunction(xRec(:,idxLastLoop:end))));
            ICsHigh(:,iFreq) = xRec(:,end);

          if check_function(xRec(1,indSwitch), xRec(2,indSwitch), xRec(4,indSwitch), forFreq, t(indSwitch),forAmp) < delta
            plot(t(indSwitch)*forFreq/(2*pi), check_function(xRec(1,indSwitch), xRec(2,indSwitch), xRec(4,indSwitch), forFreq, t(indSwitch),forAmp),'*')
            hold on 
          end

        end
        %
        figure(1000);
        h = plot([freqsLow freqsHigh(2:end)],[ampsLow ampsHigh(2:end)],'.','MarkerSize',12,'Color', colorROM,'HandleVisibility','off');
        drawnow;
    end
    delete(h)
    plot([freqsLow freqsHigh(2:end)],[ampsLow ampsHigh(2:end)],'.','MarkerSize',12,'Color', colorROM,'DisplayName','ROM');
    legend('location','NW','Interpreter','latex')
end

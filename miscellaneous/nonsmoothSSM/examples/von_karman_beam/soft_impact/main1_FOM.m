% Von Karman beam with soft impact: data generation from full-order model

clearvars
close all
format shortg
clc

%% setup 

nElements = 4;
[M, C, K, fnl, fExt, outdof, PlotFieldonDefMesh] = buildModel(nElements);
n = size(M,1);
[F, lambda, ~, G, DG] = functionFromTensors(M, C, K, fnl); d_r1 = -real(lambda(1))/abs(lambda(1))*100
[W,A,V] = linearpart(M,C,K); % figure; spy(abs(W*A*V)>1e-4)

sSSMDim = 2;
Ve = V(:,1:sSSMDim); % Mode shape
We = W(1:sSSMDim,:); % Projection to mode shape
Ae = We*A*Ve; % Reduced, linearized dynamics

% Load and displacement vector: midpoint displacement
displacementVector = zeros(1,n); displacementVector(outdof) = 1;
loadVector = zeros(n,1); loadVector(outdof) = 1;  %  could also be set as modal ones

%% Static analysis to evaluate lin. vs nonlin. response

loadCoefficients = [0.01:0.01:0.09 0.1:0.1:0.9 1:1:10 12 14 15]*1e3;
nLoads = length(loadCoefficients);
wlin = (K\loadVector);
w0 = loadCoefficients(1)*wlin; % linear initial guess
displacementLinear = loadCoefficients*(displacementVector*wlin);
IC = zeros(2*n,nLoads); displacementNonlinear = zeros(1,nLoads);
options = optimoptions('fsolve', 'MaxFunctionEvaluations', 200*n, 'Display', 'off');
for iLoad = 1:nLoads
    f_eq = @(w)([zeros(n) M]*F(0,[w; zeros(n,1)]) + loadVector*loadCoefficients(iLoad));
    [w0, ~, exitflag, output] = fsolve(f_eq, w0, options);
    if exitflag <= 0
        error('Warning: No solution found for loading configuration')
    end
    IC(:,iLoad) = [w0; zeros(n,1)]; displacementNonlinear(iLoad) = displacementVector*w0;
end

customFigure(); classicColors = colororder;
plot(displacementLinear,loadCoefficients/1e3,'k--','Linewidth',1,'DisplayName','Linear')
plot(displacementNonlinear,loadCoefficients/1e3,'Color',classicColors(1,:),'Linewidth',2,'DisplayName','Nonlinear')
legend('location','SE')
xlim(abs([displacementNonlinear([1 end])]))
xlabel('Displacement [m]','Interpreter','latex');
ylabel('Force [kN]','Interpreter','latex');
title('Static loading analysis');
drawnow;

customFigure(); classicColors = colororder;
displacementDifference = abs((displacementNonlinear-displacementLinear))./abs(displacementLinear);
plot(displacementDifference*100,loadCoefficients/1e3,'Color',classicColors(1,:),'Linewidth',2)
xlim([0 max(displacementDifference)]*100)
xlabel('Relative displacement difference $\%$','Interpreter','latex');
ylabel('Force [kN]','Interpreter','latex');
title('Linear vs. Nonlinear')
drawnow;
% Pick up two initial trajectories that has high expected nonlinear content
indNonlinIC = length(loadCoefficients) - [2 1];
ICint = IC(:,indNonlinIC);
indTest = 1;
indTrain = 2;
ICint = [ICint -ICint(:,indTrain)];

% Define the linear regime at 0.05 % relative displacement
linearDisplacementReference = displacementNonlinear(sum(displacementDifference<(0.05/100))+1);
nonlinearDisplacementReference = max(displacementNonlinear(indNonlinIC));
desiredAmpltitudeDecay = nonlinearDisplacementReference/linearDisplacementReference;

%% Setup nonsmooth system
delta = 1e5; 
K_delta = zeros(n); K_delta(outdof, outdof) = delta;

[F_minus, lambda_minus, ~, G_minus, DG_minus] = functionFromTensors(M, C, K+K_delta, fnl);
[W_minus,A_minus,V_minus] = linearpart(M,C,K+K_delta); 

forceTrain = 12e3; forceTest = 9e3;

% Minus Side
% Equilibrium
loadCoef = delta;
f_eq = @(w)([zeros(n) M]*F_minus(0,[w; zeros(n,1)]));
[w0, ~, exitflag, output] = fsolve(f_eq, wlin*loadCoef, options);
fixedpointMinus = [w0; zeros(n,1)];
AMinus = A_minus + DG_minus(fixedpointMinus);
% IC Test
ICMinus = [];
loadCoef = forceTest;
f_eq = @(w)([zeros(n) M]*F(0,[w; zeros(n,1)]) + loadVector*loadCoef);
[w0, ~, exitflag, output] = fsolve(f_eq, wlin*loadCoef, options);
ICMinus = [ICMinus [w0; zeros(n,1)]];
loadCoef = forceTrain;
f_eq = @(w)([zeros(n) M]*F(0,[w; zeros(n,1)]) + loadVector*loadCoef);
[w0, ~, exitflag, output] = fsolve(f_eq, wlin*loadCoef, options);
ICMinus = [ICMinus [w0; zeros(n,1)]];

% Plus Side
% Equilibrium
f_eq = @(w)([zeros(n) M]*F(0,[w; zeros(n,1)]));
[w0, ~, exitflag, output] = fsolve(f_eq, wlin*loadCoef, options);
fixedpointPlus = [w0; zeros(n,1)];
APlus = A + DG(fixedpointPlus);
% IC Test
ICPlus = [];
loadCoef = forceTest;
f_eq = @(w)([zeros(n) M]*F(0,[w; zeros(n,1)])+ loadVector*loadCoef);
[w0, ~, exitflag, output] = fsolve(f_eq, wlin*loadCoef, options);
ICPlus = [ICPlus [w0; zeros(n,1)]];
loadCoef = forceTrain;
f_eq = @(w)([zeros(n) M]*F(0,[w; zeros(n,1)])+ loadVector*loadCoef);
[w0, ~, exitflag, output] = fsolve(f_eq, wlin*loadCoef, options);
ICPlus = [ICPlus [w0; zeros(n,1)]];

% Combination
loadCoefficient = delta;
sFunction = @(x) x(outdof); velVec = zeros(1,2*n); velVec(outdof) = 1;
DxsFunction = @(x) velVec;
fPlus = @(t,x) F(t,x);
fMinus = @(t,x) F_minus(t,x);
fFunction = @(t,x,l) 0.5*(1+l)*fPlus(t,x)+0.5*(1-l)*fMinus(t,x);

ICNS = ICPlus(:,1);

[V_minus, D_minus] = eigSorted(full(A_minus));
lambda_minus = diag(D_minus);

%% Integration non-smooth system

newMeasurement = true;
observable = @(x) x; % Observe the full phase space
slowTimeScale = 2*pi/abs(lambda(1));
fastTimeScale = 2*pi/abs(lambda(round(sSSMDim/2)));
numberPeriodsSlow = floor(log(desiredAmpltitudeDecay)/...
    (2*pi*(-real(lambda(1))/abs(lambda(1)))));
endTime = numberPeriodsSlow*slowTimeScale;
sliceInt = [5*slowTimeScale, endTime];
if newMeasurement
    % NonSmoothSims
    [t,x,indSwitch] = odenonsmooth_1switchsurf( fFunction, sFunction, DxsFunction, [0 endTime], ICNS(:,1), 'numDofs', n);
    x = transpose(x);
    idxIni = sum(t<sliceInt(1));
    idxEnd = sum(t<t(indSwitch(end))+sliceInt(1));
    t = transpose(t(idxIni:idxEnd));
    x = x(:,idxIni:idxEnd);
    indSwitch = indSwitch - idxIni + 1;
    indSwitch = transpose(indSwitch(indSwitch>0));
    xDataNonSmooth = cell(1,3);
    xDataNonSmooth{1,1} = t; xDataNonSmooth{1,2} = x; xDataNonSmooth{1,3} = indSwitch;
    % Set the sampling time to capture approximately 100 points per period
    % on the faster time scale
    numberPeriodsFast = floor(endTime/fastTimeScale);
    numberPointsPerPeriod = 100;
    nSamp = numberPeriodsFast*numberPointsPerPeriod+1;
    dt = endTime/(nSamp-1);
    % Integrate (full model) the positive and negative systems separately
    xDataMinus = integrateTrajectories(fMinus, endTime, ICMinus, nSamp, observable);
    xDataPlus = integrateTrajectories(fPlus, endTime, ICPlus, nSamp, observable);
    xDataMinus = sliceTrajectories(xDataMinus, sliceInt);
    xDataPlus = sliceTrajectories(xDataPlus, sliceInt);
    DataInfo = struct('nElements', nElements, 'loadvector', loadVector,'delta',delta);
    save('dataVKDecay2DFriction.mat', 'DataInfo', 'xDataNonSmooth','xDataMinus', 'xDataPlus', 'dt', 'endTime', 'nSamp','M','C','K','APlus','AMinus','fixedpointPlus','fixedpointMinus','outdof','options','lambda')
else
    load dataVKDecay2DFriction.mat
    if nElements ~= DataInfo.nElements
        error('The loaded data comes from a model with a different number of elements.')
    end
end

customFigure(); plotCoord =outdof;
plot(xDataMinus{1,1},xDataMinus{1,2}(plotCoord,:),'Linewidth',2,'DisplayName','Minus ')
plot(xDataPlus{1,1},xDataPlus{1,2}(plotCoord,:),'Linewidth',2,'DisplayName','Plus ')
plot(xDataNonSmooth{1,1},xDataNonSmooth{1,2}(plotCoord,:),'Linewidth',2,'DisplayName',['\delta = ' num2str(delta,'%0.e')])
legend('location','NE','numcolumns',1)
xlabel('$t$ [s]','interpreter','latex');
ylabel('$q_{mid}$ [m]','interpreter','latex');
set(gca,'FontSize',20)

%% plot during time and physical coordinates
customFigure();
plot(t, x(outdof,:),'k-','linewidth',2)
plot(t(indSwitch), x(outdof,indSwitch), 'Marker', 'o', 'linewidth', 2)
xlabel('$time$ [s]','Interpreter','latex')
ylabel('$q_{mid}$ [m]','interpreter','latex')
xlim([t(1) t(end)])
set(gca,'FontSize',20)
        
customFigure();
plot(x(outdof,:), x(outdof+n,:),'k-','linewidth',2)
plot([0 0], [x(outdof+n,indSwitch(1)) x(outdof+n,indSwitch(2))], 'r', 'linewidth', 3)
xlabel('$q_{mid}$ [m]','interpreter','latex');
ylabel('$\dot{q}_{mid}$ [m/s]','interpreter','latex');
set(gca,'FontSize',20)

%% Non-autonomous system

loadVector_na = zeros(n,1); loadVector_na(outdof) = 1;
forcingVectors_na = [zeros(n,1); M\loadVector_na];
freq_0 = abs(imag(lambda_minus(1)));

% Parameters
p0 = [freq_0; 25];
forAmp = p0(2);
forFreq = p0(1);
nPers = 200;

ampsLin = (K-eye(n)*p0(1)^2)\forcingVectors_na(n+1:end)*p0(2);
% Linear initial condition for numerical integration
ampLin = ampsLin(1);
currPhase = -pi/2;
x0 = [ampsLin; zeros(n,1)];

% Integration
timeSpan = [0 nPers*2*pi/p0(1)];
fPlus = @(t,x) F(t,x) + forAmp * forcingVectors_na * cos(forFreq * t);
fMinus = @(t,x) F(t,x)- loadCoefficient * x(5) * [zeros(n,1); M\loadVector] + forAmp * forcingVectors_na * cos(forFreq * t);
fFunction = @(t,x,l) 0.5*(1+l)*fPlus(t,x)+0.5*(1-l)*fMinus(t,x);    
[t,x,indSwitch] = odenonsmooth_1switchsurf(fFunction, sFunction, DxsFunction, timeSpan, x0, 'numDofs', n);
x = transpose(x);

% Plot
customFigure;
plot(t*p0(1)/(2*pi),x(outdof,:),'Linewidth',1)
ylabel('$q_{out}$','Interpreter','Latex');
xlabel('$t \Omega/(2\pi)$','Interpreter','Latex');
set(gca,'FontSize',20);

%% Forced response curves
loadVector_na = zeros(n,1); loadVector_na(outdof) = 1;
forcingVectors_na = [zeros(n,1); M\loadVector_na];
forcingAmplitudes = [25 35]; % forcing amplitudes
freq_0 = (abs(imag(lambda(1))) + abs(imag(lambda_minus(1))))/2;
detunings = freq_0*[0.9 1.1];
nPers = 150;
amplitudeFunction = @(x) x(outdof,:);
ii = 0;
FRC_full = cell(2,length(forcingAmplitudes));
for forAmp = forcingAmplitudes
    ii = ii+1;
    forFreq = freq_0;
    % Initial condition
    ampsLin = (K-eye(n)*forFreq^2)\forcingVectors_na(n+1:end)*forAmp;
    x0 = [ampsLin; zeros(n,1)];
    timeSpan = [0 nPers*2*pi/forFreq];
    % Integration
    fPlus = @(t,x) F(t,x) + forAmp * forcingVectors_na * cos(forFreq * t);
    fMinus = @(t,x) F(t,x)- loadCoefficient * x(5) * [zeros(n,1); M\loadVector] + forAmp * forcingVectors_na * cos(forFreq * t);
    fFunction = @(t,x,l) 0.5*(1+l)*fPlus(t,x)+0.5*(1-l)*fMinus(t,x);    
    [t,x,indSwitch] = odenonsmooth_1switchsurf(fFunction, sFunction, DxsFunction, timeSpan, x0, 'numDofs', n);
    x = transpose(x);
    idxLastLoop = sum(t<(nPers-1)*2*pi/forFreq);
    amp = max(abs(amplitudeFunction(x(:,idxLastLoop:end))));
    IC = x(:,end);

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
            
            % Integration
            fPlus = @(t,x) F(t,x) + forAmp * forcingVectors_na * cos(forFreq * t);
            fMinus = @(t,x) F(t,x)- loadCoefficient * x(5) * [zeros(n,1); M\loadVector] + forAmp * forcingVectors_na * cos(forFreq * t);
            fFunction = @(t,x,l) 0.5*(1+l)*fPlus(t,x)+0.5*(1-l)*fMinus(t,x);    
            [t,x,indSwitch] = odenonsmooth_1switchsurf(fFunction, sFunction, DxsFunction, timeSpan, x0, 'numDofs', n);
            x = transpose(x);
            idxLastLoop = sum(t<(nPers-1)*2*pi/forFreq);
            ampsLow(iFreq) = max(abs(amplitudeFunction(x(:,idxLastLoop:end))));
            ICsLow(:,iFreq) = x(:,end);                        
    end
    for iFreq = 2:nstepsHigh
            forFreq = freqsHigh(iFreq); x0 = ICsHigh(:,iFreq-1);
            timeSpan = [0 nPers*2*pi/forFreq];

            % Integration
            fPlus = @(t,x) F(t,x) + forAmp * forcingVectors_na * cos(forFreq * t);
            fMinus = @(t,x) F(t,x)- loadCoefficient * x(5) * [zeros(n,1); M\loadVector] + forAmp * forcingVectors_na * cos(forFreq * t);
            fFunction = @(t,x,l) 0.5*(1+l)*fPlus(t,x)+0.5*(1-l)*fMinus(t,x);    
            [t,x,indSwitch] = odenonsmooth_1switchsurf(fFunction, sFunction, DxsFunction, timeSpan, x0, 'numDofs', n);
            x = transpose(x);
            idxLastLoop = sum(t<(nPers-1)*2*pi/forFreq);
            ampsHigh(iFreq) = max(abs(amplitudeFunction(x(:,idxLastLoop:end))));
            ICsHigh(:,iFreq) = x(:,end);

    end
    FRC_full{1,ii} = [ampsLow ampsHigh(2:end)];
    FRC_full{2,ii} = [freqsLow freqsHigh(2:end)];

end

save('frc_data', 'FRC_full','loadVector_na','forcingVectors_na','forcingAmplitudes','freq_0','detunings','nPers','ICNS', 'delta')

% Forced response curves reversed

loadVector_na = zeros(n,1); loadVector_na(outdof) = 1;
forcingVectors_na = [zeros(n,1); M\loadVector_na];
forcingAmplitudes = [25 35]; % forcing amplitudes
freq_0 = (abs(imag(lambda(1))) + abs(imag(lambda_minus(1))))/2;
detunings = [freq_0 * 0.97 freq_0 * 1.05];
data_label = cell(1,length(forcingAmplitudes));
nPers = 150;
amplitudeFunction = @(x) x(outdof,:);
ii = 0;
FRC_full = cell(2,length(forcingAmplitudes));
for forAmp = forcingAmplitudes
    ii = ii+1;
    forFreq = detunings(2);
    % Initial condition
    ampsLin = (K-eye(n)*forFreq^2)\forcingVectors_na(n+1:end)*forAmp;
    x0 = [ampsLin; zeros(n,1)];
    timeSpan = [0 nPers*2*pi/forFreq];
    % Integration
    fPlus = @(t,x) F(t,x) + forAmp * forcingVectors_na * cos(forFreq * t);
    fMinus = @(t,x) F(t,x)- loadCoefficient * x(5) * [zeros(n,1); M\loadVector] + forAmp * forcingVectors_na * cos(forFreq * t);
    fFunction = @(t,x,l) 0.5*(1+l)*fPlus(t,x)+0.5*(1-l)*fMinus(t,x);    
    [t,x,indSwitch] = odenonsmooth_1switchsurf(fFunction, sFunction, DxsFunction, timeSpan, x0, 'numDofs', n);
    x = transpose(x);
    idxLastLoop = sum(t<(nPers-1)*2*pi/forFreq);
    amp = max(abs(amplitudeFunction(x(:,idxLastLoop:end))));
    IC = x(:,end);

    nstepsLow = 51;
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
            
            % Integration
            fPlus = @(t,x) F(t,x) + forAmp * forcingVectors_na * cos(forFreq * t);
            fMinus = @(t,x) F(t,x)- loadCoefficient * x(5) * [zeros(n,1); M\loadVector] + forAmp * forcingVectors_na * cos(forFreq * t);
            fFunction = @(t,x,l) 0.5*(1+l)*fPlus(t,x)+0.5*(1-l)*fMinus(t,x);    
            [t,x,indSwitch] = odenonsmooth_1switchsurf(fFunction, sFunction, DxsFunction, timeSpan, x0, 'numDofs', n);
            x = transpose(x);
            idxLastLoop = sum(t<(nPers-1)*2*pi/forFreq);
            ampsLow(iFreq) = max(abs(amplitudeFunction(x(:,idxLastLoop:end))));
            ICsLow(:,iFreq) = x(:,end);            
                      
    end

    FRC_full{1,ii} = ampsLow;
    FRC_full{2,ii} = freqsLow;
    FRC_full{3,ii} = delta;

end

FRC_full_reverse = FRC_full;
save('frc_data_reverse', 'FRC_full_reverse','delta')

% Von Karman beam with friction on a moving belt: data generation from full-order model

clearvars
close all
format shortg
clc

%% setup

nElements = 4;
[M, C, K, fnl, fExt, outdof, PlotFieldonDefMesh] = buildModel(nElements);
n = size(M,1);
[F, lambda, ~, G, DG] = functionFromTensors(M, C, K, fnl); d_r1 = -real(lambda(1))/abs(lambda(1))*100;
[W,A,V] =  linearpart(M,C,K); 

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
forceTrain = 12e2; forceTest = 9e2;
% friction law
alpha = 0.5;
beta = 0.5;
v_ground = 0.2;
iota = 10;
kappa = 5;

delta = 10;

f_ns_plus = @(t,x) delta * (1 + alpha/exp(1) * exp((beta - (x(outdof + n) - v_ground))/beta)) .* sum(any(sign(sign(x(outdof + n))+1),size(x,2)),1) + (-delta * (1 + alpha/exp(1) * exp((beta - (-x(outdof + n) - v_ground))/beta)) + 2 * delta * (1 + alpha/exp(1) * exp((beta + v_ground)/beta))) .* sum(~any(sign(sign(x(outdof + n))+1),size(x,2)),1); 
f_ns_minus = @(t,x) (-delta * (1 + alpha/exp(1) * exp((beta - (v_ground - x(outdof + n)))/beta)).* sum(~any(sign(sign(x(outdof + n) - 2 * v_ground)+1),size(x,2)),1)  + (+ delta * (1 + alpha/exp(1) * exp((beta - (v_ground - (-x(outdof + n) + 2* 2 * v_ground)))/beta)) - 2 * delta * (1 + alpha/exp(1) * exp((beta + v_ground)/beta))).* sum(any(sign(sign(x(outdof + n) - 2 * v_ground)+1),size(x,2)),1));
f_delta_plus = @(t,x) [zeros(outdof-1,1); f_ns_plus(t,x); zeros(n-outdof,1)];
f_delta_minus = @(t,x) [zeros(outdof-1,1); f_ns_minus(t,x); zeros(n-outdof,1)];
fNSPlus = @(t,x) [zeros(n,1); -M\f_delta_plus(t,x)];
fNSMinus = @(t,x) [zeros(n,1); -M\f_delta_minus(t,x)];
fPlus = @(t,x) F(t,x) + fNSPlus(t,x);
fMinus = @(t,x) F(t,x) + fNSMinus(t,x);

%%%%%%%%%% Find equilibria 
options = optimoptions('fsolve', 'MaxFunctionEvaluations', 200*n, 'Display', 'off');
% postive case
f_eq = @(w)[zeros(n) M]*fPlus(0,[w; zeros(n,1)]);
[w0, ~, ~, ~] = fsolve(f_eq, zeros(n,1), options);
fixedpointPlus = [w0; zeros(n,1)];

% negative case
f_eq = @(w)[zeros(n) M]*fMinus(0,[w; zeros(n,1)]);
[w0, ~, ~, ~] = fsolve(f_eq, zeros(n,1), options);
fixedpointMinus = [w0; zeros(n,1)];

% Combination
sFunction = @(x) x(outdof+n) - v_ground; velVec = zeros(1,2*n); velVec(outdof+n) = 1;
DxsFunction = @(x) velVec;

fFunction = @(t,x,l) 0.5*(1+l)*fPlus(t,x)+0.5*(1-l)*fMinus(t,x);

% IC Test Minus
ICMinus = [];
loadCoef = forceTest;
f_eq = @(w)([zeros(n) M]*F(0,[w; zeros(n,1)]) + loadVector*loadCoef);
[w0, ~, ~, ~] = fsolve(f_eq, wlin*loadCoef, options);
ICMinus = [ICMinus [w0; zeros(n,1)]];
loadCoef = forceTrain;
f_eq = @(w)([zeros(n) M]*F(0,[w; zeros(n,1)]) + loadVector*loadCoef);
[w0, ~, ~, ~] = fsolve(f_eq, wlin*loadCoef, options);
ICMinus = [ICMinus [w0; zeros(n,1)]];


% IC Test Plus
ICPlus = [];
loadCoef = -forceTest;
f_eq = @(w)([zeros(n) M]*F(0,[w; zeros(n,1)])+ loadVector*loadCoef);
[w0, ~, ~, ~] = fsolve(f_eq, wlin*loadCoef, options);
ICPlus = [ICPlus [w0; zeros(n,1)]];
loadCoef = -forceTrain;
f_eq = @(w)([zeros(n) M]*F(0,[w; zeros(n,1)])+ loadVector*loadCoef);
[w0, ~, exitflag, output] = fsolve(f_eq, wlin*loadCoef, options);
ICPlus = [ICPlus [w0; zeros(n,1)]];
ICNS = ICPlus(:,1);

F_nl_extended = DG(fixedpointPlus);
h = 1e-2;
f_ns_plus_dq = @(dq) delta * (1 + alpha/exp(1) * exp((beta - (dq - v_ground))/beta)) .* sum(any(sign(sign(dq)+1),size(dq,2)),1) + (-delta * (1 + alpha/exp(1) * exp((beta - (-dq - v_ground))/beta)) + 2 * delta * (1 + alpha/exp(1) * exp((beta + v_ground)/beta))) .* sum(~any(sign(sign(dq)+1),size(dq,2)),1); 
f_ns_minus_dq = @(dq) (-delta * (1 + alpha/exp(1) * exp((beta - (v_ground - dq))/beta)).* sum(~any(sign(sign(dq - 2 * v_ground)+1),size(dq,2)),1)  + (+ delta * (1 + alpha/exp(1) * exp((beta - (v_ground - (-dq + 2* 2 * v_ground)))/beta)) - 2 * delta * (1 + alpha/exp(1) * exp((beta + v_ground)/beta))).* sum(any(sign(sign(dq - 2 * v_ground)+1),size(dq,2)),1));
df_plus_fp = (f_ns_plus_dq(h) - f_ns_plus_dq(-h))/(2 * h);
df_minus_fp = (f_ns_minus_dq(h) - f_ns_minus_dq(-h))/(2 * h);
C_fr_plus = -M\[zeros(outdof-1,n); [zeros(1,outdof-1) df_plus_fp zeros(1,n-outdof)]; zeros(n-outdof,n)];
C_fr_minus = -M\[zeros(outdof-1,n); [zeros(1,outdof-1) df_minus_fp zeros(1,n-outdof)]; zeros(n-outdof,n)];
C_fr_plus_extended = [zeros(n,2*n);
                      zeros(n,n) C_fr_plus];
C_fr_minus_extended = [zeros(n,2*n);
                      zeros(n,n) C_fr_minus];
APlus = A + F_nl_extended + C_fr_plus_extended;
AMinus = A + F_nl_extended + C_fr_minus_extended;

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
    [t,x,indSwitch,l0vect] = odenonsmooth_1switchsurf( fFunction, sFunction, DxsFunction, [0 endTime], ICNS(:,1), 'numDofs', n);
    x = transpose(x);
    idxIni = sum(t<sliceInt(1));
    idxEnd = sum(t<t(indSwitch(end))+sliceInt(1));
    idxEnd = length(t);
    t = transpose(t(idxIni:idxEnd));
    l0vect = transpose(l0vect(idxIni:idxEnd));
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
    save('dataVKDecay2DFriction.mat', 'DataInfo', 'xDataNonSmooth','xDataMinus', 'xDataPlus', 'dt', 'endTime', 'nSamp','M','C','K','APlus','AMinus','fixedpointPlus','fixedpointMinus','outdof','options','lambda','v_ground')
else
    load dataVKDecay2DFriction.mat
    if nElements ~= DataInfo.nElements
        error('The loaded data comes from a model with a different number of elements.')
    end
end

%%
figure
hold on 
grid on 
plot(t, x(outdof,:),'k-','linewidth',2)
xlabel('$time$ [s]','Interpreter','Latex');
ylabel('$q_{mid}$ [m]','Interpreter','Latex');
set(gca,'FontSize',20)

figure
hold on
plot(x(outdof,:), x(outdof+n,:),'k-','linewidth',2,'HandleVisibility','off')
plot(x(outdof,indSwitch(1)), x(outdof+n,indSwitch(1)), 'g*','LineWidth',2,'DisplayName','Crossing condition verified')
plot(x(outdof,indSwitch(2:end)), x(outdof+n,indSwitch(2:end)), 'g*','LineWidth',2,'HandleVisibility','off')
grid on
xlabel('$q_{mid}$ [m]','interpreter','latex');
ylabel('$\dot{q}_{mid}$ [m/s]','interpreter','latex')
set(gca,'Fontsize',20)

cross_function_1 = @(t,x) (DxsFunction(x) * fFunction(t,x,1));
cross_function_2 = @(t,x) (DxsFunction(x) * fFunction(t,x,-1));

ylim([-1 1]);

index_start = 1;
kk = 0;
for ii = 1:length(indSwitch)
    index = indSwitch(ii);
    time = t(index);
    if cross_function_1(time,x(:,index)) * cross_function_2(time,x(:,index)) < 0
        kk = kk+1;
        if kk == 1
            plot(x(outdof,index), x(outdof+n,index),'r*','linewidth',2,'DisplayName','Sticking condition verified')
        else
            plot(x(outdof,index), x(outdof+n,index),'r*','linewidth',2,'HandleVisibility','off')
        end
    end
    index_start = index;
end

legend('location','NE','Interpreter','Latex')

%% Forced System 

T_lc = t(indSwitch(end)) - t(indSwitch(end-2));
omega_lc = 2 * pi/ T_lc;
x0_lc = x(:,indSwitch(end-2));
loadVector_na = zeros(n,1); loadVector_na(outdof) = 1;
forcingVectors_na = [zeros(n,1); M\loadVector_na];
max_vel_vect = [];
freqVec = 1.125*omega_lc;

theta_vect = cell(length(freqVec),1);
phi_vect = cell(length(freqVec),1);
x_pm_vect = cell(length(freqVec),1);
dx_pm_vect = cell(length(freqVec),1);

nPers = 200;
endTime = nPers*slowTimeScale;

for ii = 1:length(freqVec)  
    % Parameters
    p0 = [freqVec(ii); 0.5*1e1];
    forAmp = p0(2);
    forFreq = p0(1);
    x0 = x0_lc;
    
    % Integration
    timeSpan = [0 endTime];
    fPlus_forced = @(t,x) fPlus(t,x) + forAmp * forcingVectors_na * cos(forFreq * t);
    fMinus_forced = @(t,x) fMinus(t,x) + forAmp * forcingVectors_na * cos(forFreq * t);
    fFunction_forced = @(t,x,l) 0.5*(1+l)*fPlus_forced(t,x)+0.5*(1-l)*fMinus_forced(t,x);    
    [t_forced,x_forced,indSwitch_forced,l0vect] = odenonsmooth_1switchsurf(fFunction_forced, sFunction, DxsFunction, timeSpan, x0, 'numDofs', n);
    x_forced = transpose(x_forced);
     
    
    figure
    hold on 
    plot(x_forced(outdof,:), x_forced(outdof+n,:),'k-','linewidth',1)
    plot(x_forced(outdof,indSwitch_forced), x_forced(outdof+n,indSwitch_forced),'g*','linewidth',2)
%     plot(fixedpointPlus(outdof), fixedpointPlus(outdof+n), 'b*', 'MarkerSize',10);
%     plot(fixedpointMinus(outdof), fixedpointMinus(outdof+n), 'r*', 'MarkerSize',10);
    plot(x_forced(outdof,1), x_forced(outdof+n,1),'r*','linewidth',2)

    grid on 
    xlabel('$q_{mid}$ [m]','interpreter','latex');
    ylabel('$\dot{q}_{mid}$ [m/s]','interpreter','latex')
    set(gca,'Fontsize',20)
    
    max_vel = max(x_forced(outdof+n,round(length(t_forced)/2):end));   
    max_vel_vect = [max_vel_vect max_vel];

    theta_local = atan2(x_forced(outdof+n,:),x_forced(outdof,:)*mean(x_forced(outdof+n,:))/mean(x_forced(outdof,:))*1e3) + pi;
    phi_local = mod(t_forced.' * p0(1),2*pi);
    theta_vect{ii,1} = theta_local;
    phi_vect{ii,1} = phi_local;

end

solution_full = cell(1,2);
solution_full{1,1} = t_forced;
solution_full{1,2} = x_forced;
save('trajectory_full','solution_full','forFreq', 'forAmp', 'delta','loadVector_na','timeSpan','forcingVectors_na');





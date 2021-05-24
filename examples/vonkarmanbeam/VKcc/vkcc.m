%% Finding a 2D SSM for a von Kármán beam
% 
% 
% This is an example of how to reconstruct a slow 2D SSM of a mechanical system 
% using synthetic measurements of a scalar quantity. In this example, we measure 
% the middle point displacement of a clamped-clamped von Kármán beam.

clearvars
close all
clc
%% Example setup

nElements = 12;
[M, C, K, fnl, f_vec, outdof, PlotFieldonDefMesh] = build_model(nElements);
n = size(M,1);    % mechanical dofs (axial def, transverse def, angle)
[F, lambda] = functionFromTensors(M, C, K, fnl);
%% Generation of Synthetic Data

loads = [2];
nTraj = size(loads, 2);
indTest = [1];
indTrain = 1;
%% 

loadvector = loads.*f_vec;
IC = getStaticResponse(K, M, F, loadvector, 1, PlotFieldonDefMesh);

%% 

new_meas = 0;
observable = @(x) x(outdof,:);
if new_meas == 1
    tEnd = 30;
    nSamp = fix(50 * tEnd * abs(imag(lambda(1))) / (2*pi));
    dt = tEnd/(nSamp-1);
    tic
    FullTrajectories = integrateTrajectories(F, @(x) x, tEnd, nSamp, nTraj, IC, 'odetol', 1e-4);
    toc
    DataInfo = struct('nElements',nElements,'loadvector',loadvector);
    save('data_VKcc.mat','DataInfo','FullTrajectories','dt','tEnd','nSamp')
else
    load data_VKcc.mat
    if nElements ~= DataInfo.nElements
       error('The loaded data comes from a model with a different number of elements.') 
    end
end

xData = cell(size(FullTrajectories,1),2);
for ii = 1:size(FullTrajectories,1)
    xData{ii,1} = FullTrajectories{ii,1};
    xData{ii,2} = observable(FullTrajectories{ii,2});
end
%% Delay embedding

SSMDim = 2;
overEmbed = 0;
yData = coordinates_embedding(xData, SSMDim, 'OverEmbedding', overEmbed);
%% 

sliceInt = [1, tEnd];
yDataTrunc = sliceTrajectories(yData, sliceInt);
%% Datadriven manifold fitting

SSMOrder = 1;
[V, SSMFunction, mfdInfo] = IMparametrization(yDataTrunc(indTrain,:), SSMDim, SSMOrder);

etaData = getProjectedTrajs(yData, V);
etaDataTrunc = getProjectedTrajs(yDataTrunc, V);
%% Reduced order model
ROMOrders = [7];
for iOrder = 1:length(ROMOrders)
ROMOrder = ROMOrders(iOrder);
[~,Tinv,N,T,NormalFormInfo] = IMdynamics_flow(etaDataTrunc(indTrain,:), ...
    'R_PolyOrd', ROMOrder, 'style', 'normalform', 'l_vals', [0,1e-4,1e-2,1e0], 'n_folds', 5);

zData = transformComplex(Tinv, etaData);
zDataTrunc = transformComplex(Tinv, etaDataTrunc);
[zRec, yRec] = integrateFlows(N, zDataTrunc, @(q) SSMFunction(T(q)));
etaRec = transformComplex(T, zRec);

%%
[reducedTrajDist, fullTrajDist] = computeRecDynErrors(etaRec, yRec, etaDataTrunc, yDataTrunc);
RRMSE_normal = mean(fullTrajDist(indTest))

plotReducedCoords(etaDataTrunc(indTest(1),:), etaRec(indTest(1),:))
legend({'Test set (truncated)', 'Prediction'})

plotReconstructedTrajectory(yData(indTest(1),:), yRec(indTest(1),:), 1, 'm')
legend({'Test set', 'Prediction'}); ylabel('$u \, [$m$]$','Interpreter','latex')
%% Backbone curves

N_info = NormalFormInfo.N;
[damp, freq] = polarnormalform(N_info.coeff, N_info.exponents, N_info.phi);
figure
rhoCal = abs(zData{indTest(1),2}(1,1));
[dmp, frq, amp, rho_plot] = backbonecurves(damp, freq, SSMFunction, T, 1, rhoCal);
subplot(121); ylabel('$u \, [$m$]$','Interpreter','latex')
subplot(122); ylabel('$u \, [$m$]$','Interpreter','latex')
%% 
% The data-driven model can now be used for forced response predictions. We 
% first compute forced response using the full model with SSMTool that are used, 
% after calibration, to validate the predictions of our data-driven reduced-order 
% model.

% Compute with SSMTool
f_full = 1e-3*[8];
w_span = abs(imag(lambda(1)))*[0.9 1.3];
FRC_full = getFRC_full(M, C, K, fnl, f_vec, f_full, outdof, w_span, ROMOrder);

%% Numerical Integration
if new_meas == 1
    f_sweep = f_full(end)*(M\f_vec);
    FRC_NI = integrateSweep(F, [33,40], n, f_sweep, observable, 'stepRevolutions', 300, 'step', 0.04);
    DataInfo = struct('nElements',nElements,'FRC_force',f_sweep,'FRC_ampcoord',outdof);
    save('data_sweep.mat','DataInfo','FRC_NI')
else
    load('data_sweep.mat','DataInfo','FRC_NI')
end


%%
% Calibration based on the maximal amplitude response of a validation FRC
integrationBasedCalibration = 0;

if integrationBasedCalibration
    Omega = 32;
    F_force = @(t,x,w) F(t,x) + [zeros(n,1); f_full(end)*(M\f_vec)*cos(w*t)];
    opts = odeset('RelTol', 1e-4);
    [t_sim,x_sim] = ode15s(@(t,x) F_force(t,x,Omega),0:dt:3*tEnd, zeros(2*n,1),opts);
    xCal = {t_sim(end-1000:end,:).', observable(x_sim(end-1000:end,:).')};
    yCalFull = coordinates_embedding(xCal, SSMDim, 'OverEmbedding', overEmbed);
    [~,indMax] = max(yCalFull{1,2}(ceil(0.5*size(V,1)),:));
    yCal = yCalFull{1,2}(:,indMax);
else
    % calibrate at FRC point
    Omega(1) = FRC_NI.omega(30);
    uCal(1) = FRC_NI.amp(30);
    yCal(:,1) = uCal*cos(Omega*dt*((1:size(V,1))-ceil(0.5*size(V,1)))).';
end

f_red = calibrateFRC(yCal, Omega, V, Tinv, damp, freq)


%% Compute FRC analytically
FRC_data = computeFRC(f_red, damp, freq, SSMFunction, T, @(y)max(abs(y)));

%% Plot
figure; hold on; colors = colororder;
plot(frq, amp,'k','DisplayName', 'Backbone - SSMlearn')
plotFRC(FRC_data, colors(1,:), 'SSMLearn')
plotFRC(FRC_full, colors(2,:), 'SSMTool')
% plot(OmegaFRC, uFRC, 'Color', colors(3,:), 'LineWidth', 2, 'DisplayName', 'SSMLearn')
scatter1 = scatter(FRC_NI.omega,FRC_NI.amp,48,'b','MarkerFaceColor','c','DisplayName','Numerical integration');
scatter1.MarkerFaceAlpha = 0.6;
scatter1.MarkerEdgeAlpha = 1.0;
xlim(w_span)
ylim([0,1.06*max(FRC_full.F1.Amp)])
xlabel('Excitation frequency $\Omega$ [rad/s]', 'Interpreter', 'latex')
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
legend

save(['FRC_data_', num2str(ROMOrder)], 'FRC_data')
save(['FRC_full_', num2str(ROMOrder)], 'FRC_full')
save(['backbone_data', num2str(ROMOrder)], 'frq', 'amp')
end
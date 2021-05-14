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

nElements = 6;
[M, C, K, fnl, f_vec, outdof, PlotFieldonDefMesh] = build_model(nElements);
n = size(M,1);    % mechanical dofs (axial def, transverse def, angle)
[F, lambda] = functionFromTensors(M, C, K, fnl);
%% Generation of Synthetic Data

loads = [1.75 2];
nTraj = size(loads, 2);
indTest = [1];
indTrain = setdiff(1:nTraj, indTest);
%% 

loadvector = loads.*f_vec;
IC = getStaticResponse(K, M, F, loadvector, 1, PlotFieldonDefMesh);

%% 

new_meas = 1;
observable = @(x) x(outdof,:);
if new_meas == 1
    tEnd = 30;
    nSamp = fix(50 * tEnd * abs(imag(lambda(1))) / (2*pi));
    dt = tEnd/(nSamp-1);
    tic
    FullTrajectories = integrateTrajectories(F, @(x) x, tEnd, nSamp, nTraj, IC);
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
overEmbed = 30;
yData = coordinates_embedding(xData, SSMDim, 'OverEmbedding', overEmbed);
%% 

sliceInt = [1, tEnd];
yDataTrunc = sliceTrajectories(yData, sliceInt);
%% Datadriven manifold fitting

SSMOrder = 1;
[V, SSMFunction, mfdInfo] = IMparametrization(yDataTrunc(indTrain,:), SSMDim, SSMOrder);
%% Plot and validation

etaData = getProjectedTrajs(yData, V);
etaDataTrunc = getProjectedTrajs(yDataTrunc, V);
%% Reduced order model

ROMOrder = 7;
[~,Tinv,N,T,NormalFormInfo] = IMdynamics_flow(etaDataTrunc(indTrain,:),'R_PolyOrd',ROMOrder,'style', 'normalform');

zData = transformComplex(Tinv, etaData);
zDataTrunc = transformComplex(Tinv, etaDataTrunc);
[zRec, yRec] = integrateFlows(N, zDataTrunc, @(q) SSMFunction(T(q)));
etaRec = transformComplex(T, zRec);

[reducedTrajDist, fullTrajDist] = computeRecDynErrors(etaRec, yRec, etaDataTrunc, yDataTrunc);
RRMSE_normal = mean(fullTrajDist(indTest))
%% Backbone curves

N_info = NormalFormInfo.N;
[damp, freq] = polarnormalform(N_info.coeff, N_info.exponents, N_info.phi);
figure
rhoCal = abs(zData{indTest(1),2}(1,1));
[dmp, frq, amp, rho_plot] = backbonecurves(damp, freq, SSMFunction, T, 1, rhoCal,'norm');
subplot(121); ylabel('$u \, [$m$]$','Interpreter','latex')
subplot(122); ylabel('$u \, [$m$]$','Interpreter','latex')
%% 
% The data-driven model can now be used for forced response predictions. We 
% first compute forced response using the full model with SSMTool that are used, 
% after calibration, to validate the predictions of our data-driven reduced-order 
% model.

% Compute with SSMTool
f_full = 1e-3*[2 8];
w_span = abs(imag(lambda(1)))*[0.9 1.2];
FRC_full = getFRC_full(M, C, K, fnl, f_vec, f_full, outdof, w_span, ROMOrder);

%%
% Calibration based on the maximal amplitude response of a validation FRC
idx_f_full = 1;
amp_max = max(FRC_full.(['F' num2str(idx_f_full)]).Amp);
[~,pos] = min(abs(amp-amp_max));
rho_max = rho_plot(pos);
f_red = abs(damp(rho_max)*rho_max) * f_full/f_full(idx_f_full);
ratio_force_data_full = f_red/f_full(idx_f_full)

calibrationLoads = [0.000001];
calloadvector = calibrationLoads.*f_vec;
ICCal = getStaticResponse(K, M, F, calloadvector, 0, 0);
zCal = Tinv(V.'*repmat(ICCal(outdof,indTest),size(yData{1,2},1),1));
rhoCal = abs(zCal(1));
fpsi = fsolve(@(fpsi) [damp(rhoCal)*rhoCal + fpsi(1)*sin(fpsi(2)); ...
    freq(rhoCal) - 0 + fpsi(1)/rhoCal*cos(fpsi(2))],...
    [0; 0]);
ratio = 0.5*abs(fpsi(1)/calibrationLoads(indTest));
imag(lambda(1))*abs(NormalFormInfo.iT.lintransf*sum(V).'*(f_vec'*(K\f_vec)))
% f_red = ratio*f_full;

%%
% Compute with data-driven model
ddROM = struct('Dim', SSMDim, 'Param', SSMFunction, 'CCtoNormal', Tinv, ...
    'ReducedDynNormal', N, 'CCfromNormal', T);
FRC_data = getFRC_ddROM(ddROM,f_red,w_span,1);

%%
% Plot
figure(100); hold on;
plot(frq, amp,'k','DisplayName', 'Backbone - SSMlearn')
plot_FRC(FRC_full, 'r', 'SSMTool')
plot_FRC(FRC_data, 'b', 'SSMLearn')

%% Numerical Integration
if new_meas == 1
    f_sweep = f_full(end)*(M\f_vec);
    FRC_NI = integrateSweep(F, w_span, n, f_sweep, observable);
    DataInfo = struct('nElements',nElements,'FRC_force',f_sweep,'FRC_ampcoord',outdof);
    save('data_sweep.mat','DataInfo','FRC_NI')
else
    load('data_sweep.mat','DataInfo','FRC_NI')
end
%%
figure(100);
plot(FRC_NI.omega,FRC_NI.amp,'c.','MarkerSize',24,'DisplayName','Numerical integration')

%% Learning a 4D SSM with fastSSMplus from experimental resonant vibration data
% We analyze the dynamics along a slow 4D resonant SSM. This example uses velocity 
% measurements from laser scanner vibrometry [1].

clearvars
close all
%
% [1] M. Cenedese, J. Axås, H. Yang, M. Eriten & G. Haller. (2022). Data-driven nonlinear 
% model reduction to spectral submanifolds in mechanical systems. Phil.
% Trans. R. Soc. A. 380:20210194. https://doi.org/10.1098/rsta.2021.0194

% Example setup
% The data consists of twelve decaying velocity signals measured on the tip 
% of the inner beam, excite via hammer impacts. These occur on three positions 
% (4 repetitions per impact) along the inner beam, spaced from the joint to the 
% free end. An initial data inspection reveals two fundamental resonant frequencies, 
% hereby justifying the identification of the slow 4D SSM of the system. Additional 
% information can be found in [1].
% 

load('resonantbeamdata.mat')
nTraj = size(xData,1);
for iTraj = 1:nTraj
    dt = mean(xData{iTraj,1}(2:end)-xData{iTraj,1}(1:end-1));
    xData{iTraj,1} = dt*(0:(size(xData{iTraj,1},2)-1));
end
%% Plot data
trajNames = {'(1,1)', '(1,2)', '(1,3)', '(1,4)', '(2,1)', '(2,2)', '(2,3)', '(2,4)', '(3,1)', '(3,2)', '(3,3)', '(3,4)'};
customFigure;
for iTraj = 1:nTraj
    plot(xData{iTraj,1},xData{iTraj,2},'Linewidth',1,'displayname',['Trajectory ',trajNames{iTraj}])
end
xlabel('time [s]', 'interpreter','latex')
ylabel('velocity [m/s]', 'interpreter','latex')
legend('interpreter','latex')

% Train and test
indTest = [2 6 10];
indTrain = setdiff(1:nTraj,indTest);

%% Delay embedding
dt = xData{1,1}(2)-xData{1,1}(1);
SSMDim = 4;

lag = round(pi/2/1531/dt);
overemb = 12 - 2*SSMDim-1;
[yData, opts_embd] = coordinatesEmbedding(xData, SSMDim, 'OverEmbedding', overemb, 'ShiftSteps', lag);

plotinds = [1,4,7];
customFigure();
colors = [0.1,0.1,0.7; 0.7,0.1,0.1; 0.1,0.7,0.1];
iColor = 1;
for iTraj = indTest
    plot3(yData{iTraj,2}(plotinds(1),:), yData{iTraj,2}(plotinds(2),:), yData{iTraj,2}(plotinds(3),:), 'Color', colors(iColor,:));
    iColor = iColor + 1;
end
xlabel('$s(t)$', 'interpreter', 'latex'); 
ylabel(['$s(t+',sprintf('%6.4f',lag*dt*(plotinds(2)-1)),')$'], 'interpreter', 'latex'); 
zlabel(['$s(t+',sprintf('%6.4f',lag*dt*(plotinds(3)-1)),')$'], 'interpreter', 'latex'); 
view(137,22); legend('off'); colororder(jet(nTraj))

%% Train model
SSMOrder = 3;
ROMOrder = 3;
NFOrder = 3;
[Mmap, iMmap, Tmap, iTmap, Nflow, Yrec] = fastSSMplus(yData(indTrain,:), SSMDim, SSMOrder, ROMOrder, NFOrder);

zData = transformTrajectories(@(y) iTmap(iMmap(y)), yData);
zRec = integrateFlows(@(z) Nflow(0,z), zData);
yRec = transformTrajectories(@(z) Mmap(Tmap(z)), zRec);

trajDist = computeTrajectoryErrors(yRec, yData, [1]);
NMTE = mean(trajDist(indTest))
%% Evaluation, prediction
indPlot = indTest(1);
customFigure();
plot(yData{indPlot,1}(1,:), yData{indPlot,2}(1,:), 'Color', [0,0,0], 'LineWidth', 0.5, 'DisplayName', ['Test trajectory ', trajNames{indPlot}]);
plot(yRec{indPlot,1}(1,:), yRec{indPlot,2}(1,:), '--', 'Color', [0.2,0.9,0.2], 'LineWidth', 0.5, 'DisplayName', 'Reconstruction');
xlim([yData{indPlot,1}(1), yData{indPlot,1}(end)])
ylim([-2.4, 2.4])
xlabel('time [s]', 'Interpreter', 'latex');
ylabel('velocity [m/s]','Interpreter','latex');
legend

customFigure();
plot(yData{indPlot,1}(1,:), yData{indPlot,2}(1,:), 'Color', [0,0,0], 'LineWidth', 1.5, 'DisplayName', ['Test trajectory ', trajNames{indPlot}]);
plot(yRec{indPlot,1}(1,:), yRec{indPlot,2}(1,:), '--', 'Color', [0.2,0.9,0.2], 'LineWidth', 1.5, 'DisplayName', 'Reconstruction');
xlim([0.2, 0.3])
ylim([-1.2, 1.2])
xlabel('time [s]', 'Interpreter', 'latex');
ylabel('velocity [m/s]','Interpreter','latex');

%% Decays along the SSM and instantaneous damping
% We plot the trajectories decay on the $(\rho_1,\rho_2)$ plane, splitting hence 
% the slow and fast directions.

zRec = transformTrajectories(@(y)iTmap(iMmap(y)), yRec);
customFigure; colororder(jet(nTraj))
for iTraj = 1:nTraj
    plot(abs(zRec{iTraj,2}(1,:)),abs(zRec{iTraj,2}(3,:)),'Linewidth',1.5,'displayname',['Traj. ',trajNames{iTraj}])
end
xlabel('$\rho_1$','interpreter','latex')
ylabel('$\rho_2$','interpreter','latex')
legend('interpreter','latex','NumColumns',3)
set(gca,'xscale','log')
set(gca,'yscale','log')
ylim([1e-5,1e1])
%% 
% We now plot the instantaneous damping on the decaying trajectories excited 
% from the third impact position.

customFigure;
for iTraj = 9:nTraj-1
    instDamping = real(Nflow(0,zRec{iTraj,2})./exp(1i*angle(zRec{iTraj,2})))./abs(zRec{iTraj,2});
    yyaxis left
    plot(zRec{iTraj,1},instDamping(1,:),'Linewidth',1.5,'displayname',['Mode 1, Traj. ',trajNames{iTraj}])
    yyaxis right
    plot(zRec{iTraj,1},instDamping(3,:),'Linewidth',1.5,'displayname',['Mode 2, Traj. ',trajNames{iTraj}])
end
xlabel('time [s]','interpreter','latex')
xlim([zRec{iTraj,1}(1), zRec{iTraj,1}(end)])
yyaxis left
ylabel('$c_1$ $\left[\frac{1}{\mathrm{s}}\right]$','interpreter','latex')
ylim([-1.5,2])
yyaxis right
ylabel('$c_2$ $\left[\frac{1}{\mathrm{s}}\right]$','interpreter','latex')
ylim([-30,40])
legend('interpreter','latex','NumColumns',2)
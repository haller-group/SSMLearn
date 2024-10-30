%% Identifying nonlinear resonances on a 4D SSM from decaying sloshing data
% See J. Axas, B. Bauerlein, K. Avila, and G. Haller. Data-driven modeling 
% of subharmonic forced response due to nonlinear resonance. Scientific 
% Reports, 2024.

clearvars
close all
clc

%% Example setup
nicered = [0.7,0.1,0.1];
niceblue = [0.1,0.1,0.7];
nicegreen = [0.1,0.9,0.1];
nicegray = [0.6,0.6,0.6];

%% Delay embedding
load sloshingResonanceData.mat
% Note: Data dimensionality has been scaled down 20 times to save repo 
% space. To retrieve the exact results from Axas et al (2024), download the
% full dataset from the ETH library (url pending)

nTraj = size(xData,1);
width = 500;
nSurfacePoints = size(xData{1,2}, 1);
dt = 0.01;
indDecay = 1:3;
indForced = 4;

shiftsteps = 5;
embedDim = 47;
outdof = 1; % leftmost point

yData = embedCoordinates(xData, embedDim, shiftsteps);

sliceInt = [1.2, Inf];
yDataTrunc = sliceTrajectories(yData, sliceInt);
clear yData;

%% Delay embedded tangent space (see Axas & Haller. Nonlinear dyn (2023))
modes = [1, 2];

SSMDim = 2*length(modes);
SSMOrder = 3;
ROMOrder = 5;
NFOrder = 7;

c = [-0.052,-0.07,-0.075,-0.08,-0.1]; % estimated real eigenvalue parts
frequencies = sqrt(9.81*pi/(width/1000)*(1:5).*tanh(pi*(1:5)*400/width)); % eigenfrequencies from potential theory
c = c(modes); frequencies = frequencies(modes);
lambda = reshape([c;c], 1, []) + 1i*reshape([frequencies;-frequencies], 1, []);

W = cos((1/2+linspace(-1/2, 1/2, nSurfacePoints)')*modes*pi); % sloshing mode shapes from theory

Vn = delayTangentSpace(shiftsteps*dt, embedDim, lambda); % Vandermonde matrix

clear V
for mode = modes
    V(:,2*mode-1) = reshape(W(:,mode)*Vn(:,2*mode-1).', [], 1);
    V(:,2*mode) = reshape(W(:,mode)*Vn(:,2*mode).', [], 1);
end
V = V./vecnorm(V); % SSM tangent space

%% fastSSM
[Mmap, iMmap, Tmap, iTmap, Nflow] = fastSSMplus(yDataTrunc(indDecay,:), SSMDim, SSMOrder, ROMOrder, NFOrder, 'mini', V, 0.05);

%% Transformation to normal form
clear zDataTrunc
for iTraj = 1:nTraj
    fprintf('transforming traj. %d out of %d to normal form...\n\n', iTraj, nTraj)
    zDataTrunc{iTraj,1} = yDataTrunc{iTraj,1};
    for iN = 1:size(yDataTrunc{iTraj,1},2)
        if mod(iN-1, 250) == 0
            fprintf('computed %5d points out of %5d... %5.1f %% done\n', iN-1, size(yDataTrunc{iTraj,1},2), 100*((iN-1)/size(yDataTrunc{iTraj,1},2)/nTraj+(iTraj-1)/nTraj))
        end
        zDataTrunc{iTraj,2}(:,iN) = iTmap(iMmap(yDataTrunc{iTraj,2}(:,iN)));
    end
    fprintf('computed %5d points out of %5d... %5.1f %% done\n\n', iN, size(yDataTrunc{iTraj,1},2), 100*iTraj/nTraj)
end

%% Reconstruction
zRec = integrateFlows(@(z)Nflow(0,z), zDataTrunc(indDecay,:));
yRec = transformTrajectories(@(z)Mmap(Tmap(z)), zRec);
NMTE = computeTrajectoryErrors(yRec, yDataTrunc(indDecay,:));

%% Plot nonlinear frequency relations
[X,Y] = meshgrid(1e-6:0.02:1,1e-6:0.02:0.62);
XY = 1e-6*ones(SSMDim,length(X(:)));
XY(1:2,:) = [X(:),X(:)]';
XY(3:4,:) = [Y(:),Y(:)]';
Nsurf = imag(Nflow(0,XY)./XY);
Nsurf1 = reshape(Nsurf(1,:)',size(X));
Nsurf2 = reshape(Nsurf(3,:)',size(X));
paperFigure('x', '$\rho_1$', 'y', '$\rho_2$', 'legendcols', 1);
for iTraj = indDecay
    rho = abs(zDataTrunc{iTraj,2});
    plot(rho(1,:),rho(3,:), 'linewidth', 4, 'DisplayName', ['Traj. ', num2str(iTraj)])
end
rho = abs(zDataTrunc{indForced,2});
plot(rho(1,:),rho(3,:), 'linewidth', 4, 'DisplayName', 'Forced data')
[c,h] = contour(X, Y, Nsurf2./Nsurf1, 'DisplayName', '$\dot{\theta}_2/\dot{\theta}_1$');
clabel(c,h)

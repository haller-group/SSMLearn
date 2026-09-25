clearvars
clc
%%
% The data loaded in this code is the velocity data of the lid-driven
% cavity simulation projected to the leading 10 POD modes for compression
% purposes.
% The same data projected to the leading 100 POD modes along with the
% reconstruction basis can be found under:
% https://polybox.ethz.ch/index.php/s/R5M4BJdpEyTfznB
% In order to run this code on the larger dataset, download the data files 
% and place them into the directory of this file. Then set
is_compressed_data = false; % to false, and run.
% The experiments reported in the paper were run on the full phase space
% data.

% In order to ensure a consistent basis ordering as the Reynolds number
% changes, this code uses munkres:
% "Munkres' Assignment Algorithm, Modified for Rectangular Matrices", 
% http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html
% version 2.3 by Yi Cao at Cranfield University on 11th September 2011


% James King, 17.03.2026

%% --- Setup ---

Re_train = [7900 7950 8200 8300 8400 8500]; % Parameters to train parametric model on
Re_target = 8450; % Parameter to evaluate parametric model on. Choose from [7925, 7975, 8100, 8350, 8450]
SSMDim = 2; % Dimension of the SSM
SSMOrder = 4; % Polynomial order of the SSM parameterization
ROMOrder = 5; % Polynomial order of the reduced dynamics
indTrain = 1; % Train on just the first trajectory

base_train_dir = 'IMInfoForTraining';

%% --- Training --

for i = 1:numel(Re_train)
    Re = Re_train(i);
    fprintf('Processing Re = %d\n', Re);
    
    % Load training data:
    if is_compressed_data
        data_file = fullfile(sprintf('Re%d_data_compressed.mat', Re));
        load(data_file, 'xData');
    else
        data_file = fullfile(sprintf('Re%d_data.mat', Re));
        load(data_file, 'xData', 'compression_basis', 'dof_coords', 'steady_state');
    end

    [yData, opts_embd] = coordinatesEmbedding(xData, SSMDim);

    dataMatrix = horzcat(yData{indTrain,2});
    [u, s, v] = svds(dataMatrix, SSMDim); % calculate leading 2 POD modes to graph SSM over
    
    % Ensure consistent ordering of u using munkres:
    if exist('u_ref', 'var')
        corrmat = abs(u_ref' * u);
        [assignment,~] = munkres(-corrmat);
        u = u(:,assignment);
        for j = 1:size(u,2)
            if dot(u(:,j), u_ref(:,j)) < 0
                u(:,j) = -u(:,j);
            end
        end
    end
    u_ref = u; % Update reference

    IMInfo = IMGeometry(yData(indTrain,:), SSMDim, SSMOrder, 'chart', @(x) u'*x); % Fit SSM geometry
    IMInfo.dual_basis = u'; % Add dual basis attribute for interpolation
    IMInfo.chart.map = @(x) IMInfo.dual_basis * x; % update chart for consistency
    etaData = projectTrajectories(IMInfo, yData); % Project to reduced coordinates (here called eta)

    RDInfo = IMDynamicsFlow(etaData,'R_PolyOrd', ROMOrder); % Fit SSM-reduced dynamics
    % Other reduced dynamics styles also possible, change in interpolation
    % accordingly.
    % modal:
    % RDInfo = IMDynamicsFlow(etaData, ...
    % 'R_PolyOrd', ROMOrder, ...
    % 'style', 'modal');
    % normalform:
    % RDInfo = IMDynamicsFlow(etaData, ...
    % 'R_PolyOrd', ROMOrder, ...
    % 'style', 'normalform');

    output_dir = fullfile(base_train_dir, sprintf('Re%d_coeffs_poly_%d', Re, ROMOrder));
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    save(fullfile(output_dir, 'IMInfo.mat'), 'IMInfo');
    save(fullfile(output_dir, 'RDInfo.mat'), 'RDInfo');
end

%% --- Interpolation ---
config.Re_values = Re_train;
config.folderPattern = sprintf('Re%%d_coeffs_poly_%d', ROMOrder);
config.interpMethod = 'spline'; % or 'linear'
config.style = 'polynomial'; % or 'modal' or 'normalform'
config.baseTrainDir = base_train_dir;
config.baseOutputDir = 'IMInfoRDInfoInterpolated';

createInterpolatedIMInfo_Re(Re_target, config); % Interpolate SSM geometry coefficients (tangent space, dual space, polynomial coefficients)
createInterpolatedRDInfo_Re(Re_target, config); % Interpolate SSM-reduced dynamics coefficients

interp_dir = fullfile(config.baseOutputDir, sprintf('Re%d_coeffs_poly_%d', Re_target, ROMOrder));
IMInfo = load(fullfile(interp_dir, 'IMInfo.mat'), 'IMInfo'); IMInfo = IMInfo.IMInfo;
RDInfo = load(fullfile(interp_dir, 'RDInfo.mat'), 'RDInfo'); RDInfo = RDInfo.RDInfo;

%% --- Testing ---
% Load testing data:
if is_compressed_data
    load(fullfile(sprintf('Re%d_data_compressed.mat', Re_target)), ...
        'xData');
else
    load(fullfile(sprintf('Re%d_data.mat', Re_target)), ...
        'xData', 'compression_basis', 'dof_coords', 'steady_state');
end

[yDataTest, opts_embd] = coordinatesEmbedding(xData, SSMDim);
etaDataTest = projectTrajectories(IMInfo, yDataTest); % project using interpolated chart
[yRecTest, etaRecTest, zRecTest] = advect(IMInfo, RDInfo, yDataTest); % advect and reconstruct using interpolated coefficients

%% --- Calculate Errors ---

[normedTrajDist, ampErrors] = computeTrajectoryErrors(yRecTest, yDataTest);

NMTE = mean(normedTrajDist)*100      % Normalized Mean Trajectory Error (%)
NMAE = mean(ampErrors)*100           % Normalized Mean Amplitude Error (%)

RDError = mean(computeTrajectoryErrors(etaRecTest, etaDataTest))*100 % NMTE of reduced dynamics (%)

yDataTestLifted = liftTrajectories(IMInfo, etaDataTest);
ReconstructionError = mean(computeTrajectoryErrors(yDataTestLifted, yDataTest))*100

%% --- Plotting ---

figure('Position', [100, 100, 900, 400]);

% --- Phase Portrait: eta_1 vs eta_2 ---
subplot(1,2,1);
hold on; grid on;
plotIdx = 1;
% Interpolate for smoother curves
x = real(etaDataTest{plotIdx,2}(1,:));
y = real(etaDataTest{plotIdx,2}(2,:));
t = 1:length(x);
t_fine = linspace(1, length(x), 10*length(x));
x_smooth = interp1(t, x, t_fine, 'spline');
y_smooth = interp1(t, y, t_fine, 'spline');
plot(x_smooth, y_smooth, 'k-', 'LineWidth', 2.0, 'DisplayName', 'Ground Truth');

x_pred = real(etaRecTest{plotIdx,2}(1,:));
y_pred = real(etaRecTest{plotIdx,2}(2,:));
t_pred = 1:length(x_pred);
t_pred_fine = linspace(1, length(x_pred), 10*length(x_pred));
x_pred_smooth = interp1(t_pred, x_pred, t_pred_fine, 'spline');
y_pred_smooth = interp1(t_pred, y_pred, t_pred_fine, 'spline');
plot(x_pred_smooth, y_pred_smooth, '-', 'Color', [214 101 2 255*0.8]/255, 'LineWidth', 1.5, 'DisplayName', 'Prediction');
xlabel('$u_1$', 'Interpreter', 'latex');
ylabel('$u_2$', 'Interpreter', 'latex');
legend('Location', 'northwest');
axis equal;

% --- Time Series: eta_1 vs Time ---
subplot(1,2,2);
hold on; grid on;
t = etaDataTest{plotIdx,1};
y = real(etaDataTest{plotIdx,2}(1,:));
t_fine = linspace(t(1), t(end), 10*length(t));
y_smooth = interp1(t, y, t_fine, 'spline');
plot(t_fine, y_smooth, 'k-', 'LineWidth', 2.0, 'DisplayName', 'Ground Truth');

y_pred = real(etaRecTest{plotIdx,2}(1,:));
y_pred_smooth = interp1(t, y_pred, t_fine, 'spline');
plot(t_fine, y_pred_smooth, '-', 'Color', [214 101 2 255*0.8]/255, 'LineWidth', 1.5, 'DisplayName', 'Prediction');
xlabel('Time [s]');
ylabel('$u_1$', 'Interpreter', 'latex');
legend('Location', 'northwest');

% --- 3D plot against shift mode ---
plot_shift_mode_surface(IMInfo, yDataTest, yRecTest, 1, 1)

%% --- Optional: Full Phase space plot ---
% % Reconstruct full phase space data
if ~is_compressed_data
    yData_full = yDataTest;
    yRec_full = yRecTest;
    for k = 1:size(yDataTest,1)
        yData_full{k,2} = compression_basis * yDataTest{k,2};
        yRec_full{k,2} = compression_basis * yRecTest{k,2};
    end
    % --- Full Phase Space Plots ---
    snapshotIdx = 2600;
    grid_resolution = 2000;
    plot_velocity_fields(real(yRec_full{1,2} + 0*steady_state), dof_coords, snapshotIdx, grid_resolution);
end
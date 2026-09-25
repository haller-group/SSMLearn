function plot_shift_mode_surface(IMInfo, yDataTest, yDataRec, idx, start_idx, show_predictions)
    % plot_shift_mode_surface visualizes a 2D SSM with trajectories
    % overlayed as a surface over the reduced coordinates against the shift
    % mode.
    %   Inputs:
    %     IMInfo           - SSM geometry and parametrization structure.
    %     yDataTest        - Cell array of test trajectories.
    %     yDataRec         - Cell array of predicted trajectories.
    %     idx              - Index of trajectory to plot.
    %     start_idx        - Starting time index (optional, default: 1).
    %     show_predictions - If true, also plots prediction (optional, default: true).
    
    if nargin < 6
        show_predictions = true;
    end
    if nargin < 5
        start_idx = 1;
    end
    yDataTest{idx,1} = yDataTest{idx,1}(:, start_idx:end);
    yDataTest{idx,2} = yDataTest{idx,2}(:, start_idx:end);
    yDataRec{idx,1} = yDataRec{idx,1}(:, start_idx:end);
    yDataRec{idx,2} = yDataRec{idx,2}(:, start_idx:end);
    
    % 1. Take the mean of the tail end of the trajectory
    tail_length = 800; % Number of time steps to average
    u_delta_a = mean(yDataTest{idx,2}(:, end-tail_length+1:end), 2); % [n_dof x 1]
    
    % 2. Project onto reduced basis (V) and reconstruct
    encoded = IMInfo.chart.map(u_delta_a);              % [SSMDim x 1]
    u_delta_a_proj = IMInfo.parametrization.map(encoded); % [n_dof x 1]
    
    % 3. Compute shift mode (difference, normalized)
    u_delta_b = u_delta_a - u_delta_a_proj;
    shift_mode = u_delta_b / norm(u_delta_b);
    
    % 4. Automatically set grid radius based on trajectory
    eta_traj = IMInfo.chart.map(yDataRec{idx,2}); % [2 x T]
    margin = 0.0001;
    max_radius = real(max(sqrt(sum(eta_traj.^2, 1))));
    radius = max_radius * (1 + margin);
    
    % 5. Create grid in eta space
    n_radial = 50;
    n_angular = 100;
    r = linspace(0, radius, n_radial);
    theta = linspace(0, 2*pi, n_angular);
    [R, Theta] = meshgrid(r, theta);
    eta1 = R(:) .* cos(Theta(:));
    eta2 = R(:) .* sin(Theta(:));
    eta_pts = [eta1'; eta2'];
    
    % 6. Decode grid points to full state space
    decoded = IMInfo.parametrization.map(eta_pts);  % [n_dof x N]
    
    % 7. Compute shift mode coefficient at each grid point
    shift_coeff_surface = shift_mode' * decoded; % [1 x N]
    shift_coeff_surface = reshape(shift_coeff_surface, size(R));
    
    % 8. Get shift mode coefficient along test and predicted trajectory
    eta_traj_test = IMInfo.chart.map(yDataTest{idx,2}); % [2 x T]
    shift_coeff_traj_test = shift_mode' * yDataTest{idx,2}; % [1 x T]
    
    % 9. Plot surface and trajectories
    figure;
    surf(R .* cos(Theta), R .* sin(Theta), real(shift_coeff_surface), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    plot3(eta_traj_test(1,:), real(eta_traj_test(2,:)), shift_coeff_traj_test, 'k-', 'LineWidth', 2);
    
    if show_predictions
        eta_traj_pred = IMInfo.chart.map(yDataRec{idx,2}); % [2 x T]
        shift_coeff_traj_pred = shift_mode' * yDataRec{idx,2}; % [1 x T]
        h_pred = plot3(eta_traj_pred(1,:), real(eta_traj_pred(2,:)), shift_coeff_traj_pred, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 2);
        legend({'SSM', 'Ground Truth', 'Prediction'}, 'Location', 'northeast');
    else
        legend({'SSM', 'Ground Truth'}, 'Location', 'northeast');
    end
    
    xlabel('$u_1$'); ylabel('$u_2$'); zlabel('Shift Mode Coefficient');
    grid on;
end
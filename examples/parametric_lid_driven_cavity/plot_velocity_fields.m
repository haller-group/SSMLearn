function plot_velocity_fields(U_data, dof_coords, snapshot_idx, grid_resolution, figsize)
% plot_velocity_fields visualizes 2D velocity fields from FEniCS data.
%   Plots the velocity magnitude, u, and v components as images over a regular grid.
%   Inputs:
%     U_data         - [n_dofs x n_snapshots] or [n_dofs x 1] (single snapshot)
%     dof_coords     - [n_dofs x 2] velocity DOF coordinates
%     snapshot_idx   - Snapshot to plot (optional, default: 1)
%     grid_resolution- Grid size for interpolation (optional, default: 1000)
%     figsize        - [width, height] in pixels (optional, default: [2400, 800])

if nargin < 3 || isempty(snapshot_idx)
    snapshot_idx = 1;
end
if nargin < 4 || isempty(grid_resolution)
    grid_resolution = 1000;
end
if nargin < 5 || isempty(figsize)
    figsize = [2400, 800];
end

% Extract snapshot
if ndims(U_data) == 3
    U_snapshot = squeeze(U_data(:, snapshot_idx, 1));
elseif ndims(U_data) == 2
    U_snapshot = U_data(:, snapshot_idx);
else
    U_snapshot = U_data;
end

% Extract u and v values (interleaved DOF format)
u_vals = U_snapshot(1:2:end);
v_vals = U_snapshot(2:2:end);

% Extract coordinates (every other for velocity DOFs)
coords = dof_coords(1:2:end, :);

% Calculate velocity magnitude
vel_mag = sqrt(u_vals.^2 + v_vals.^2);

% Create regular grid for interpolation
x_min = min(coords(:,1)); x_max = max(coords(:,1));
y_min = min(coords(:,2)); y_max = max(coords(:,2));
xi = linspace(x_min, x_max, grid_resolution);
yi = linspace(y_min, y_max, grid_resolution);
[XI, YI] = meshgrid(xi, yi);

% Interpolate data onto grid
vel_mag_grid = griddata(coords(:,1), coords(:,2), vel_mag, XI, YI, 'linear');
u_vals_grid  = griddata(coords(:,1), coords(:,2), u_vals,  XI, YI, 'linear');
v_vals_grid  = griddata(coords(:,1), coords(:,2), v_vals,  XI, YI, 'linear');

% Plot
figure('Position', [100, 100, figsize(1), figsize(2)]);
subplot(1,3,1);
imagesc(xi, yi, vel_mag_grid); axis xy equal tight;
title('Velocity Magnitude');
xlabel('x'); ylabel('y');
cb = colorbar;
ticks = cb.Ticks;
cb.TickLabels = arrayfun(@(x) sprintf('%.3f', x), ticks, 'UniformOutput', false);

subplot(1,3,2);
imagesc(xi, yi, u_vals_grid); axis xy equal tight;
title('u-velocity');
xlabel('x'); ylabel('y');
cb = colorbar;
ticks = cb.Ticks;
cb.TickLabels = arrayfun(@(x) sprintf('%.3f', x), ticks, 'UniformOutput', false);

subplot(1,3,3);
imagesc(xi, yi, v_vals_grid); axis xy equal tight;
title('v-velocity');
xlabel('x'); ylabel('y');
cb = colorbar;
ticks = cb.Ticks;
cb.TickLabels = arrayfun(@(x) sprintf('%.3f', x), ticks, 'UniformOutput', false);

sgtitle(sprintf('Velocity fields (Snapshot %d)', snapshot_idx));
end

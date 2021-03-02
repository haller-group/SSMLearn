% Test Code for data driven invarian manifolds
%
% Reconstruction of slow 2D SSM of a 3 Dof mechanical system observing
% a scalar measurement. Do not change the variables in the first section
% otherwise the SSM initial condition in the .mat file won't be valid
% anymore

clear all
close all
format shortg
clc

%% Input

% Define equations as:
% M qdd + C qd + K q + fq(q,qd) = 0
% Number of mechanical degrees of freedom
N = 3;
% Mass Matrix
M = eye(N);
% Stiffness Matrix
K = diag(2*ones(1,N))-diag(ones(1,N-1),-1)-diag(ones(1,N-1),+1);
% Damping Matrix
Kc = K; Kc(1,1)=4; Kc(end,end)=8;
C = 0.006*Kc; C(end,end) = 3*0.006;
% Nonlinearities
f_q1 = @(x) 0.33*x(4).^2+2*x(1).^3+0.3*x(1).^2.*x(4)+0.5*x(4).^3;
f_q1_eval = @(x) 0.33*x(4,:).^2+2*x(1,:).^3+0.3*x(1,:).^2.*x(4,:)+0.5*x(4,:).^3;

% Linear Dynamics
A = [zeros(N) eye(N); -M\K -M\C];
[V,D] = eig(A);
[~,pos] = sort(real( diag (D) ),'descend');
V = V(:,pos); l = diag(D(pos,pos))
% Define non complex modal coordinates
V_real = [];
for ii=1:N
    V_real = [V_real real(V(:,2*ii-1)) imag(V(:,2*ii-1))];
end
V_ort = V_real(:,1)/norm(V_real(:,1));
for ii = 2:2*N
    vadd = V_real(:,ii) - V_ort * transpose(V_ort) * V_real(:,ii);
    vadd = vadd/norm(vadd);
    V_ort = [V_ort vadd];
end

% Vector field of equivalent first order system
f = @(t,x) A*x - [zeros(N,1); f_q1(x); zeros(N-1,1)];
f_eval = @(x) A*x - [zeros(N,size(x,2)); f_q1_eval(x); zeros(N-1,size(x,2))];

% Load Initial Conditions on the slowest SSM
load SlowSSM_InitialConditions.mat

%% Construction of Synthetic Data

% Set time instances and sampling time
Tend = 150; Nsamp = .5e4; Tscale = 2*pi/abs(imag(l(1)));
t_dd = linspace(0,Tend*Tscale,Nsamp); t_samp = t_dd(2);
N_traj = size(SSM_IC,2);
t_int = [0 Tend*Tscale];
X_sim = cell(N_traj,2);

% Integrate and save trajectory data
figure(1); clf; hold on; grid on; box on;
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
colors = winter(size(SSM_IC,2) ); % jet, spring, autumn, winter, summer cool
cc = 0;
for ii = 1:N_traj
    [t,x] = ode45(f,t_int,SSM_IC(:,ii),options);
        X_sim{ii,1} = t_dd;
        X_sim{ii,2} = transpose(interp1(t,x,t_dd,'spline'));
    x = transpose(x);
    q1 = x(1,:); q2 = x(4,:); q3 = x(6,:);
    plot3(q1,q2,q3,'Linewidth',.1,'Color',colors(ii,:))
end
xlabel('$q_1$','Interpreter','latex')
ylabel('$\dot{q}_1$','Interpreter','latex')
zlabel('$\dot{q}_3$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)


%% Datadriven SSM
u_v_a = 3; % 1 for displacement, for 2 velocity, 3 or acceleration
dof = 1; 
coord_to_embd = (u_v_a-1)*N+dof; % Acceleration of the first mass
X_obs_traj = cell(N_traj,2);

for ii = 1:N_traj
    x_i = [X_sim{ii,2}];
    xd_i = f_eval([X_sim{ii,2}]);
    x_ext_i = [x_i; xd_i(N+1:end,:)]; % Rows are disapl., vel. and acceleration
    X_obs_traj{ii,1} = t_dd;
    X_obs_traj{ii,2} = x_ext_i(coord_to_embd,:);
end

%%
% Data embedding
SSM_dim = 2; over_embd = 50 ;
[X_embd_traj,opts_embd] = coordinates_embedding(X_obs_traj,SSM_dim,'OverEmbedding',over_embd);
opts_embd

% Split in training and testing dataset
ind_test = [1 5];
ind_train = setdiff(1:N_traj, ind_test);
X_train = X_embd_traj(ind_train,:);
X_test = X_embd_traj(ind_test,:);

%%
% Polynomial degree of the parametrization
ParamDeg = 5; 
% Eigenspace
V = [];
% Optimization
[V_ort_data,SSM_func] = IMparametrization(X_train,SSM_dim,ParamDeg);

% Plot and validation

coordplot = 1;

figure(2); clf; hold on; grid on; box on;
figure(3); clf; hold on; grid on; box on;
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
colors = winter(size(SSM_IC,2) );
Err_rec = []; rho_max = 0; Y = [];

for ii = 1:length(ind_test)
    % Simulate synthetic data, get observables and reconstruct test data
    t = X_test{ii,1}; x = X_test{ii,2};
    y = transpose(V_ort_data)*x; Y = [Y y]; x_rec = SSM_func(y);
    rho_max = max([ rho_max max(sqrt(sum(y.^2)))]);
    % Error evaluation
    Err_rec = [Err_rec mean(sqrt(sum( (x-x_rec).^2 )))];
    if ii ==1
       figure(7); clf; hold on; grid on; box on;
       plot(t,x(2,:),'k','Linewidth',2)
       plot(t,x_rec(2,:),'r:','Linewidth',2)
       xlabel('$t$','Interpreter','latex')
       ylabel('$x$','Interpreter','latex')
       set(gca,'fontname','times')
       set(gca,'fontsize',18)
    end
    figure(2)
    q1 = y(1,:); q2 = y(2,:); q3 = x(coordplot,:);
    plot3(q1,q2,q3,'Linewidth',2,'Color',colors(ind_test(ii),:))
    figure(3)
    q1 = y(1,:); q2 = y(2,:);
    plot(q1,q2,'Linewidth',2,'Color',colors(ind_test(ii),:))
end
xlabel('$\eta_1$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
ARSE = mean(Err_rec)


% Plot 3D SSM Surface
figure(2);
h = plot_2dSSM_surf(coordplot,Y,SSM_func,10,50,0);
xlabel('$\eta_2$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
zlabel('$x$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)

%% Reduced Dynamics


% Arrange trajectories
Y_train = cell(N_traj-length(ind_test),2);
for ii = 1:size(X_train,1)
    Y_train{ii,1} = X_train{ii,1};
    Y_train{ii,2} = transpose(V_ort_data)*X_train{ii,2};
end
Y_test  = cell(length(ind_test),2);
for ii = 1:size(X_test,1)
    Y_test{ii,1} = X_test{ii,1};
    Y_test{ii,2} = transpose(V_ort_data)*X_test{ii,2};
end

% Dynamics identification
PolyOrd = 7;
[R,~,~,~,Maps_info] = IMdynamics_map(Y_train,PolyOrd);
R_info = Maps_info.R

%% Error of the dynamics

% Error on training trajectory
ii = 1;
t_i = Y_train{ii,1};
y_i = Y_train{ii,2};
y_sim = iterate_map(R,length(t_i),y_i(:,1));
figure(41); clf;
subplot(211); hold on; grid on; box on;
plot(t_i,y_i(1,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,y_sim(1,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
legend
subplot(212); hold on; grid on; box on;
plot(y_i(1,:),y_i(2,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(y_sim(1,:),y_sim(2,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
RMSE = mean(sqrt(sum( (y_i-y_sim).^2 )))

% Error on testing trajectory
t_i = Y_test{ii,1};
y_i = Y_test{ii,2};
y_sim = iterate_map(R,length(t_i),y_i(:,1));
figure(42); clf;
subplot(211); hold on; grid on; box on;
plot(t_i,y_i(1,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,y_sim(1,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
legend
subplot(212); hold on; grid on; box on;
plot(y_i(1,:),y_i(2,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(y_sim(1,:),y_sim(2,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$\eta_1$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
RMSE = mean(sqrt(sum( (y_i-y_sim).^2 )))

%% Overall ROM performance

coordplot = 1; i_test = 1;
% Error on testing trajectory
t_i = X_test{ii,1};
x_i = X_test{ii,2};
y_i = transpose(V_ort_data)*x_i;
y_ROM = iterate_map(R,length(t_i),y_i(:,1));
x_ROM = SSM_func(y_ROM);
figure(42); clf;
subplot(211); hold on; grid on; box on;
plot(t_i,x_i(coordplot,:),'k','Linewidth',2,'DisplayName','Testing Trajectory')
plot(t_i,x_ROM(coordplot,:),'r:','Linewidth',2,'DisplayName','ROM Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
legend
subplot(212); hold on; grid on; box on;
plot(y_i(1,:),y_i(2,:),'k','Linewidth',2,'DisplayName','Testing Trajectory')
plot(y_ROM(1,:),y_ROM(2,:),'r:','Linewidth',2,'DisplayName','ROM Trajectory')
xlabel('$\eta_1$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
RRMSE = mean(sqrt(sum( (x_i-x_ROM).^2 )))/max(sqrt(sum(x_i.^2)))*100 % Percentage error based on the max ||x||

%% Sanity check on eigenvalues
R_coeff = R_info.coeff;
R_coeff_lin = R_coeff(:,1:SSM_dim);
l_est = eig(R_coeff_lin);
[log(l_est)/(t_i(2)-t_i(1)) l(1:SSM_dim)]
% Analyze Transient 4D dataset

clear all
close all
clc

% Load dataset
load('SSM4D.mat')
% Plot Data
plot(X_traj{1,1},X_traj{1,2},X_traj{2,1},X_traj{2,2})
xlabel('$t \, [$s$]$','Interpreter','latex')
ylabel('$v \, [$m/s$]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
grid on

% Data Embedding (takens)
SSM_dim = 4; over_embd = 491;%191;
[X_traj,opts_embd] = coordinates_embedding(X_traj,SSM_dim,'OverEmbedding',over_embd,'TimeStepping',1);
opts_embd

% Split in training and testing dataset
ind_test = [];
X_train = X_traj;
X_test = X_traj(ind_test,:);
X_train(ind_test,:)=[];
%%
% Manifold learning
% Polynomial degree of the parametrization
ParamDeg = 2; 
% Optimization
[V_ort_data,SSM_func,IM_para_info] = IMparametrization(X_traj,SSM_dim,ParamDeg);

% Plot SSM Data

figure(3); clf; hold on; grid on; box on; colororder(winter(2))
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
Err_rec = []; rho_max = 0;
coordplot = [1];
for ii = 1:2%:length(ind_test)
    % Get reduced coordinates of testing data
    x = X_train{ii,2}; t = X_train{ii,1};
    y = transpose(V_ort_data)*x; 
    % Reconstruction of testing Trajectory
    x_rec = SSM_func(y);
    rho_max = max([ rho_max max(sqrt(sum(y.^2)))]);
    % Get Error
    Err_rec = [Err_rec mean(sqrt(mean( (x-x_rec).^2 )))];
    if ii==1
        % Plot error vs time
        figure(7); clf; hold on; grid on; box on;
        plot(t,x(coordplot(1),:),'k','Linewidth',2,'DisplayName','Full Trajectory')
        plot(t,x_rec(coordplot(1),:),'r:','Linewidth',2,'DisplayName','Reconstructed Trajectory')
        xlabel('$ t$','Interpreter','latex')
        ylabel('$x_1$','Interpreter','latex')
        set(gca,'fontname','times')
        set(gca,'fontsize',18)
        xlim([t(1) t(end)])
        legend
    end
    figure(3)
    q1 = y(1,:); q2 = y(2,:); q3 = y(3,:);
    plot3(q1,q2,q3,'Linewidth',2)
end
figure(3)
xlabel('$\eta_1$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
zlabel('$\eta_3$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
ARSE = mean(Err_rec)
view(3)
drawnow;
%%
% Arrange trajectories
N_traj = size(X_traj,1);
Y_train = cell(N_traj-length(ind_test),2);
for ii = 1:size(X_train,1)
    Y_train{ii,1} = [X_train{ii,1}];
    Y_train{ii,2} = transpose(V_ort_data)*[X_train{ii,2}];
end
Y_test  = cell(length(ind_test),2);
for ii = 1:size(X_test,1)
    Y_test{ii,1} = [X_test{ii,1}];
    Y_test{ii,2} = transpose(V_ort_data)*[X_test{ii,2}];
end

% Dynamics identification
PolyOrd = 4; 
[R,iT,Nf,T,Maps_info] = IMdynamics_flow(Y_train,'R_PolyOrd',PolyOrd,'style','normalform','IC_nf',1);
R_info = Maps_info.R;
N_info = Maps_info.N;
T_info = Maps_info.T;
iT_info = Maps_info.iT;
%%
% Error on training trajectories
ii = 1;
t_i = Y_train{ii,1};
y_i = Y_train{ii,2};
[t_sim,y_sim] = ode45(@(t,x) R(x),t_i,y_i(:,1)); y_sim = transpose(y_sim);
figure(41); clf;
subplot(211); hold on; grid on; box on;
plot(t_i,y_i(1,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_sim,y_sim(1,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
ylim([-1 1]*max(abs(y_i(1,:))))
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(t_i,y_i(3,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_sim,y_sim(3,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_3$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
RMSE = mean(sqrt(sum( (y_i-y_sim).^2 )))

ii = 2;
t_i = Y_train{ii,1};
y_i = Y_train{ii,2};
[~,y_sim] = ode45(@(t,x) R(x),t_i,y_i(:,1)); y_sim = transpose(y_sim);
figure(42); clf;
subplot(211); hold on; grid on; box on;
plot(t_i,y_i(1,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,y_sim(1,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
ylim([-1 1]*max(abs(y_i(1,:))))
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
subplot(212); hold on; grid on; box on;
plot(t_i,y_i(3,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,y_sim(3,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_3$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
RMSE = mean(sqrt(sum( (y_i-y_sim).^2 )))

% Linearized System
ii = 1;
t_i = Y_train{ii,1};
y_i = Y_train{ii,2};
W_r = R_info.coeff;
R_lin= @(x) -W_r(:,1:size(W_r,1))*x;
[~,y_sim] = ode45(@(t,x) R_lin(x),t_i,y_i(:,end)); y_sim =fliplr( transpose(y_sim));
figure(43); clf;
subplot(211); hold on; grid on; box on;
plot(t_i,y_i(1,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,y_sim(1,:),'c:','Linewidth',2,'DisplayName','Simulated Trajectory (linearized dynamics)')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
ylim([-1 1]*max(abs(y_i(1,:))))
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(t_i,y_i(3,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,y_sim(3,:),'c:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_3$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
RMSE = mean(sqrt(sum( (y_i-y_sim).^2 )))

%% Normal Form
% Error on training trajectory
ii = 1;
t_i = Y_train{ii,1}; 
y_i = iT(Y_train{ii,2}); 
[t_sim,y_sim] = ode45(@(t,x) Nf(x),t_i,y_i(:,1)); y_sim = transpose(y_sim);
figure(51); clf;
subplot(211); hold on; grid on; box on;
plot(t_i,abs(y_i(1,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,abs(y_sim(1,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\rho_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(t_i,abs(y_i(2,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,abs(y_sim(2,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\rho_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
RMSE = mean(sqrt(sum( (y_i-y_sim).*conj(y_i-y_sim) )))

figure(52); clf;
subplot(211); hold on; grid on; box on;
plot(t_i,real(y_i(1,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,real(y_sim(1,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\mathrm{Re}(z_1)$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
legend
xlim([t_i(1) t_i(end)])
subplot(212); hold on; grid on; box on;
plot(t_i,real(y_i(2,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,real(y_sim(2,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\mathrm{Re}(z_2)$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])

%% Overall ROM performance - default

coordplot = 1; ii = 1;
% Error on testing trajectory
t_i = X_train{ii,1}; 
x_i = X_train{ii,2}; 
y_i = (transpose(V_ort_data)*x_i);
[~,y_ROM] = ode45(@(t,x) R(x),t_i,y_i(:,1)); y_ROM = transpose(y_ROM);
x_ROM = SSM_func((y_ROM));
figure(91); clf; hold on; grid on; box on; 
plot(t_i,x_i(coordplot,:),'k','Linewidth',2,'DisplayName','Trajectory 1')
plot(t_i,x_ROM(coordplot,:),'r:','Linewidth',2,'DisplayName','ROM Trajectory')
xlabel('$t \, [$s$]$','Interpreter','latex')
ylabel('$v \, [$m/s$]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
RRMSE = mean(sqrt( (x_i(coordplot,:)-x_ROM(coordplot,:)).^2 ))/max(sqrt(x_i(coordplot,:).^2))*100 % Percentage error based on the max ||x||

% Overall ROM performance - normal form
coordplot = 1; ii = 1;
% Error on testing trajectory
t_i = X_train{ii,1}; 
x_i = X_train{ii,2}; 
y_i = iT(transpose(V_ort_data)*x_i);
[t_sim,y_ROM] = ode45(@(t,x) Nf(x),t_i,y_i(:,1)); y_ROM = transpose(y_ROM);
x_ROM = SSM_func(T(y_ROM));
figure(92); clf; hold on; grid on; box on; %title('Overall ROM performance - flow')
plot(t_i,x_i(coordplot,:),'k','Linewidth',2,'DisplayName','Trajectory 1')
plot(t_i,x_ROM(coordplot,:),'r:','Linewidth',2,'DisplayName','ROM Trajectory')
xlabel('$t \, [$s$]$','Interpreter','latex')
ylabel('$v \, [$m/s$]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
RRMSE = mean(sqrt( (x_i(coordplot,:)-x_ROM(coordplot,:)).^2 ))/max(sqrt(x_i(coordplot,:).^2))*100 % Percentage error based on the max ||x||

%% Frequencies and damping rations each trajectory

[damp,freq] = nonres_normalform(N_info.coeff,N_info.exponents);
colors_m = [0 0.45 0.74; 0.85 0.33 0.1];
ii = 1;
t_1 = Y_train{ii,1}; 
y_i = iT(Y_train{ii,2}); 
[~,y_sim] = ode45(@(t,x) Nf(x),t_1,y_i(:,1)); y_sim = transpose(y_sim);
r_1 = abs(y_sim(1:2,:));
insta_damp_1 = damp(r_1); insta_freq_1 = freq(r_1);
insta_dr_1 = -insta_damp_1./insta_freq_1;
ii = 2;
t_2 = Y_train{ii,1}; 
y_i = iT(Y_train{ii,2}); 
[~,y_sim] = ode45(@(t,x) Nf(x),t_2,y_i(:,1)); y_sim = transpose(y_sim);
r_2 = abs(y_sim(1:2,:));
insta_damp_2 = damp(r_2); insta_freq_2 = freq(r_2);
insta_dr_2 = -insta_damp_2./insta_freq_2;
figure(101); clf;
subplot(211); hold on; grid on; box on;
yyaxis left
ylabel('$f_1 \, [$Hz$]$','Interpreter','latex')
plot(t_1,insta_freq_1(1,:)/2/pi,'Linewidth',2,'DisplayName','Mode 1 Traj 1')
plot(t_2,insta_freq_2(1,:)/2/pi,'Linewidth',2,'DisplayName','Mode 1 Traj 2')
yyaxis right
plot(t_1,insta_freq_1(2,:)/2/pi,'Linewidth',2,'DisplayName','Mode 2 Traj 1')
plot(t_2,insta_freq_2(2,:)/2/pi,'Linewidth',2,'DisplayName','Mode 2 Traj 2')
xlabel('$t \, [s]$','Interpreter','latex')
ylabel('$f_2 \, [$Hz$]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
subplot(212); hold on; grid on; box on;
yyaxis left
ylabel('$\zeta_1 \, [$-$]$','Interpreter','latex')
plot(t_1,insta_dr_1(1,:),'Linewidth',2,'DisplayName','Mode 1 Traj 1')
plot(t_2,insta_dr_2(1,:),'Linewidth',2,'DisplayName','Mode 1 Traj 2')
yyaxis right
plot(t_1,insta_dr_1(2,:),'Linewidth',2,'DisplayName','Mode 2 Traj 1')
plot(t_2,insta_dr_2(2,:),'Linewidth',2,'DisplayName','Mode 2 Traj 2')
xlabel('$t \, [s]$','Interpreter','latex')
ylabel('$\zeta_2 \, [$-$]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend


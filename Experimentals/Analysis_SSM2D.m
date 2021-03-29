clear all
close all
format shortg
clc

Resc = 1;
load('SSM2D.mat')
N_traj = size(X_traj,1);
for ii = 1:N_traj
    X_traj{ii,2} = X_traj{ii,2}/Resc;
end
figure(1); clf; hold on; grid on;  box on;
ii = 1; plot(X_traj{ii,1},X_traj{ii,2}*Resc,'Linewidth',1,'DisplayName','Trajectory 1')
ii = 2; plot(X_traj{ii,1},X_traj{ii,2}*Resc,'Linewidth',1,'DisplayName','Trajectory 2')
xlabel('$t \, [$s$]$','Interpreter','latex')
ylabel('$A \, [$m/s$^2]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
%%

% Data Embedding (takens)
SSM_dim = 2; over_embd = 195;
[X_traj,opts_embd] = coordinates_embedding(X_traj,SSM_dim,'OverEmbedding',over_embd,'TimeStepping',1);
opts_embd

% Split in training and testing dataset
ind_test = [];
X_train = X_traj;
X_test = X_traj(ind_test,:);
X_train(ind_test,:)=[];
%
% Manifold learning
% Polynomial degree of the parametrization
ParamDeg = 1;
% Optimization
[V_ort_data,SSM_func,IM_para_info] = IMparametrization(X_traj,SSM_dim,ParamDeg);
%%
% Plot SSM Data
figure(2); clf; hold on; grid on; box on; colororder(winter(2))
figure(3); clf; hold on; grid on; box on; colororder(winter(2))
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
Err_rec = []; Y = [];
coordplot = 1;
for ii = 1:2%:length(ind_test)
    % Get reduced coordinates of testing data
    x = X_train{ii,2}; t = X_train{ii,1};
    y = transpose(V_ort_data)*x; Y = [Y y];
    % Reconstruction of testing Trajectory
    x_rec = SSM_func(y);
    % Get Error
    Err_rec = [Err_rec mean(sqrt(mean( (x-x_rec).^2 )))];
    if ii==2
        % Plot error vs time
        figure(7); clf; hold on; grid on; box on;
        plot(t,x(coordplot(1),:),'k','Linewidth',2,'DisplayName','Full Trajectory')
        plot(t,x_rec(coordplot(1),:),'r:','Linewidth',2,'DisplayName','Reconstructed Trajectory')
        xlabel('$ t$','Interpreter','latex')
        ylabel('$q_1$','Interpreter','latex')
        set(gca,'fontname','times')
        set(gca,'fontsize',18)
        xlim([t(1) t(end)])
        legend
        figure(2)
    q1 = y(1,:); q2 = y(2,:); q3 = x_rec(coordplot,:);
    plot3(q1,q2,q3,'Linewidth',2)
    end
    figure(3)
    q1 = y(1,:); q2 = y(2,:);
    plot(q1,q2,'Linewidth',2)
end
figure(3)
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
PolyOrd = 8; Nfolds = 5; Nregvals = 30; 
[R,iT,N,T,Maps_info] = IMdynamics_flow(Y_train,'R_PolyOrd',PolyOrd,'n_folds',Nfolds,'l_vals',logspace(-6,0,Nregvals),'style','normalform');
R_info = Maps_info.R
%%
% Error on training trajectory
ii = 2;
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
xlim([t_i(1) t_i(end)])
ylim([-1 1]*max(abs(y_i(1,:))))
legend
subplot(212); hold on; grid on; box on;
plot(t_i,y_i(2,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_sim,y_sim(2,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
ylim([-1 1]*max(abs(y_i(2,:))))
RMSE = mean(sqrt(sum( (y_i-y_sim).^2 )))

coordplot = 1;
x_i = X_traj{ii,2};
x_ROM = SSM_func((y_sim));
figure(91); clf; hold on; grid on; box on; 
plot(t_i,x_i(coordplot,:),'k','Linewidth',2,'DisplayName','Trajectory 1')
plot(t_sim,x_ROM(coordplot,:),'r:','Linewidth',2,'DisplayName','ROM Trajectory')
xlabel('$t \, [$s$]$','Interpreter','latex')
ylabel('$a \, [$m/s$^2]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
RRMSE = mean(sqrt( (x_i(coordplot,:)-x_ROM(coordplot,:)).^2 ))/max(sqrt(x_i(coordplot,:).^2))*100
%% Normal form

% Modal
% Error on training trajectory
t_i = Y_train{ii,1};
y_i = iT(Y_train{ii,2});
[t_sim,y_sim] = ode45(@(t,x) N(x),t_i,y_i(:,1)); y_sim = transpose(y_sim);
figure(81); clf;
subplot(211); hold on; grid on; box on;
title('Reduced dynamics (flow): normal form - training trajectory')
plot(t_i,abs(y_i(1,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_sim,abs(y_sim(1,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\rho$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(real(y_i(1,:)),imag(y_i(1,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(real(y_sim(1,:)),imag(y_sim(1,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$z_1$','Interpreter','latex')
ylabel('$z_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
RMSE = mean(sqrt(sum( (y_i-y_sim).*conj(y_i-y_sim) )))
xlim([-1 1]*max(abs(y_i(1,:))))
ylim([-1 1]*max(abs(y_i(2,:))))

% Overall ROM performance
coordplot = 1; 
x_ROM = SSM_func(T(y_sim));
figure(92); clf;
subplot(211); hold on; grid on; box on; title('Overall ROM performance - flow')
plot(t_i,x_i(coordplot,:),'k','Linewidth',2,'DisplayName','Testing Trajectory')
plot(t_sim,x_ROM(coordplot,:),'r:','Linewidth',2,'DisplayName','ROM Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$x_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(real(y_i(1,:)),imag(y_i(2,:)),'k','Linewidth',2,'DisplayName','Testing Trajectory')
plot(real(y_sim(1,:)),imag(y_sim(2,:)),'r:','Linewidth',2,'DisplayName','ROM Trajectory')
xlabel('$z_1$','Interpreter','latex')
ylabel('$z_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([-1 1]*max(abs(real(y_i(1,:)))))
ylim([-1 1]*max(abs(imag(y_i(1,:)))))
RRMSE = mean(sqrt(sum( (x_i-x_ROM).^2 )))/max(sqrt(sum(x_i.^2)))*100 % Percentage error based on the max ||x||

%% Normal form and backbone curves
N_info = Maps_info.N; T =  Maps_info.T; T =  T.Map; 
[damp,freq] = nonres_normalform(N_info.coeff,N_info.exponents);
figure(101); clf; hold on; grid on; box on;
backbonecurves(damp,freq,SSM_func,T,coordplot,abs(y_i(1,1)),'norm')
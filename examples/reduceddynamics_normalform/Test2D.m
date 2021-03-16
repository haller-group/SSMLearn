% Test Code for data driven invarian manifolds
%
% Reconstruction of slow 2D SSM of a 5 Dof mechanical system observing
% the full state space. Do not change the variables in the first section
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
ndof = 5;
% Mass Matrix
M = eye(ndof);
% Stiffness Matrix
K = diag(2*ones(1,ndof))-diag(ones(1,ndof-1),-1)-diag(ones(1,ndof-1),+1);
% Damping Matrix
Kc = K; Kc(1,1)=4; Kc(end,end)=0;
C = 0.009*Kc; % Purely proportional damping
% Nonlinearities
f_q1 = @(x) 0.66*x(ndof+1).^2+3*x(1).^3+0.75*x(1).^2.*x(ndof+1)+1*x(ndof+1).^3;

% Linear Dynamics
A = [zeros(ndof) eye(ndof); -M\K -M\C];
[V,D] = eig(A);
[~,pos] = sort(real( diag (D) ),'descend');
V = V(:,pos); l = diag(D(pos,pos))
% Define non complex modal coordinates
V_real = [];
for ii=1:ndof
    V_real = [V_real real(V(:,2*ii-1)) imag(V(:,2*ii-1))];
end
V_ort = V_real(:,1)/norm(V_real(:,1));
for ii = 2:2*ndof
    vadd = V_real(:,ii) - V_ort * transpose(V_ort) * V_real(:,ii);
    vadd = vadd/norm(vadd);
    V_ort = [V_ort vadd];
end

% Vector field of equivalent first order system
f = @(t,x) A*x - [zeros(ndof,1); f_q1(x); zeros(ndof-1,1)];
% Load Initial Conditions on the slowest SSM
load Slow2DSSM_InitialConditions.mat

%% Construction of Synthetic Data

% Set time instances and sampling time
Tend = 150; Nsamp = 1e4; Tscale = 2*pi/abs(imag(l(1)));
t_dd = linspace(0,Tend*Tscale,Nsamp); t_samp = t_dd(2);
N_traj = size(SSM_IC,2);
t_int = [0 Tend*Tscale];
% Select eventual testing trajectories
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
    q1 = x(1,:); q2 = x(ndof+1,:); q3 = x(2*ndof,:);
    plot3(q1,q2,q3,'Linewidth',.1,'Color',colors(ii,:))
end
xlabel('$q_1$','Interpreter','latex')
ylabel('$\dot{q}_1$','Interpreter','latex')
zlabel('$\dot{q}_5$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
view(3)

%% Datadriven SSM

% Data embedding
SSM_dim = 2;
[X_traj,opts_embd] = coordinates_embedding(X_sim,SSM_dim);
opts_embd

% Split in training and testing dataset
ind_test = [1 5];
X_train = X_traj;
X_test = X_traj(ind_test,:);
X_train(ind_test,:)=[];

%%
% Polynomial degree of the parametrization
ParamDeg = 7; 
% Eigenspace
V = V_real(:,1:2);
% Optimization
[V_ort_data,SSM_func,IM_para_info] = IMparametrization(X_traj,SSM_dim,ParamDeg,V);

% Plot and validation
%sqrt(sum((V_ort_data'*V_ort).^2))
figure(2); clf; hold on; grid on; box on;
figure(3); clf; hold on; grid on; box on;
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
colors = winter(size(SSM_IC,2) );
Err_rec = []; rho_max = 0;
coordplot = [1 ndof+1 ndof*2];
for ii = 1:length(ind_test)
    % Get reduced coordinates of testing data
    x = X_test{ii,2}; t = X_test{ii,1};
    y = transpose(V_ort_data)*x; 
    % Reconstruction of testing Trajectory
    x_rec = SSM_func(y);
    rho_max = max([ rho_max max(sqrt(sum(y.^2)))]);
    % Get Error
    Err_rec = [Err_rec mean(sqrt(sum( (x-x_rec).^2 )))];
    if ii==1
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
    end
    figure(2)
    q1 = x(coordplot(1),:); q2 = x(coordplot(2),:); q3 = x(coordplot(3),:);
    plot3(q1,q2,q3,'Linewidth',.1,'Color',colors(ii,:))
    figure(3)
    q1 = y(1,:); q2 = y(2,:);
    plot(q1,q2,'Linewidth',2,'Color',colors(ii,:))
end
xlabel('$\eta_1$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
ARSE = mean(Err_rec)

% Plot 3D SSM Surface
figure(2);
h = plot_2dSSM_surf(coordplot,y,SSM_func,25,50,[0.39 0.83 0.07]);
xlabel('$q_1$','Interpreter','latex')
ylabel('$\dot{q}_1$','Interpreter','latex')
zlabel('$\dot{q}_3$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
%% Reduced Dynamics

% Arrange trajectories
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
PolyOrd = 10; 
[R,iT,N,T,Maps_info] = IMdynamics_map(Y_train,'R_PolyOrd',PolyOrd','style','normalform');
%% Error of the dynamics

% Error on training trajectory
t_i = Y_train{ii,1};
y_i = Y_train{ii,2};
y_sim = iterate_map(R,length(t_i),y_i(:,1));
figure(41); clf;
subplot(211); hold on; grid on; box on;
title('Reduced dynamics (map): default - training trajectory')
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
xlim([-1 1]*max(abs(y_i(1,:))))
ylim([-1 1]*max(abs(y_i(2,:))))
RMSE = mean(sqrt(sum( (y_i-y_sim).^2 )))

% Error on testing trajectory
t_i = Y_test{ii,1};
y_i = Y_test{ii,2};
y_sim = iterate_map(R,length(t_i),y_i(:,1));
figure(42); clf;
subplot(211); hold on; grid on; box on;
title('Reduced dynamics (map): default - testing trajectory')
plot(t_i,y_i(1,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,y_sim(1,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(y_i(1,:),y_i(2,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(y_sim(1,:),y_sim(2,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$\eta_1$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([-1 1]*max(abs(y_i(1,:))))
ylim([-1 1]*max(abs(y_i(2,:))))
RMSE = mean(sqrt(sum( (y_i-y_sim).^2 )))

%% Normal form
% Error on training trajectory
t_i = Y_train{ii,1}; 
y_i = iT(Y_train{ii,2}); 
t_i=t_i(4:end); y_i=y_i(:,4:end);
y_sim = iterate_map(N,length(t_i),y_i(:,1)); 
figure(51); clf;
subplot(211); hold on; grid on; box on;
title('Reduced dynamics (map): normal form - training trajectory')
plot(t_i,real(y_i(1,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,real(y_sim(1,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\mathrm{Re}(z_1)$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(real(y_i(1,:)),imag(y_sim(1,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(real(y_sim(1,:)),imag(y_sim(1,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$\mathrm{Re}(z_1)$','Interpreter','latex')
ylabel('$\mathrm{Im}(z_1)$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([-1 1]*max(abs(y_i(1,:))))
ylim([-1 1]*max(abs(y_i(2,:))))
RMSE = mean(sqrt(sum( (y_i-y_sim).*conj(y_i-y_sim) )))
%
% Error on testing trajectory
t_i = Y_test{ii,1};
y_i = iT(Y_test{ii,2});
t_i=t_i(4:end); y_i=y_i(:,4:end);
y_sim = iterate_map(N,length(t_i),y_i(:,1)); 
figure(52); clf;
subplot(211); hold on; grid on; box on; 
title('Reduced dynamics (map): normal form - testing trajectory')
plot(t_i,real(y_i(1,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,real(y_sim(1,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\mathrm{Re}(z_1)$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(real(y_i(1,:)),imag(y_sim(1,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(real(y_sim(1,:)),imag(y_sim(1,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$\mathrm{Re}(z_1)$','Interpreter','latex')
ylabel('$\mathrm{Im}(z_1)$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
RMSE = mean(sqrt(sum( (y_i-y_sim).*conj(y_i-y_sim) )))
xlim([-1 1]*max(abs(y_i(1,:))))
ylim([-1 1]*max(abs(y_i(2,:))))
%% Overall ROM performance

coordplot = 1; i_test = 1;
% Error on testing trajectory
t_i = X_test{i_test,1};
x_i = X_test{i_test,2};
y_i = iT(transpose(V_ort_data)*x_i);
y_ROM = iterate_map(N,length(t_i),y_i(:,1));
x_ROM = SSM_func(T(y_ROM));
figure(42); clf;
subplot(211); hold on; grid on; box on; 
title('Overall ROM performance (map)')
plot(t_i,x_i(coordplot,:),'k','Linewidth',2,'DisplayName','Testing Trajectory')
plot(t_i,x_ROM(coordplot,:),'r:','Linewidth',2,'DisplayName','ROM Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$q_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(real(y_i(1,:)),imag(y_i(2,:)),'k','Linewidth',2,'DisplayName','Testing Trajectory')
plot(real(y_ROM(1,:)),imag(y_ROM(2,:)),'r:','Linewidth',2,'DisplayName','ROM Trajectory')
xlabel('$\mathrm{Re}(z_1)$','Interpreter','latex')
ylabel('$\mathrm{Im}(z_1)$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([-1 1]*max(abs(y_i(1,:))))
ylim([-1 1]*max(abs(y_i(2,:))))
RRMSE = mean(sqrt(sum( (x_i-x_ROM).^2 )))/max(sqrt(sum(x_i.^2)))*100 % Percentage error based on the max ||x||


%% Nromal form and backbone curves
N_info = Maps_info.N;
[damp,freq] = nonres_normalform(N_info.coeff,N_info.exponents,t_i(2)-t_i(1));
figure(100); clf;
backbonecurves(damp,freq,SSM_func,T,coordplot,abs(y_i(1,1)),'norm');

%% Vector Field Dynamics
% Dynamics identification

PolyOrd = 10; Nfolds = 5; Nregvals = 20;
[R,iT,N,T,Maps_info] = IMdynamics_flow(Y_train,'R_PolyOrd',PolyOrd','style','normalform');

%%
% Error of the dynamics

% Error on training trajectory
t_i = Y_train{ii,1}; 
y_i = Y_train{ii,2}; 
t_i=t_i(4:end); y_i=y_i(:,4:end);
[t_sim,y_sim] = ode45(@(t,x) R(x),t_i,y_i(:,1)); y_sim = transpose(y_sim);
figure(61); clf;
subplot(211); hold on; grid on; box on;
title('Reduced dynamics (flow): default - training trajectory')
plot(t_i,y_i(1,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,y_sim(1,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(y_i(1,:),y_i(2,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(y_sim(1,:),y_sim(2,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$\eta_1$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([-1 1]*max(abs(y_i(1,:))))
ylim([-1 1]*max(abs(y_i(2,:))))
RMSE = mean(sqrt(sum( (y_i-y_sim).^2 )))

% Error on testing trajectory
t_i = Y_test{ii,1};
y_i = Y_test{ii,2};
t_i=t_i(4:end); y_i=y_i(:,4:end);
[~,y_sim] = ode45(@(t,x) R(x),t_i,y_i(:,1)); y_sim = transpose(y_sim);
figure(72); clf;
subplot(211); hold on; grid on; box on;
title('Reduced dynamics (flow): default - testing trajectory')
plot(t_i,y_i(1,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,y_sim(1,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(y_i(1,:),y_i(2,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(y_sim(1,:),y_sim(2,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$\eta_1$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([-1 1]*max(abs(y_i(1,:))))
ylim([-1 1]*max(abs(y_i(2,:))))
RMSE = mean(sqrt(sum( (y_i-y_sim).^2 )))

%% Modal
% Error on training trajectory
t_i = Y_train{ii,1}; 
y_i = iT(Y_train{ii,2}); 
[~,y_sim] = ode45(@(t,x) N(x),t_i,y_i(:,1)); y_sim = transpose(y_sim);
figure(81); clf;
subplot(211); hold on; grid on; box on; 
title('Reduced dynamics (flow): normal form - training trajectory')
plot(t_i,real(y_i(1,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,real(y_sim(1,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\mathrm{Re}(z_1)$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(real(y_i(1,:)),imag(y_sim(1,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(real(y_sim(1,:)),imag(y_sim(1,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$\mathrm{Re}(z_1)$','Interpreter','latex')
ylabel('$\mathrm{Im}(z_1)$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
RMSE = mean(sqrt(sum( (y_i-y_sim).*conj(y_i-y_sim) )))
xlim([-1 1]*max(abs(y_i(1,:))))
ylim([-1 1]*max(abs(y_i(2,:))))

% Error on testing trajectory
t_i = Y_test{ii,1};
y_i = iT(Y_test{ii,2});
[t_sim,y_sim] = ode45(@(t,x) N(x),t_i,y_i(:,1)); y_sim = transpose(y_sim);
figure(82); clf;
subplot(211); hold on; grid on; box on; 
title('Reduced dynamics (flow): normal form - testing trajectory')
plot(t_i,real(y_i(1,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,real(y_sim(1,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\mathrm{Re}(z_1)$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(real(y_i(1,:)),imag(y_sim(1,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(real(y_sim(1,:)),imag(y_sim(1,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$\mathrm{Re}(z_1)$','Interpreter','latex')
ylabel('$\mathrm{Im}(z_1)$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
RMSE = mean(sqrt(sum( (y_i-y_sim).*conj(y_i-y_sim) )))
xlim([-1 1]*max(abs(y_i(1,:))))
ylim([-1 1]*max(abs(y_i(2,:))))

%% Overall ROM performance
coordplot = 1; i_test = 1;

% Error on testing trajectory
t_i = X_test{i_test,1};
x_i = X_test{i_test,2};
y_i = iT(transpose(V_ort_data)*x_i);
[~,y_ROM] = ode45(@(t,x) N(x),t_i,y_i(:,1)); y_ROM = transpose(y_ROM);
x_ROM = SSM_func(T(y_ROM));
figure(92); clf;
subplot(211); hold on; grid on; box on; title('Overall ROM performance - flow')
plot(t_i,x_i(coordplot,:),'k','Linewidth',2,'DisplayName','Testing Trajectory')
plot(t_i,x_ROM(coordplot,:),'r:','Linewidth',2,'DisplayName','ROM Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$q_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(real(y_i(1,:)),imag(y_i(2,:)),'k','Linewidth',2,'DisplayName','Testing Trajectory')
plot(real(y_ROM(1,:)),imag(y_ROM(2,:)),'r:','Linewidth',2,'DisplayName','ROM Trajectory')
xlabel('$\mathrm{Re}(z_1)$','Interpreter','latex')
ylabel('$\mathrm{Im}(z_1)$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([-1 1]*max(abs(real(y_i(1,:)))))
ylim([-1 1]*max(abs(imag(y_i(1,:)))))
RRMSE = mean(sqrt(sum( (x_i-x_ROM).^2 )))/max(sqrt(sum(x_i.^2)))*100 % Percentage error based on the max ||x||

%% Normal form and backbone curves
N_info = Maps_info.N;
[damp,freq] = nonres_normalform(N_info.coeff,N_info.exponents);
figure(100);
backbonecurves(damp,freq,SSM_func,T,coordplot,abs(y_sim(1,1)),'norm');
legend('Map RD','Flow RD','Location','SE')

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
load Slow4DSSM_InitialConditions.mat

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

% Takens
X_traj = X_sim; opt_tak = 0; obs_coord = ndof+[1 2];%:2:N
if opt_tak == 1
for ii = 1:size(X_sim,1)
    X_i = X_sim{ii,2}; X_traj{ii,2} = X_i(obs_coord,:);
end
    over_embd = 500;
else
    over_embd = 0;
end
% Data embedding
SSM_dim = 4;
[X_traj,opts_embd] = coordinates_embedding(X_traj,SSM_dim,'OverEmbedding',over_embd);
opts_embd

% Split in training and testing dataset
ind_test = [1 5];
X_train = X_traj;
X_test = X_traj(ind_test,:);
X_train(ind_test,:)=[];
%%
% Polynomial degree of the parametrization
ParamDeg = 4; 
% Eigenspace
V = V_real(:,1:SSM_dim);
% Optimization
[V_ort_data,SSM_func,IM_para_info] = IMparametrization(X_traj,SSM_dim,ParamDeg,V);
%%
%ProjC = sum((V_ort'*V_ort_data).^2,2)
%sum(ProjC(SSM_dim+1:end))/length(ProjC(SSM_dim+1:end))

% Plot SSM Data

figure(2); clf; hold on; grid on; box on; colororder(winter(2))
figure(3); clf; hold on; grid on; box on; colororder(winter(2))
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
Err_rec = []; rho_max = 0;
coordplot = [1 ndof+1 9];
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
    plot3(q1,q2,q3,'Linewidth',.1)
    figure(3)
    q1 = y(1,:); q2 = y(2,:); q3 = y(3,:);
    plot3(q1,q2,q3,'Linewidth',2)
end
figure(2)
xlabel('$q_1$','Interpreter','latex')
ylabel('$\dot{q}_1$','Interpreter','latex')
zlabel('$q_5$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
figure(3)
xlabel('$\eta_1$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
zlabel('$\eta_3$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
ARSE = mean(Err_rec)
figure(7)

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
PolyOrd = 10; % Nfolds = 5; Nregvals = 20; % 'N_PolyOrd',7,'T_PolyOrd',9,)%,'n_folds',Nfolds,'l_vals',logspace(-6,0,Nregvals)
[R,iT,Nf,T,Maps_info] = IMdynamics_map(Y_train,'R_PolyOrd',PolyOrd,'style','normalform');
R_info = Maps_info.R;
N_info = Maps_info.N;
T_info = Maps_info.T;
iT_info = Maps_info.iT;
%%

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
ylim([-1 1]*max(abs(y_i(1,:))))
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(t_i,y_i(2,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,y_sim(2,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
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
ylim([-1 1]*max(abs(y_i(1,:))))
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
subplot(212); hold on; grid on; box on;
plot(t_i,y_i(2,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,y_sim(2,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
RMSE = mean(sqrt(sum( (y_i-y_sim).^2 )))

%% Modal
% Error on training trajectory
t_i = Y_train{ii,1}; 
y_i = iT(Y_train{ii,2}); 
y_sim = iterate_map(Nf,length(t_i),y_i(:,1)); 
figure(51); clf;
subplot(211); hold on; grid on; box on;
title('Reduced dynamics (map): default - training trajectory')
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
%
% Error on testing trajectory
t_i = Y_test{ii,1};
y_i = iT(Y_test{ii,2});
y_sim = iterate_map(Nf,length(t_i),y_i(:,1)); 
figure(52); clf;
subplot(211); hold on; grid on; box on;
title('Reduced dynamics (map): default - testing trajectory')
plot(t_i,real(y_i(1,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,real(y_sim(1,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\mathrm{Re}(z_1)$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
legend
xlim([t_i(1) t_i(end)])
subplot(212); hold on; grid on; box on;
% plot3(real(y_i(1,:)),imag(y_sim(1,:)),real(y_i(2,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
% plot3(real(y_sim(1,:)),imag(y_sim(1,:)),real(y_sim(2,:)),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
plot(t_i,real(y_i(2,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,real(y_sim(2,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\mathrm{Re}(z_2)$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
RMSE = mean(sqrt(sum( (y_i-y_sim).*conj(y_i-y_sim) )))
%view(3)

%% Overall ROM performance

coordplot = 1; i_test = 1;
% Error on testing trajectory
t_i = X_test{i_test,1};
x_i = X_test{i_test,2};
y_i = iT(transpose(V_ort_data)*x_i);
y_ROM = iterate_map(Nf,length(t_i),y_i(:,1));
x_ROM = SSM_func(T(y_ROM));
figure(42); clf;
subplot(211); hold on; grid on; box on;
title('Reduced dynamics (map): normal form - training trajectory')
plot(t_i,x_i(coordplot,:),'k','Linewidth',2,'DisplayName','Testing Trajectory')
plot(t_i,x_ROM(coordplot,:),'r:','Linewidth',2,'DisplayName','ROM Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$q_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(t_i,x_i(2,:),'k','Linewidth',2,'DisplayName','Testing Trajectory')
plot(t_i,x_ROM(2,:),'r:','Linewidth',2,'DisplayName','ROM Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])

RRMSE = mean(sqrt(sum( (x_i-x_ROM).^2 )))/max(sqrt(sum(x_i.^2)))*100 % Percentage error based on the max ||x||

%% Surface backbone plot
[damp,freq] = nonres_normalform(N_info.coeff,N_info.exponents,t_i(2)-t_i(1));
backbonesurfaces(damp,freq,abs(y_i([1 2],1)),'norm')

%% Reduced Dynamics

% Dynamics identification
PolyOrd = 10; 
[R,iT,Nf,T,Maps_info] = IMdynamics_flow(Y_train,'R_PolyOrd',PolyOrd,'style','normalform');
N_info = Maps_info.N;
%%

% Error on training trajectory
t_i = Y_train{ii,1};
y_i = Y_train{ii,2};
[t_sim,y_sim] = ode45(@(t,x) R(x),t_i,y_i(:,1)); y_sim = transpose(y_sim);
figure(71); clf;
subplot(211); hold on; grid on; box on;
title('Reduced dynamics (flow): default - training trajectory')
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
plot(t_i,y_i(2,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,y_sim(2,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
RMSE = mean(sqrt(sum( (y_i-y_sim).^2 )))

% Error on testing trajectory
t_i = Y_test{ii,1};
y_i = Y_test{ii,2};
[t_sim,y_sim] = ode45(@(t,x) R(x),t_i,y_i(:,1)); y_sim = transpose(y_sim);
figure(72); clf;
subplot(211); hold on; grid on; box on;
title('Reduced dynamics (flow): default - testing trajectory')
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
plot(t_i,y_i(2,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,y_sim(2,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
RMSE = mean(sqrt(sum( (y_i-y_sim).^2 )))

%% Normal form
% Error on training trajectory
t_i = Y_train{ii,1}; 
y_i = iT(Y_train{ii,2}); 
[t_sim,y_sim] = ode45(@(t,x) Nf(x),t_i,y_i(:,1)); y_sim = transpose(y_sim);
figure(81); clf;
subplot(211); hold on; grid on; box on;
title('Reduced dynamics (flow): normal form - training trajectory')
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
%
% Error on testing trajectory
t_i = Y_test{ii,1};
y_i = iT(Y_test{ii,2});
[t_sim,y_sim] = ode45(@(t,x) Nf(x),t_i,y_i(:,1)); y_sim = transpose(y_sim);
figure(82); clf;
subplot(211); hold on; grid on; box on;
title('Reduced dynamics (flow): normal form - testing trajectory')
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
RMSE = mean(sqrt(sum( (y_i-y_sim).*conj(y_i-y_sim) )))
%view(3)

%% Overall ROM performance

coordplot = 1; i_test = 1;
% Error on testing trajectory
t_i = X_test{i_test,1};
x_i = X_test{i_test,2};
y_i = iT(transpose(V_ort_data)*x_i);
[t_sim,y_ROM] = ode45(@(t,x) Nf(x),t_i,y_i(:,1)); y_ROM = transpose(y_ROM);
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
plot(t_i,x_i(2,:),'k','Linewidth',2,'DisplayName','Testing Trajectory')
plot(t_i,x_ROM(2,:),'r:','Linewidth',2,'DisplayName','ROM Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$q_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
RRMSE = mean(sqrt(sum( (x_i-x_ROM).^2 )))/max(sqrt(sum(x_i.^2)))*100 % Percentage error based on the max ||x||

%% Surface backbone plot
[damp,freq] = nonres_normalform(N_info.coeff,N_info.exponents);
backbonesurfaces(damp,freq,abs(y_i([1 2],1)),'norm')

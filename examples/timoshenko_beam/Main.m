% 2SSM from static deflection of a Timoschenko beam
% Vector field for data generattion taken from SSM Tool 1.0, see paper
% Ponsioen, Pedergnana, Haller - Automated Computation of Autonomous 
% Spectral Submanifolds for Nonlinear Modal Analysis

clear all
close all
format shortg
clc

% Input MCK matrices and dynamical model
ndof = 16;
MCKmatrices;
f_dyn = @(t,y) nonlinsys(y);

% Impose forcing on the tip (last dof)
f_vec = zeros(ndof,1); f_vec(end) = 1;
iMf_vec_ext = [zeros(ndof,1); M\f_vec];

% Linear Dynamics
A=[zeros(16) eye(16); -M\K -M\C];
[Vc,Dc]=eig(A); l=diag(Dc); [ictemp,indtemp]=sort(abs(imag(l)));
l=l(indtemp); %  Eigenvalues
mshape=Vc(:,indtemp); % Mode shapes

% Inetgration Setup
options = odeset('RelTol',1e-4,'AbsTol',1e-8);

%% Static simulations
% force amplitudes
f_static = [1e2:1e2:1e3 2e3:1e3:1e4 2e4:1e4:9e5];

U_static = zeros(ndof,length(f_static));
u0 = zeros(ndof,1);
fsolve_options = optimoptions('fsolve','Display','none');
figure(1); clf; hold on; grid on; box on;
ylabel('$f \, [$N$]$','Interpreter','latex')
xlabel('$u \, [$mm$]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
ylim([f_static(1) f_static(end)])
for ii = 1:length(f_static)
    f_stat_fun = @(x) M*[zeros(ndof) eye(ndof)]*f_dyn(0,[x; zeros(ndof,1)]) + f_vec*f_static(ii);
    U_static(:,ii) = fsolve(f_stat_fun,u0,fsolve_options);
    u0 = U_static(:,ii);
    U_static(:,ii) = fsolve(f_stat_fun,u0,fsolve_options);
    u0 = U_static(:,ii);
    plot(U_static(16,ii),f_static(ii),'.k','MarkerSize',12)
    drawnow;
end
plot(U_static(16,:),f_static,'k','Linewidth',1)

%% Decaying Trajectories from static deflections

% Initial conditions
X_static = [U_static; zeros(ndof,size(U_static,2))];
% Time settings: sampling and duration
t_scale =0.6; t_samp = t_scale/278; t_pers = 300; 
t = [0:t_samp:t_scale*t_pers];
[~,x1] = ode45(f_dyn,t,X_static(:,end),options);
[~,x2] = ode45(f_dyn,t,X_static(:,end-2),options);
figure(3); clf; 
subplot(211); hold on; grid on;
title('Tip Displacement')
plot(t,x1(:,ndof)); 
plot(t,x2(:,ndof)); 
xlabel('$t \, [$s$]$','Interpreter','latex')
ylabel('$u \, [$mm$]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
subplot(212); hold on; grid on;
title('Tip Velocity')
plot(t,x1(:,2*ndof)); 
plot(t,x2(:,2*ndof))
xlabel('$t \, [$s$]$','Interpreter','latex')
ylabel('$v \, [$mm/s$]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)

% Store trajectories in cells
observable_idx = ndof; % set on tip displacement, tip velocity has idx 2*ndof
X_traj_obs = cell(2,2);
X_traj_obs{1,1} = t; X_traj_obs{2,1} = t; 
X_traj_obs{1,2} = transpose(x1(:,observable_idx)); 
X_traj_obs{2,2} = transpose(x2(:,observable_idx)); 

%% SSMLearn

% clear previous variables for data generation
clearvars -except X_traj_obs
X_traj = X_traj_obs;
% Plot Data
figure(1); clf; hold on; grid on; box on;
plot(X_traj{1,1},X_traj{1,2},X_traj{2,1},X_traj{2,2})
xlabel('$t \, [$s$]$','Interpreter','latex')
ylabel('$u \, [$mm$]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)

%% Preprocessing

% Spectrogram
figure(10); clf; hold on; box on;
t = X_traj{1,1}; x = X_traj{1,2};
Nwin = round(length(t)/25); 
[Ss,Ws,Ts]=spectrogram(x,Nwin,round(Nwin*0.75),[],1./(t(2)-t(1)));
surf(Ts,Ws,log10(abs(Ss)))
ylim([0 40])
xlim([t(1) t(end)])
xlabel('$t \, [$s$]$','Interpreter','latex')
ylabel('$f \, [$Hz$]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
view(2)
shading interp

% DMD
X_traj_DMD = coordinates_embedding(X_traj,0,'OverEmbedding',29);
r = 2; Nwins = 100; ii = 1;
t_i = X_traj_DMD{ii,1};  t_i = t_i-t_i(1); X_traj_i = X_traj_DMD{ii,2}; 
[omega,lambda,t_eig] = sigprocess_DMD(X_traj_i(1:r*10,1:end-1),t_i,round(length(t_i)/Nwins),.99,r);
figure(11); clf; 
mmean_step = Nwins*2;
subplot(211); hold on; grid on; box on;
plot(t_eig,movmean(-real(omega(1:2:end,:))./abs(omega(1:2:end,:)),mmean_step,2)*100,'Linewidth',2)
xlabel('$t$','Interpreter','latex')
ylabel('$\xi \,[$\%$]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
subplot(212); hold on; grid on; box on;
plot(t_eig,movmean(imag(omega(1:2:end,:)),mmean_step,2),'Linewidth',2)
xlabel('$t$','Interpreter','latex')
ylabel('$f\,[$rad/s$]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])

%% Cut initial interval of the signals
t_cut = 4;
for ii = 1:2
    t = X_traj{ii,1}; x = X_traj{ii,2}; idx = sum(t<t_cut);
    X_traj{ii,1} = t(idx:end)-t(idx);  X_traj{ii,2} = x(:,idx:end); 
end
figure(2); hold on; grid on; box on;
plot(X_traj{1,1},X_traj{1,2},X_traj{2,1},X_traj{2,2})
xlabel('$t \, [$s$]$','Interpreter','latex')
ylabel('$u \, [$mm$]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)

%% Data Embedding (takens)
SSM_dim = 2; over_embd = 95;
[X_traj,opts_embd] = coordinates_embedding(X_traj,SSM_dim,'OverEmbedding',over_embd,'TimeStepping',1);
opts_embd

% Split in training and testing dataset
ind_test = [2];
X_train = X_traj;
X_test = X_traj(ind_test,:);
X_train(ind_test,:)=[];

%%
% Polynomial degree of the parametrization
ParamDeg = 1;
% Optimization
[V_ort_data,SSM_func] = IMparametrization(X_train,SSM_dim,ParamDeg);

% Plot and validation

coordplot = 1;

figure(2); clf; hold on; grid on; box on;
figure(3); clf; hold on; grid on; box on;
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
colors = winter(size(X_traj,1) );
Err_rec = []; Y = [];

for ii = 1:size(X_traj,1)
    % Simulate synthetic data, get observables and reconstruct test data
    t = X_traj{ii,1}; x = X_traj{ii,2};
    y = transpose(V_ort_data)*x; Y = [Y y]; x_rec = SSM_func(y);
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
    plot3(q1,q2,q3,'Linewidth',2,'Color',colors(ii,:))
    figure(3)
    q1 = y(1,:); q2 = y(2,:);
    plot(q1,q2,'Linewidth',2,'Color',colors(ii,:))
end
xlabel('$\eta_1$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
ARSE = mean(Err_rec)/max(sqrt(sum(x.^2)))*100

% Plot 3D SSM Surface
figure(2);
h = plot_2dSSM_surf(coordplot,Y,SSM_func,10,50,0);
xlabel('$\eta_1$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
ylabel('$u \, [$mm$]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)

%% Reduced Dynamics

% Arrange trajectories
Y_train = cell(size(X_traj,1)-length(ind_test),2);
for ii = 1:size(X_train,1)
    Y_train{ii,1} = X_train{ii,1};
    Y_train{ii,2} = transpose(V_ort_data)*X_train{ii,2};
end
Y_test  = cell(length(ind_test),2);
for ii = 1:size(X_test,1)
    Y_test{ii,1} = X_test{ii,1};
    Y_test{ii,2} = transpose(V_ort_data)*X_test{ii,2};
end

[R,iT,N,T,Maps_info] = IMdynamics_flow(Y_train,'R_PolyOrd',6,'style','normalform');

%%
% Error of the dynamics
ii = 1;
% Error on training trajectory
t_i = Y_train{ii,1};
y_i_nf = Y_train{ii,2};
[t_sim,y_sim] = ode45(@(t,x) R(x),t_i,y_i_nf(:,1)); y_sim = transpose(y_sim);
figure(61); clf;
subplot(211); hold on; grid on; box on;
title('Reduced dynamics (flow): default - training trajectory')
plot(t_i,y_i_nf(1,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_sim,y_sim(1,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(y_i_nf(1,:),y_i_nf(2,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(y_sim(1,:),y_sim(2,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$\eta_1$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([-1 1]*max(abs(y_i_nf(1,:))))
ylim([-1 1]*max(abs(y_i_nf(2,:))))
RMSE = mean(sqrt(sum( (y_i_nf-y_sim).^2 )))/max(sqrt(sum(y_i_nf.^2)))*100

% Error on testing trajectory
t_i = Y_test{ii,1};
y_i_nf = Y_test{ii,2};
[~,y_sim] = ode45(@(t,x) R(x),t_i,y_i_nf(:,1)); y_sim = transpose(y_sim);
figure(72); clf;
subplot(211); hold on; grid on; box on;
title('Reduced dynamics (flow): default - testing trajectory')
plot(t_i,y_i_nf(1,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_sim,y_sim(1,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\eta_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(y_i_nf(1,:),y_i_nf(2,:),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(y_sim(1,:),y_sim(2,:),'g:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$\eta_1$','Interpreter','latex')
ylabel('$\eta_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([-1 1]*max(abs(y_i_nf(1,:))))
ylim([-1 1]*max(abs(y_i_nf(2,:))))
RMSE = mean(sqrt(sum( (y_i_nf-y_sim).^2 )))/max(sqrt(sum(y_i_nf.^2)))*100

%% Normal Form
% Error on training trajectory
t_i = Y_train{ii,1};
y_i_nf = iT(Y_train{ii,2});
[~,y_sim] = ode45(@(t,x) N(x),t_i,y_i_nf(:,1)); y_sim = transpose(y_sim);
figure(81); clf;
subplot(211); hold on; grid on; box on;
title('Reduced dynamics (flow): normal form - training trajectory')
plot(t_i,real(y_i_nf(1,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,real(y_sim(1,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\rho$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(real(y_i_nf(1,:)),imag(y_sim(1,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(real(y_sim(1,:)),imag(y_sim(1,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$z_1$','Interpreter','latex')
ylabel('$z_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
RMSE = mean(sqrt(sum( (y_i_nf-y_sim).*conj(y_i_nf-y_sim) )))
xlim([-1 1]*max(abs(y_i_nf(1,:))))
ylim([-1 1]*max(abs(y_i_nf(2,:))))

% Error on testing trajectory
t_i = Y_test{ii,1};
y_i_nf = iT(Y_test{ii,2});
[t_sim,y_sim] = ode45(@(t,x) N(x),t_i,y_i_nf(:,1)); y_sim = transpose(y_sim);
figure(82); clf;
subplot(211); hold on; grid on; box on;
title('Reduced dynamics (flow): normal form - testing trajectory')
plot(t_i,real(y_i_nf(1,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(t_i,real(y_sim(1,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$\rho$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
subplot(212); hold on; grid on; box on;
plot(real(y_i_nf(1,:)),imag(y_sim(1,:)),'k','Linewidth',2,'DisplayName','Projected Trajectory')
plot(real(y_sim(1,:)),imag(y_sim(1,:)),'m:','Linewidth',2,'DisplayName','Simulated Trajectory')
xlabel('$z_1$','Interpreter','latex')
ylabel('$z_2$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
RMSE = mean(sqrt(sum( (y_i_nf-y_sim).*conj(y_i_nf-y_sim) )))
xlim([-1 1]*max(abs(y_i_nf(1,:))))
ylim([-1 1]*max(abs(y_i_nf(2,:))))

%% Overall ROM performance
coordplot = 1; i_test = 1;

% Error on testing trajectory
t_i = X_test{i_test,1};
x_i = X_test{i_test,2};
y_i_nf = iT(transpose(V_ort_data)*x_i);
[~,y_ROM_nf] = ode45(@(t,x) N(x),t_i,y_i_nf(:,1)); 
y_ROM_nf = transpose(y_ROM_nf);
x_ROM_nf = SSM_func(T(y_ROM_nf));

figure(92); clf; hold on; grid on; box on; 
title('Overall ROM performance - flow')
plot(t_i,x_i(coordplot,:),'k','Linewidth',2,'DisplayName','Testing Trajectory')
plot(t_i,x_ROM_nf(coordplot,:),'r:','Linewidth',2,'DisplayName','ROM Trajectory')
xlabel('$t$','Interpreter','latex')
ylabel('$u \, [$mm$]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_i(1) t_i(end)])
legend
RRMSE = mean(sqrt(sum( (x_i-x_ROM_nf).^2 )))/max(sqrt(sum(x_i.^2)))*100 % Percentage error based on the max ||x||

% Amplitude and phase errors
y_i_reg = transpose(V_ort_data)*x_i;
[~,y_ROM_reg] = ode45(@(t,x) R(x),t_i,y_i_reg(:,1)); 
y_ROM_reg = transpose(y_ROM_reg);
x_ROM_reg = SSM_func(y_ROM_reg);

t_ref = t_i;
x_ref = x_i(coordplot,:);
x_rec = [x_ROM_reg(coordplot,:); x_ROM_nf(coordplot,:)];
t_scale = 2*pi/11.02;
x_scale = max(abs(x_ref));

Err = ampph_error(x_rec,x_ref,t_ref);
figure(93); clf;
subplot(211); hold on; grid on; box on;
title('Relative Amplitude Error')
plot(Err.t,abs(Err.AmpError)/x_scale*100,'Linewidth',2)
xlabel('$t$','Interpreter','latex')
ylabel('$Err\,[$\%$]$','Interpreter','latex')
set(gca,'fontname','times')
%set(gca,'yscale','log')
set(gca,'fontsize',18)
xlim([t_ref(1) t_ref(end)])
ylim([1e-4 inf])

subplot(212); hold on; grid on; box on;
title('Relative Phase Error')
plot(Err.t,movmean(abs(Err.PhError),5,2)/t_scale*100,'Linewidth',2)
xlabel('$t$','Interpreter','latex')
ylabel('$Err\,[$\%$]$','Interpreter','latex')
set(gca,'fontname','times')
%set(gca,'yscale','log')
set(gca,'fontsize',18)
xlim([t_ref(1) t_ref(end)])
legend('Regression','Normal Form','Location','SE')

%% Normal form and backbone curves
N_info = Maps_info.N;
[damp,freq] = nonres_normalform(N_info.coeff,N_info.exponents);
figure(100); clf; hold on; grid on;
backbonecurves(damp,freq,SSM_func,T,coordplot,abs(y_i_nf(1,1)),'norm');
subplot(121); ylabel('$u \, [$mm$]$','Interpreter','latex')
subplot(122); ylabel('$u \, [$mm$]$','Interpreter','latex')
% Compare with DMD results
ii = 1;
t_i = X_traj{ii,1};
x_i = X_traj{ii,2};
y_i_nf = iT(transpose(V_ort_data)*x_i);
[~,y_ROM_nf] = ode45(@(t,x) N(x),t_i,y_i_nf(:,1)); 
rho_ROM = transpose(abs(y_ROM_nf(:,1)));
inst_freq = freq(rho_ROM); inst_damp = damp(rho_ROM);
inst_damp_rat = -inst_damp./inst_freq*100;
figure(11); 
subplot(211); plot(t_i+t_cut,inst_damp_rat,'Linewidth',2)
subplot(212); plot(t_i+t_cut,inst_freq,'Linewidth',2)
legend('DMD','Normal Form')

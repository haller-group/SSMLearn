function backbonesurfaces(a,w,max_rho,varargin)
% Plot instantaneous amplitudes and frequency/damping surfaces. The 
% amplitude are the parametrizing amplitude from 0 to max_rho. If an
% additional argument is given, the values of frequency and damping are
% nrmalized with respect to the linear limit

% Compute instantaneous damping and frequency
rho1_plot = linspace(0,(max_rho(1)),40); rho2_plot = linspace(0,(max_rho(2)),40);
[inst_amp1,inst_amp2] = meshgrid(rho1_plot,rho2_plot);
rr1 = transpose(inst_amp1(:)); rr2 = transpose(inst_amp2(:));
inst_damp_curve = a([rr1; rr2]); inst_freq_curve = w([rr1; rr2]);
if isempty(varargin)==0
inst_vars = [inst_freq_curve./inst_freq_curve(:,1); ...
    inst_damp_curve./inst_damp_curve(:,1)];
L_tick = min(inst_vars,[],2); U_tick = max(inst_vars,[],2); 
tick_v = 1+max(abs([L_tick U_tick]-1),[],2)*[-1 1];
color_map = pinkgreen(51);
figure; clf;
subplot(121); hold on; grid on; box on;
ZZ = reshape(transpose(inst_vars(1,:)),size(inst_amp1));
surf(inst_amp1,inst_amp2,ZZ)
shading interp
xlabel('$\rho_1$','Interpreter','latex')
ylabel('$\rho_2$','Interpreter','latex')
zlabel('$\omega_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
caxis(tick_v(1,:))
c = colorbar;
c.Location = 'northoutside';
c.Label.Interpreter = 'latex';
c.Label.String = '$\omega_1(\rho_1,\rho_2)/\omega_1(0,0)$';
colormap(color_map)
xlim([0 max_rho(1)])
ylim([0 max_rho(2)])
view(2)
subplot(122); hold on; grid on; box on;
ZZ = reshape(transpose(inst_vars(3,:)),size(inst_amp1));
surf(inst_amp1,inst_amp2,ZZ)
shading interp
xlabel('$\rho_1$','Interpreter','latex')
ylabel('$\rho_2$','Interpreter','latex')
zlabel('$\zeta_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
caxis(tick_v(3,:))
c = colorbar;
c.Location = 'northoutside';
c.Label.Interpreter = 'latex';
c.Label.String = '$\zeta_1(\rho_1,\rho_2)/\zeta_1(0,0)$';
colormap(color_map)
xlim([0 max_rho(1)])
ylim([0 max_rho(2)])
view(2)
figure; clf;
subplot(121); hold on; grid on; box on;
ZZ = reshape(transpose(inst_vars(2,:)),size(inst_amp1));
surf(inst_amp1,inst_amp2,ZZ)
shading interp
xlabel('$\rho_1$','Interpreter','latex')
ylabel('$\rho_2$','Interpreter','latex')
zlabel('$\omega_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
caxis(tick_v(2,:))
c = colorbar;
c.Location = 'northoutside';
c.Label.Interpreter = 'latex';
c.Label.String = '$\omega_2(\rho_1,\rho_2)/\omega_2(0,0)$';
colormap(color_map)
xlim([0 max_rho(1)])
ylim([0 max_rho(2)])
view(2)
subplot(122); hold on; grid on; box on;
ZZ = reshape(transpose(inst_vars(4,:)),size(inst_amp1));
surf(inst_amp1,inst_amp2,ZZ)
shading interp
xlabel('$\rho_1$','Interpreter','latex')
ylabel('$\rho_2$','Interpreter','latex')
zlabel('$\zeta_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
caxis(tick_v(4,:))
c = colorbar;
c.Location = 'northoutside';
c.Label.Interpreter = 'latex';
c.Label.String = '$\zeta_2(\rho_1,\rho_2)/\zeta_2(0,0)$';
colormap(color_map)
xlim([0 max_rho(1)])
ylim([0 max_rho(2)])
view(2)
else
inst_vars = [inst_freq_curve; inst_damp_curve];
figure; clf;
subplot(121); hold on; grid on; box on;
ZZ = reshape(transpose(inst_vars(1,:)),size(inst_amp1));
surf(inst_amp1,inst_amp2,ZZ)
shading interp
xlabel('$\rho_1$','Interpreter','latex')
ylabel('$\rho_2$','Interpreter','latex')
zlabel('$\omega_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
c = colorbar;
c.Location = 'northoutside';
c.Label.Interpreter = 'latex';
c.Label.String = '$\omega_1(\rho_1,\rho_2)$';
xlim([0 max_rho(1)])
ylim([0 max_rho(2)])
view(2)
subplot(122); hold on; grid on; box on;
ZZ = reshape(transpose(inst_vars(3,:)),size(inst_amp1));
surf(inst_amp1,inst_amp2,ZZ)
shading interp
xlabel('$\rho_1$','Interpreter','latex')
ylabel('$\rho_2$','Interpreter','latex')
zlabel('$\zeta_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
c = colorbar;
c.Location = 'northoutside';
c.Label.Interpreter = 'latex';
c.Label.String = '$\zeta_1(\rho_1,\rho_2)$';
xlim([0 max_rho(1)])
ylim([0 max_rho(2)])
view(2)
figure; clf;
subplot(121); hold on; grid on; box on;
ZZ = reshape(transpose(inst_vars(2,:)),size(inst_amp1));
surf(inst_amp1,inst_amp2,ZZ)
shading interp
xlabel('$\rho_1$','Interpreter','latex')
ylabel('$\rho_2$','Interpreter','latex')
zlabel('$\omega_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
c = colorbar;
c.Location = 'northoutside';
c.Label.Interpreter = 'latex';
c.Label.String = '$\omega_2(\rho_1,\rho_2)$';
xlim([0 max_rho(1)])
ylim([0 max_rho(2)])
view(2)
subplot(122); hold on; grid on; box on;
ZZ = reshape(transpose(inst_vars(4,:)),size(inst_amp1));
surf(inst_amp1,inst_amp2,ZZ)
shading interp
xlabel('$\rho_1$','Interpreter','latex')
ylabel('$\rho_2$','Interpreter','latex')
zlabel('$\zeta_1$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
c = colorbar;
c.Location = 'northoutside';
c.Label.Interpreter = 'latex';
c.Label.String = '$\zeta_2(\rho_1,\rho_2)$';
xlim([0 max_rho(1)])
ylim([0 max_rho(2)])
view(2) 
end
end
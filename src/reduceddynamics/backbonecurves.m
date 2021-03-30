function backbonecurves(a,w,SSM_func,T,coordplot,max_rho,varargin)
% Plot instantaneous amplitude and frequency/damping curves. The amplitude
% is defined as the maximum value reached by the physical coordinate in
% coordplot along a full rotation in the angle theta for each value of the
% parametrizing amplitude from 0 to max_rho. If an additional argument is
% given, the values of frequency and damping are normalized with respect to
% the linear limit

% Compute instantaneous damping and frequency
rho_plot = linspace(0,max_rho,100);
inst_damp_curve = a(rho_plot); inst_freq_curve = w(rho_plot);
% Compute instantaneous amplitude
theta_plot = linspace(0,2*pi,51); theta_plot = theta_plot(1:end-1);
[RR,TT] = meshgrid(rho_plot,theta_plot);
ZZ = RR.*exp(1i*TT);
z_eval = ZZ(:); z_eval = transpose(z_eval); z_eval =[z_eval; conj(z_eval)];
x_eval = SSM_func(T(z_eval)); amp_eval = transpose(x_eval(coordplot,:));
AA = reshape(amp_eval,size(RR));
% Amplitude as maximum value along theta
inst_amp = max(AA,[],1);
% Plot results
if isempty(varargin)==0
    if strcmp(varargin{:},'Hz')==1
        subplot(121); hold on; grid on; box on;
        plot(-inst_damp_curve./inst_freq_curve*100,inst_amp,'Linewidth',2)
        xlabel('$\xi \, [$\%$]$','Interpreter','latex')
        ylabel(['$x_{' num2str(coordplot) '}$'],'Interpreter','latex')
        set(gca,'fontname','times')
        set(gca,'fontsize',18)
        subplot(122); hold on; grid on; box on;
        plot(inst_freq_curve/2/pi,inst_amp,'Linewidth',2)
        xlabel('$\omega \, [$Hz$]$','Interpreter','latex')
        ylabel(['$x_{' num2str(coordplot) '}$'],'Interpreter','latex')
        set(gca,'fontname','times')
        set(gca,'fontsize',18)
    else
        subplot(121); hold on; grid on; box on;
        plot(inst_damp_curve./inst_damp_curve(1),inst_amp,'Linewidth',2)
        xlabel('$\zeta/\zeta(0)$','Interpreter','latex')
        ylabel(['$x_{' num2str(coordplot) '}$'],'Interpreter','latex')
        set(gca,'fontname','times')
        set(gca,'fontsize',18)
        subplot(122); hold on; grid on; box on;
        plot(inst_freq_curve/inst_freq_curve(1),inst_amp,'Linewidth',2)
        xlabel('$\omega/\omega(0)$','Interpreter','latex')
        ylabel(['$x_{' num2str(coordplot) '}$'],'Interpreter','latex')
        set(gca,'fontname','times')
        set(gca,'fontsize',18)
    end
else
    subplot(121); hold on; grid on; box on;
    plot(inst_damp_curve,inst_amp,'Linewidth',2)
    xlabel('$\zeta\, [$1/s$]$','Interpreter','latex')
    ylabel(['$x_{' num2str(coordplot) '}$'],'Interpreter','latex')
    set(gca,'fontname','times')
    set(gca,'fontsize',18)
    subplot(122); hold on; grid on; box on;
    plot(inst_freq_curve,inst_amp,'Linewidth',2)
    xlabel('$\omega\, [$rad/s$]$','Interpreter','latex')
    ylabel(['$x_{' num2str(coordplot) '}$'],'Interpreter','latex')
    set(gca,'fontname','times')
    set(gca,'fontsize',18)
end
end
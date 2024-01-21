function redSamp = ep_reduced_results(runid,sampStyle,ispolar,isomega,args1,args2,varargin)
% EP_REDUCED_RESULTS This function extract results of reduced dynamics at
% sampled forcing frequencies/amplitudes in continuation of equilibrium points
%
% REDSAMP =
% EP_REDUCED_RESULTS(RUNID,SAMPSTYLE,ISPOLAR,ISOMEGA,ARGS1,ARGS2,VARARGIN)
%
% runid:     id of continuation run of equilibrium points
% sampStyle: samping style of excitation frequencies
% ispolar:   coordinate representation of reduced dynamics
% isomega:   is the continuation parameter forcing frequency
% args1:     a cell arry of rho, namely {'rho1,'rho2',...,'rhom'} or a cell
%            array of real parts of z, {'Rez1','Rez2',...,'Rezm'}
% args2:     a cell arry of th, namely {'th1,'th2',...,'thm'} or a cell
%            array of imaginary parts of z, {'Imz1','Imz2',...,'Imzm'}
% varargin:  'plot-off', turn off plotting figures

%% extract results of reduced dynamics at sampled frequency/forcing
bd = coco_bd_read(runid);
m  = numel(args1);

if isempty(isomega)
    figure;
    coco_plot_bd(runid, 'om', 'eps');
    grid on; box on; 
    set(gca,'LineWidth',1.2);
    set(gca,'FontSize',14);
    xlabel('$\Omega$','interpreter','latex','FontSize',16);
    ylabel('$\epsilon$','interpreter','latex','FontSize',16);    
end

if strcmp(sampStyle,'cocoBD')
    om  = coco_bd_col(bd, 'om')';
    rho = zeros(numel(om),m);
    th  = zeros(numel(om),m);
    zRe = zeros(numel(om),m);
    zIm = zeros(numel(om),m);  
    if ispolar
        for k=1:m
            rho(:,k) = coco_bd_col(bd, args1{k})'; 
            th(:,k)  = coco_bd_col(bd, args2{k})';
        end
        state = rho.*exp(1j*th);
    else
        for k=1:m
            zRe(:,k) = coco_bd_col(bd, args1{k})'; 
            zIm(:,k) = coco_bd_col(bd, args2{k})';
        end
        rho = sqrt(zRe.^2+zIm.^2);
        th  = atan2(zIm,zRe);
        state = zRe+1j*zIm;
    end
    epsf = coco_bd_col(bd, 'eps')';
    if isempty(isomega)
        stab = nan;
    else
        stab = coco_bd_col(bd, 'eigs')'; 
        stab = real(stab);
        stab = all(stab<0,2);
    end
    SNidx = coco_bd_idxs(bd, 'SN');
    HBidx = coco_bd_idxs(bd, 'HB');
else
    if strcmp(sampStyle, 'uniform')
        sampLabs = coco_bd_labs(bd, 'UZ');     
    elseif strcmp(omegaSampStyle, 'cocoOut')
        sampLabs = coco_bd_labs(bd);
    end
    om  = coco_bd_vals(bd, sampLabs, 'om');
    rho = zeros(numel(om),m);
    th  = zeros(numel(om),m);
    zRe = zeros(numel(om),m);
    zIm = zeros(numel(om),m);
    if ispolar
        for k=1:m
            rho(:,k) = coco_bd_vals(bd, sampLabs, args1{k})'; 
            th(:,k)  = coco_bd_vals(bd, sampLabs, args2{k})';
        end
        state = rho.*exp(1j*th);
    else
        for k=1:m
            zRe(:,k) = coco_bd_vals(bd, sampLabs, args1{k})'; 
            zIm(:,k) = coco_bd_vals(bd, sampLabs, args2{k})';
        end
        rho = sqrt(zRe.^2+zIm.^2);
        th  = atan2(zIm,zRe);
        state = zRe+1j*zIm;
    end        
    epsf = coco_bd_vals(bd, sampLabs, 'eps');
    if isempty(isomega)
        stab = nan;
    else
        stab = coco_bd_vals(bd, sampLabs, 'eigs');
        stab = real(stab);
        stab = all(stab<0,2);
    end
    SNidx = [];
    HBidx = [];
end
% covert th such that it is in [0,2pi] and rho such it is postive
th             = mod(th,2*pi);
negRhoIdx      = find(rho<0);
rho(negRhoIdx) = -rho(negRhoIdx);
th(negRhoIdx)  = th(negRhoIdx)-pi;

redSamp = struct();
redSamp.om  = om;
redSamp.z   = state;
redSamp.rho = rho;
redSamp.th  = th;
redSamp.ep  = epsf;
redSamp.st  = stab;
redSamp.SNidx = SNidx;
redSamp.HBidx = HBidx;

if numel(varargin)==0
    
%% plot continuation results in normal coordinates
thm = struct( 'special', {{'SN', 'HB'}});
thm.SN = {'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'cyan', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'cyan', 'MarkerFaceColor', 'white'};
thm.HB = {'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 's', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white'};
figure;
for k=1:m
    subplot(m,1,k)
    if ispolar
        rhok = strcat('rho',num2str(k));
        if isempty(isomega)
            coco_plot_bd(thm, runid, 'om', 'eps', rhok, @(x) abs(x)); hold on
        else
            if isomega
                coco_plot_bd(thm, runid, 'om', rhok, @(x) abs(x)); hold on
            else
                coco_plot_bd(thm, runid, 'eps', rhok, @(x) abs(x)); hold on
            end
        end
    else
        Rezk = strcat('Rez',num2str(k));
        Imzk = strcat('Imz',num2str(k));
        if isempty(isomega)
            coco_plot_bd(thm, runid, 'om', 'eps', {Rezk,Imzk}, @(x,y) sqrt(x.^2+y.^2)); hold on
        else
            if isomega
                coco_plot_bd(thm, runid, 'om', {Rezk,Imzk}, @(x,y) sqrt(x.^2+y.^2)); hold on
            else
                coco_plot_bd(thm, runid, 'eps', {Rezk,Imzk}, @(x,y) sqrt(x.^2+y.^2)); hold on
            end
        end
    end
    grid on; box on; 
    set(gca,'LineWidth',1.2);
    set(gca,'FontSize',14);
    if isempty(isomega)
        xlabel('$\Omega$','interpreter','latex','FontSize',16);
        ylabel('$\epsilon$','interpreter','latex','FontSize',16);
        zlabel(strcat('$\rho_',num2str(k),'$'),'interpreter','latex','FontSize',16);
    else
        if isomega
            xlabel('$\Omega$','interpreter','latex','FontSize',16);
        else
            xlabel('$\epsilon$','interpreter','latex','FontSize',16);
        end
        ylabel(strcat('$\rho_',num2str(k),'$'),'interpreter','latex','FontSize',16);
    end
end


figure;
for k=1:m
    subplot(m,1,k)
    if ispolar
        rhok = strcat('rho',num2str(k));
        thk  = strcat('th',num2str(k));
        if isempty(isomega)
            coco_plot_bd(thm, runid, 'om', 'eps', {thk, rhok}, @(x,y) x+0.5*(1-sign(y))*pi); hold on
        else
            if isomega
                coco_plot_bd(thm, runid, 'om', {thk, rhok}, @(x,y) x+0.5*(1-sign(y))*pi); hold on
            else
                coco_plot_bd(thm, runid, 'eps', {thk, rhok}, @(x,y) x+0.5*(1-sign(y))*pi); hold on
            end
        end
    else
        Rezk = strcat('Rez',num2str(k));
        Imzk = strcat('Imz',num2str(k));
        if isempty(isomega)
            coco_plot_bd(thm, runid, 'om', 'eps', {Rezk,Imzk}, @(x,y) atan2(y,x)); hold on
        else
            if isomega
                coco_plot_bd(thm, runid, 'om', {Rezk,Imzk}, @(x,y) atan2(y,x)); hold on
            else
                coco_plot_bd(thm, runid, 'eps', {Rezk,Imzk}, @(x,y) atan2(y,x)); hold on
            end
        end
    end
    grid on; box on; 
    set(gca,'LineWidth',1.2);
    set(gca,'FontSize',14);
    if isempty(isomega)
        xlabel('$\Omega$','interpreter','latex','FontSize',16);
        ylabel('$\epsilon$','interpreter','latex','FontSize',16);
        zlabel(strcat('$\theta_',num2str(k),'$'),'interpreter','latex','FontSize',16);
    else
        if isomega
            xlabel('$\Omega$','interpreter','latex','FontSize',16);
        else
            xlabel('$\epsilon$','interpreter','latex','FontSize',16);
        end
        ylabel(strcat('$\theta_',num2str(k),'$'),'interpreter','latex','FontSize',16);
    end
end

end

end

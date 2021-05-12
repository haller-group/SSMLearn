%% Geometrically nonlinear von Karman beam
% Finite element model from the following reference:
% 
% Jain, S., Tiso, P., & Haller, G. (2018). Exact nonlinear model reduction for 
% a von Kármán beam: slow-fast decomposition and spectral submanifolds. _Journal 
% of Sound and Vibration_, _423_, 195–211. <https://doi.org/10.1016/J.JSV.2018.01.049 
% https://doi.org/10.1016/J.JSV.2018.01.049>
% 
% Finite element code taken from the following package:
% 
% Jain, S., Marconi, J., Tiso P. (2020). YetAnotherFEcode (Version v1.1). Zenodo. 
% <http://doi.org/10.5281/zenodo.4011282 http://doi.org/10.5281/zenodo.4011282>
% 
clear all
close all
clc

nElementsv = [4 8 16 32 64 128 256];
forcingv = 5e-4*[4 8 12]; 
figure(100); hold on; grid on; 
set(gca,'fontname','helvetica')
set(gca,'fontsize',18)
xlabel('$\Omega \, [$rad/s$]$','Interpreter','latex')
ylabel('$u \, [$m$]$','Interpreter','latex')
colors = colororder; colororder(colors(1:length(nElementsv),:))
% *system parameters*
for f_coeff = forcingv
for nElements = nElementsv
    
%figures to keep
figs2keep = [100]; all_figs = findobj(0, 'type', 'figure');
delete(setdiff(all_figs, figs2keep));
clearvars -except nElementsv forcingv f_coeff nElements

%% 

epsilon = 1;
%% generate model

[M,C,K,fnl,f_0,outdof] = build_model(nElements);
n = length(M);
%% Dynamical system setup 
% We consider the forced system
% 
% $$\mathbf{M}\ddot{\mathbf{x}}+\mathbf{C}\dot{\mathbf{x}}+\mathbf{K}\mathbf{x}+\mathbf{f}(\mathbf{x},\dot{\mathbf{x}})=\epsilon\mathbf{f}^{ext}(\mathbf{\Omega}t),$$
% 
% which can be written in the first-order form as 
% 
% $$\mathbf{B}\dot{\mathbf{z}}	=\mathbf{Az}+\mathbf{F}(\mathbf{z})+\epsilon\mathbf{F}^{ext}(\mathbf{\phi}),\\\dot{\mathbf{\phi}}	
% =\mathbf{\Omega}$$
% 
% where
% 
% $\mathbf{z}=\left[\begin{array}{c}\mathbf{x}\\\dot{\mathbf{x}}\end{array}\right],\quad\mathbf{A}=\left[\begin{array}{cc}-\mathbf{K} 
% & \mathbf{0}\\\mathbf{0} & \mathbf{M}\end{array}\right],\mathbf{B}=\left[\begin{array}{cc}\mathbf{C} 
% & \mathbf{M}\\\mathbf{M} & \mathbf{0}\end{array}\right],\quad\quad\mathbf{F}(\mathbf{z})=\left[\begin{array}{c}\mathbf{-\mathbf{f}(\mathbf{x},\dot{\mathbf{x}})}\\\mathbf{0}\end{array}\right],\quad\mathbf{F}^{ext}(\mathbf{z},\mathbf{\phi})=\left[\begin{array}{c}\mathbf{f}^{ext}(\mathbf{\phi})\\\mathbf{0}\end{array}\right]$.

DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
% set(DS.Options,'Emax',5,'Nmax',10,'notation','tensor')
%% 
% We assume periodic forcing of the form
% 
% $$\mathbf{f}^{ext}(\phi) = \mathbf{f}_0\cos(\phi)=\frac{\mathbf{f}_0}{2}e^{i\phi} 
% + \frac{\mathbf{f}_0}{2}e^{-i\phi}  $$
% 
% Fourier coefficients of Forcing

kappas = [-1; 1];
coeffs = f_coeff*[f_0 f_0]/2; % 7.5*
DS.add_forcing(coeffs, kappas,epsilon);
%% Linear Modal analysis and SSM setup

[V,D,W] = DS.linear_spectral_analysis();
%% 
% *Choose Master subspace (perform resonance analysis)*

S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
% set(S.Options, 'reltol', 0.1,'notation','tensor')
masterModes = [1,2]; 
S.choose_E(masterModes);
%% Forced response curves using SSMs
% Obtaining *forced response curve* in reduced-polar coordinate

order = 9; % Approximation order
%% 
% setup options

set(S.Options, 'reltol', 1,'IRtol',0.02,'notation', 'multiindex','contribNonAuto',true)
set(S.FRCOptions, 'nt', 2^7, 'nRho', 200, 'nPar', 200, 'nPsi', 100, 'rhoScale', 2 )
% set(S.FRCOptions, 'method','level set')
set(S.FRCOptions, 'method','continuation ep', 'z0', 1e-4*[1; 1]) % 'level set' 
set(S.FRCOptions, 'outdof',outdof)
%% 
% choose frequency range around the first natural frequency

omega0 = imag(S.E.spectrum(1));
omegaRange = omega0*[0.9 1.1];%[30 37];%
%% 
% extract forced response curve

FRC = S.extract_FRC('freq',omegaRange,order);
figFRC = gcf;
%%
figure(100);
%open('FRC_el.fig'); 
h = plot([FRC.Omega],[FRC.Aout],'Linewidth',1,'DisplayName',['Nr. Elements ' num2str(nElements)]);
if f_coeff~= forcingv(1); h.Annotation.LegendInformation.IconDisplayStyle = 'off';  end
legend('location','NW')
xlim(omegaRange)
%drawnow;
end
end

saveas(gcf, ['FRC_ElementsAnalysis_O' num2str(order)], 'fig');

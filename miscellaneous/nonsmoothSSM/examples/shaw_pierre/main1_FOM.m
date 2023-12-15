% Shaw-Pierre example: data generation from the full-order model 

clearvars
close all
clc

%% setup

% Parameters
ndof = 2;     % number of masses
mass = 1;
stiffness = 1;
gamma = 0.5;    
damping = 0.3;
delta = 1e-3;      % friction 

% Mass, stiffness and damping matrices
M = mass*diag(ones(1,ndof));
K = stiffness*(2*eye(ndof) - diag(ones(1,ndof-1),-1) - ...
    diag(ones(1,ndof-1),1));
C = damping * (2*eye(ndof) - diag(ones(1,ndof-1),-1) - ...
    diag(ones(1,ndof-1),1));
C(1,1) = damping * 1;

% Conservative mode shapes and natural frequencies
[Vcons,Dcons] = eig(full(K),full(M));
dfull = diag(Dcons);
[~,pos] = sort(dfull); Vcons = Vcons(:,pos);
Vcons = Vcons*diag( 1./sqrt(diag( transpose(Vcons) * M * Vcons ) ) );

% Phase space modes: definition and sorting
V = [Vcons zeros(ndof); zeros(ndof) Vcons];
W = [transpose(Vcons)*full(M) zeros(ndof); ...
        zeros(ndof) transpose(Vcons)*full(M)];

% Damped eigenvalues
C_M = M\C; K_M = M\K; 
Amat = [zeros(ndof) eye(ndof); -K_M -C_M];
eigs_val = eig(Amat); 
[~,pos] = sort(abs(eigs_val));
damp0_vals = -real(eigs_val(pos(1:2:end)));
freq0_vals = abs(eigs_val(pos(1:2:end)));

% Setup nonsmooth system
% Switching conditions
sFunction = @(x) x(ndof+1); velVec = zeros(1,2*ndof); velVec(ndof+1) = 1;
DxsFunction = @(x) velVec;
% Vector field
f_smooth = @(t,x) [x(ndof+1:2*ndof,:); -K_M*x(1:ndof,:)-C_M*x(ndof+1:2*ndof,:)-[gamma*x(1,:).^3; zeros(ndof-1,size(x,2))]];
fPlus = @(t,x) f_smooth(t,x)-transpose(velVec)*delta;
fMinus = @(t,x) f_smooth(t,x)+transpose(velVec)*delta;
fFunction = @(t,x,l) 0.5*(1+l)*fPlus(t,x)+0.5*(1-l)*fMinus(t,x);

%% Fixed points due to friction (finite delta approach)

q1Equilibrium = roots([gamma 0 3/2 -delta]); %equilibrium for the negative case
q1Equilibrium(imag(q1Equilibrium)~=0) = [];
equilibria = [q1Equilibrium -q1Equilibrium]; 
equilibria = [equilibria; equilibria/2]; % first col for +delta (q1dot<0)
fixedpointMinus = [equilibria(:,1); 0; 0];
fixedpointPlus = [equilibria(:,2); 0; 0];

% Additional stiffness terms when linearizing to new equilibria
KIII = [3*gamma*q1Equilibrium^2 0; 0 0];

%% Numerical integration of a nonsmooth trajectory

% Set initial condition and functions
x0 = [1 1 1 1].';
slowTime = 2*pi/freq0_vals(1);
nPers = ceil(freq0_vals(1)/damp0_vals(1));
pxPer = 100; 
tSpan = linspace(0,slowTime*nPers,1+pxPer*nPers);

% Integrate
[t,x,indSwitch] = odenonsmooth_1switchsurf(fFunction, sFunction, DxsFunction, [0 slowTime*nPers], x0, 'numDofs', ndof);
x = transpose(x);
% Plot in phase space
figure
hold on 
cPlot = [1 2 3]; % coordinates to plot
classicColors = colororder;
plot3(x(cPlot(1),:),x(cPlot(2),:),x(cPlot(3),:),'c','Linewidth',2)
plot3(x(cPlot(1),indSwitch),x(cPlot(2),indSwitch),x(cPlot(3),indSwitch),'.','Color',classicColors(3,:),'MarkerSize',18)
plot3(x(cPlot(1),1),x(cPlot(2),1),x(cPlot(3),1),'.','Color',[0 0 0],'MarkerSize',18)

xlabel('$q_1$ [m]','interpreter','latex')
ylabel('$q_2$ [m]','interpreter','latex')
zlabel('$\dot{q}_1$ [m/s]','interpreter','latex');
view(-10,11)
set(gca,'FontSize',20)
grid on 
tLims = [0 ceil(t(indSwitch(end))/100)*100];

% Plots in time
customFigure('subPlot',[2 1]); 
subplot(211)
plot(t,x(1,:),'Linewidth',2)
plot(t(indSwitch),x(1,indSwitch),'.','Color',classicColors(3,:),'MarkerSize',18)
xlabel('$t$ [s]','interpreter','latex')
ylabel('$q_1$ [m]','interpreter','latex')
xlim(tLims)
subplot(212)
plot(t,x(3,:),'Linewidth',2)
plot(t(indSwitch),x(3,indSwitch),'.','Color',classicColors(3,:),'MarkerSize',18)
xlabel('$t$ [s]','interpreter','latex')
ylabel('$\dot{q}_1$ [m/s]','interpreter','latex');
xlim(tLims)

customFigure('subPlot',[2 1]); 
subplot(211)
plot(t,x(2,:),'Linewidth',2)
plot(t(indSwitch),x(2,indSwitch),'.','Color',classicColors(3,:),'MarkerSize',18)
xlabel('$t$ [s]','interpreter','latex')
ylabel('$q_2$ [m]','interpreter','latex')
xlim(tLims)
subplot(212)
plot(t,x(4,:),'Linewidth',2)
plot(t(indSwitch),x(4,indSwitch),'.','Color',classicColors(3,:),'MarkerSize',18)
xlabel('$t$ [s]','interpreter','latex')
ylabel('$\dot{q}_2$ [m/s]','interpreter','latex');
xlim(tLims)

%% Numerical integration smooth systems
[~,x_minus] = ode45(fMinus,tSpan, x0);
[~,x_plus] = ode45(fPlus,tSpan, x0);

% Save data for comparison with reduced order modeling

x_full_DataNonSmooth = cell(1,3);
x_full_DataNonSmooth{1,1} = transpose(t); 
x_full_DataNonSmooth{1,2} = x; 
x_full_DataNonSmooth{1,3} = indSwitch;

xData{1,1} = tSpan; xData{1,2} = transpose(x_minus); xData{1,3} = '-';
xData{2,1} = tSpan; xData{2,2} = transpose(x_plus); xData{2,3} = '+';

save('sim_data','x_full_DataNonSmooth','xData','fixedpointMinus',...
    'fixedpointPlus','M','C','K','delta','KIII','gamma','mass','stiffness','damping')

%% Non-autonomous system

cont_name = 'dynsystem_dry_friction_f';
loadVector = [1; 1]/sqrt(2);
forcingVectors = [zeros(ndof,1); M\loadVector];
forcingAmplitudes = [0.15 0.2 0.25]; % forcing amplitudes
% detunings = [.9 1.1]; % frequency span
detunings = [.5 1.5];
data_label = cell(1,length(forcingAmplitudes));

% Loop for each forcing amplitude
for iAmp = 1:length(forcingAmplitudes)

amp_cont_name = [cont_name num2str(iAmp)];
data_label{iAmp} = ['F = ' num2str(forcingAmplitudes(iAmp))];

% Initial guess for continuation
% Parameters
p0 = [.95*freq0_vals(1); forcingAmplitudes(iAmp)];
ampsLin = (K-eye(ndof)*p0(1)^2)\forcingVectors(ndof+1:end)*p0(2);

% Linear initial condition for numerical integration
ampLin = ampsLin(1);
if ampLin<delta+p0(2)
    disp('No forced solution expected...')
end
%
currPhase = -pi/2;
x0 = [ampsLin*cos(currPhase); -ampsLin*p0(1)*sin(currPhase)*1.1; 1; 0];

% Switching conditions
sFunction = @(x) x(ndof+1); velVec = zeros(1,2*ndof+2); velVec(ndof+1) = 1;
DxsFunction = @(x) velVec;

% Numerical integration of a decaying trajectory
forcedSystem = @(y, p, mode) dynsystem(y, p, mode,@(x) fMinus(0,x),@(x) fPlus(0,x),forcingVectors);  
fFunction = @(t,x,l) 0.5*(1+l)*forcedSystem(x,p0,'right')+...
                     0.5*(1-l)*forcedSystem(x,p0,'left');
[t,x,indSwitch] = odenonsmooth_1switchsurf( fFunction, sFunction, DxsFunction, [0 nPers*2*pi/p0(1)], x0, 'numDofs', ndof);
x = transpose(x);
% Plot
customFigure;
plot(t*p0(1)/(2*pi),x(ndof+1,:),'Linewidth',1)
ylabel('q̇_1 [m/s]');
xlabel('tΩ/(2π)');

%% Obtain orbit segment

% Left
x1 = x(:,indSwitch(end-3):indSwitch(end-2));
t1 = t(indSwitch(end-3):indSwitch(end-2));
% Right
x2 = x(:,indSwitch(end-2):indSwitch(end-1));
t2 = t(indSwitch(end-2):indSwitch(end-1));
customFigure('subPlot', [2 2]);
subplot(2,2,1)
plot(x1(1,:),x1(ndof+1,:),'Linewidth',2)
plot(x2(1,:),x2(ndof+1,:),'Linewidth',2)
xlabel('q_1 [m]');
ylabel('q̇_1 [m/s]');
subplot(2,2,2)
plot(t1-t1(1),x1(end-1,:),'Linewidth',2)
plot(t2-t1(1),x2(end-1,:),'Linewidth',2)
ylabel('cos(Ωt)');
xlabel('t [s]');
xlim([0 2*pi/p0(1)])

title(['\rm Ω = ' num2str(ceil(100*2*pi/t2(end))/100)])
subplot(2,2,[3,4]); hold on; grid on; box on;
plot(t1-t1(1),x1(ndof+1,:),'Linewidth',2)
plot(t2-t1(1),x2(ndof+1,:),'Linewidth',2)
xlabel('t [s]');
ylabel('q̇_1 [m/s]');
xlim([0 2*pi/p0(1)])
set(gca, 'fontname', 'helvetica')
set(gca, 'fontsize', 16)
legend('Continuation - Seg. 1','Continuation - Seg. 2','location','NW')

%% Numerical continuation of FRC

modes  = {'left' 'right'};
events = {'boundary' 'boundary'};
resets = {'switch' 'switch'};
t0 = {t1-t1(1) t2-t2(1)};
x0 = {transpose(x1) transpose(x2)};

% Event and reset functions

stop    = @(x,p,e) x(ndof+1);
jump    = @(x,p,r) x;

% Use segment-specific settings for the initial values of the
% discretization parameters NTST and NCOL of the COLL toolbox.

prob = coco_prob();
prob = coco_set(prob, 'hspo.orb.bvp.seg1.coll', 'NTST', 10, 'NCOL', 6);
prob = coco_set(prob, 'hspo.orb.bvp.seg2.coll', 'NTST', 20, 'NCOL', 4);
prob = ode_isol2hspo(prob, '', ...
    {forcedSystem, stop, jump}, ...
    modes, events, resets, t0, x0, {'ffre' 'famp'}, p0);
prob = coco_set(prob, 'cont','NAdapt', 5, 'PtMX', [200 300],'NPR', 5);
%     
%     
fprintf('\n Run=''%s'': Continue family of two-segment periodic orbits.\n', ...
    amp_cont_name);
% Add amplitudes and initial conditions
prob = coco_add_slot(prob, 'hspoFRC', @add_slot_hspo_FRC, [], 'bddat');
% Run continuation
bd = coco(prob, amp_cont_name, [], 1, 'ffre', detunings);

% Extract continuation results
perSols=coco_bd_col(bd, 'hspo.period');
plot_coord = 1;
amp_seg1 = coco_bd_col(bd, 'seg1.MAX|X|');
amp_seg2 = coco_bd_col(bd, 'seg2.MAX|X|');
ampSols = max([amp_seg1(plot_coord,:); amp_seg2(plot_coord,:) ]);
phaSols = coco_bd_col(bd, 'phase'); phaSols = phaSols(plot_coord,:);
staSols = sum(abs(coco_bd_col(bd, 'eigs'))>1)==0;
nSwitches = coco_bd_col(bd, 'switches');

FRC.(['F' num2str(iAmp)]) = struct('Freq',2*pi./perSols,'Amp',...
    ampSols,'Nf_Phs',phaSols,'Stab',staSols,'NSwitch',nSwitches);

%%
customFigure('figNumber',1000);
colors = colororder;
plotFRC(FRC,'color', [0 0 0]);
xlim(detunings)
ylim([0 max(max(ampSols))])
legend('location','NE')
xlabel('forcing frequency Ω [rad/s]');
ylabel('amplitude max|q_1| [m]');
freqs = [0:0.01:5]; figure(1000)
fre_unphys = [];
amp_unphys = [];
for iiAmp = 1:iAmp
    indPhys = (FRC.(['F' num2str(iiAmp)]).NSwitch~=2);
    fre_unphys = [fre_unphys FRC.(['F' num2str(iiAmp)]).Freq(indPhys)];
    amp_unphys = [amp_unphys FRC.(['F' num2str(iiAmp)]).Amp(indPhys)];
end
plot(fre_unphys,amp_unphys,'x','color', colors(7,:),'Linewidth',2,'DisplayName','Unphysical')
% Compute linear respose at which we have sticking
if plot_coord ==1
    area(freqs,delta./(freqs.^2),'FaceColor',[1 1 1]*0.4,'DisplayName','Sticking','FaceAlpha',.3)
end
ylim([0 0.5]);
drawnow; 
% saveas(gcf,[imageFolder num2str(-log10(delta)) 'FRCamp.pdf'])
% ylim([0 0.2]); saveas(gcf,[imageFolder num2str(-log10(delta)) 'FRCampUnp.pdf'])
% ylim([0 0.5]);
% ylim([0 1])
%%
customFigure('figNumber',1001);
colors = colororder;
plotFRC(FRC,'color', [0 0 0],'y', 'Phase');
xlim(detunings)
legend('location','NE')
xlabel('forcing frequency Ω [rad/s]');
ylabel('phase q_1');
drawnow;
if plot_coord == 1
    ylim([-pi 0])
    set(gca,'YTick',[-pi -3*pi/4 -pi/2 -pi/4 0])
    set(gca,'YTickLabel',{'-\pi', '-3\pi/4', '-\pi/2', '-\pi/4', '0'})
end
% saveas(gcf,[imageFolder num2str(-log10(delta)) 'FRCphase.pdf'])
%% Plot of solutions and compare with numerical integration
plotDof = 1;
labs_EP = coco_bd_labs(bd, 'EP');
labs_FP = coco_bd_labs(bd, 'FP');
for solNr = [min(labs_FP)-1 max(labs_EP)]
    % Extract continuation solution
    sol = hspo_read_solution([], amp_cont_name, solNr);
    x1 = transpose(sol.xbp{1}); x2 = transpose(sol.xbp{2});
    t1 = transpose(sol.tbp{1}); t2 = transpose(sol.tbp{2})+t1(end);
    
    % Numerical integration
    iSeg = 1; iPhase = 2;
    p0 = [2*pi/t2(end); forcingAmplitudes(iAmp)];
    fFunction = @(t,x,l) 0.5*(1+l)*forcedSystem(x,p0,'right')+...
                     0.5*(1-l)*forcedSystem(x,p0,'left');
    x = transpose(sol.xbp{iSeg}(iPhase,:));
    [t,x,indSwitch] = odenonsmooth_1switchsurf( fFunction, sFunction, DxsFunction, [0 100*2*pi/p0(1)], x, 'numDofs', ndof);

    x = transpose(x);
    switch length(indSwitch)
        case 1
            iniPoint = 1;
        case 2
            iniPoint = indSwitch(1);
        case 3
            if x(2,indSwitch(end-2)+1)<0
                iniPoint = indSwitch(end-2);
            else
                iniPoint = 1;
            end
        otherwise
            if x(2,indSwitch(end)-1)>0
                iniPoint = indSwitch(end-2);
            else
                iniPoint = indSwitch(end-3);
            end
    end
    
    % Plot
    customFigure('subPlot', [2 2]);
    subplot(2,2,1)
    plot(x1(plotDof,:),x1(ndof+plotDof,:),'Linewidth',2)
    plot(x2(plotDof,:),x2(ndof+plotDof,:),'Linewidth',2)
    xlabel(['q_' num2str(plotDof) ' [m]']);
    ylabel(['q̇_' num2str(plotDof) ' [m/s]']);
    subplot(2,2,2)
    plot(t1,x1(end-1,:),'Linewidth',2)
    plot(t2,x2(end-1,:),'Linewidth',2)
    ylabel('cos(Ωt)');
    xlabel('t [s]');
    xlim([0 2*pi/p0(1)])
    subplot(2,2,[3,4]); hold on; grid on; box on;
    
    title(['\rm Ω = ' num2str(ceil(100*2*pi/t2(end))/100)])
    plot(t-t(iniPoint),x(ndof+plotDof,:),'Color',colors(5,:),'Linewidth',3)
    plot(t1,x1(ndof+plotDof,:),'Color',colors(1,:),'Linewidth',2)
    plot(t2,x2(ndof+plotDof,:),'Color',colors(2,:),'Linewidth',2)
    xlabel('t [s]');
    ylabel(['q̇_' num2str(plotDof) ' [m/s]']);
    set(gca, 'fontname', 'helvetica')
    set(gca, 'fontsize', 16)
    xlim([0 2*pi/p0(1)])
    legend('Integration','Continuation - Seg. 1','Continuation - Seg. 2','location','NW')
%     if iAmp == length(forcingAmplitudes); saveas(gcf,[imageFolder num2str(-log10(delta)) 'Trajectory' num2str(solNr) '.pdf']); end
%     if iAmp == 1 && solNr == max(labs_EP); saveas(gcf,[imageFolder num2str(-log10(delta)) 'TrajectoryUnp' num2str(solNr) '.pdf']); end
end
end
FRC_FOM = FRC;


save('frc_data','FRC_FOM','forcingVectors','forcingAmplitudes','detunings','delta');

%% specific forced trajectory
cont_name = 'dynsystem_dry_friction_f';
loadVector = [1; 1]/sqrt(2);
forcingVectors = [zeros(ndof,1); M\loadVector];
% Parameters
p0 = [0.95; 0.2];

ampsLin = (K-eye(ndof)*p0(1)^2)\forcingVectors(ndof+1:end)*p0(2);
nPers = 10;

% Linear initial condition for numerical integration
ampLin = ampsLin(1);
if ampLin<delta+p0(2)
    disp('No forced solution expected...')
end
%
currPhase = -pi/2;
x0 = [ampsLin*cos(currPhase); -ampsLin*p0(1)*sin(currPhase)*1.1; 1; 0];

% Switching conditions
sFunction = @(x) x(ndof+1); velVec = zeros(1,2*ndof+2); velVec(ndof+1) = 1;
DxsFunction = @(x) velVec;

% Numerical integration of a decaying trajectory
forcedSystem = @(y, p, mode) dynsystem(y, p, mode,@(x) fMinus(0,x),@(x) fPlus(0,x),forcingVectors);  
fFunction = @(t,x,l) 0.5*(1+l)*forcedSystem(x,p0,'right')+...
                     0.5*(1-l)*forcedSystem(x,p0,'left');
[t,x,indSwitch] = odenonsmooth_1switchsurf( fFunction, sFunction, DxsFunction, [0 nPers*2*pi/p0(1)], x0, 'numDofs', ndof);
x = transpose(x);
% Plot
customFigure;
plot(t*p0(1)/(2*pi),x(1,:),'Linewidth',1)
ylabel('q_1');
xlabel('tΩ/(2π)');


x_full_forced = cell(1,4);
x_full_forced{1,1} = p0;
x_full_forced{1,2} = transpose(t); 
x_full_forced{1,3} = x; 
x_full_forced{1,4} = indSwitch;

save('full_forced_data','x_full_forced');

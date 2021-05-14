function FRC_data = getFRC_ddROM(ddROM,f_red,Wspan,coordplot)

% Definitions
N = ddROM.ReducedDynNormal;
SSM_func = ddROM.Param;
T = ddROM.CCfromNormal;
n = ddROM.Dim;
fr_coco = @(t,x,p) reim_transf(N(cc_transf(x(1,:)+1i*x(2,:))) - ...
          1i*p(2).*[exp(1i*2*pi*t./p(1)); exp(-1i*2*pi*t./p(1))]);
Tspan = sort(2*pi./Wspan);

FRC_data = struct();
for ii = 1:length(f_red)
p0 = [Tspan(end)*0.99; f_red(ii)];
bd_name = ['FRC_SSMLearn_F' num2str(ii) ];
% Initial guess
[~, x0] = ode45(@(t,x) fr_coco(t,x,p0), [0 350*p0(1)], zeros(n,1)); % Transients
[t0, x0] = ode45(@(t,x) fr_coco(t,x,p0), [0 p0(1)], x0(end,:)'); % Approximate periodic orbit

% Numerical continuation
prob = coco_prob();
prob = coco_set(prob, 'ode', 'autonomous', false);
coll_args = [{fr_coco}, {t0, x0, {'T' 'e'}, p0}];
prob = ode_isol2po(prob, '', coll_args{:});
prob = coco_set(prob, 'corr', 'ItMX', 80);
prob = coco_set(prob, 'coll', 'NTST', 120,'NCOL',4);
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'PtMX', [1 0]*40000,'NPR',10,'h_max',0.1);
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_glue(prob, 'glue', uidx(maps.T_idx), uidx(maps.p_idx(1)));
prob = po_mult_add(prob, 'po.orb');  % Store Floquet multipliers with bifurcation data
prob = coco_add_slot(prob, 'ode_cons_bddat', @add_slot_IP, [], 'bddat');
cont_args = {1, {'po.period' 'T' 'f_red'}, Tspan};

% Run coco
fprintf('\n Run=''%s'': Continue primary family of periodic orbits.\n', ...
  bd_name);
bd  = coco(prob,bd_name, [], cont_args{:});

% Extracts simulation results
T_vec=coco_bd_col(bd, 'po.period');
Z0=coco_bd_col(bd, 'X0');
ph_vec = unwrap(angle(Z0(1,:)+1i*Z0(2,:)));
rho_vec = abs(Z0(1,:)+1i*Z0(2,:));
stability = sum(abs(coco_bd_col(bd, 'eigs'))>1)==0;
    
% Compute instantaneous amplitude
theta_plot = linspace(0,2*pi,51); theta_plot = theta_plot(1:end-1);
[RR,TT] = meshgrid(rho_vec,theta_plot);
ZZ = RR.*exp(1i*TT);
z_eval = ZZ(:); z_eval = transpose(z_eval); z_eval =[z_eval; conj(z_eval)];
x_eval = SSM_func(T(z_eval)); amp_eval = transpose(x_eval(coordplot,:));
AA = reshape(amp_eval,size(RR)); amp_vec = max(abs(AA),[],1);

% Save and output
FRC_data.(['F' num2str(ii)]) = struct('Freq',2*pi./T_vec,'Amp',...
                amp_vec,'Nf_Amp',rho_vec,'Nf_Phs',ph_vec,'Stab',stability);
end

end

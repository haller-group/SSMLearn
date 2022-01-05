function FRC_data = computeFRC(IMInfo, RDInfo, reducedForcing, ...
    frequencySpan, amplitudeFunction)
% Computes forced frequency sweeps for every forcing parameter in
% reducedForcing within the frequency span in frequencySpan
% The non-autonomous dynamics has the form \dot{z} = N(z) + f(z,t)
% where f(z,t) = [ [f_1 f_2 ... f_{ndof}] * e^{\Omega t}, compl. conj.]
% where ndof = k/2 and k is dimension of reduced order model. The values
% f_1, f_2 can be complex, whose magnitude is the forcing amplitude on a
% specific mode. The matrix reducedForcing has dimension k/2 x nSweeps,
% whose rows store the forcing for each mode for each required sweep.

% Definitions and Inputs
N = RDInfo.conjugateDynamics.map;
T = RDInfo.transformation.map;
SSM_func = IMInfo.parametrization.map;
k = length(RDInfo.eigenvaluesLinPartFlow);
Tspan = sort(2*pi./frequencySpan);
if strcmp(RDInfo.dynamicsType,'flow') ~= 1
    error('Error. Forced responses are only available for flows.')
else
    ndof = k/2;
end
if size(reducedForcing,1) ~= ndof
    error(['Error. The number of forcing values for a single sweep is' ...
        'not equal to the number of modes.'])
end

% Run continuations
FRC_data = struct(); nSweeps = size(reducedForcing,2);
initialT = Tspan(end); Tspan(end) = Tspan(end)*1.05;
for idxSweep = 1:nSweeps
    forceVector = reducedForcing(:,idxSweep);
    fr_coco = @(t,x,p) transformationReIm( N( ...
        transformationComplexConj(x(1:ndof,:)+1i*x(ndof+1:end,:)) ) + ...
        transformationComplexConj(1i*forceVector*exp(1i*2*pi*t./p(1))) );
    
    p0 = initialT;
    bd_name = ['FRC_SSMLearn_F' num2str(idxSweep) ];
    
    % Initial guess
    % Transients
    [~, x0] = ode45(@(t,x) fr_coco(t,x,p0), [0 500*p0(1)], zeros(k,1)); 
    % Approximate periodic orbit
    [t0, x0] = ode45(@(t,x) fr_coco(t,x,p0), [0 p0(1)], x0(end,:)'); 
    
    % Numerical continuation
    prob = coco_prob();
    prob = coco_set(prob, 'ode', 'autonomous', false);
    coll_args = [{fr_coco}, {t0, x0, {'T'}, p0}];
    prob = ode_isol2po(prob, '', coll_args{:});
    prob = coco_set(prob, 'corr', 'ItMX', 80);
    prob = coco_set(prob, 'coll', 'NTST', 120,'NCOL',4);
    prob = coco_set(prob, 'cont', 'NAdapt', 1, 'PtMX', [1 0]*20000,'NPR',10);
    [data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
    maps = data.coll_seg.maps;
    prob = coco_add_glue(prob, 'glue', uidx(maps.T_idx), uidx(maps.p_idx(1)));
    prob = po_mult_add(prob, 'po.orb');  % Store Floquet multipliers
    prob = coco_add_slot(prob, 'ode_cons_bddat', @add_slot_IP, [], 'bddat');
    cont_args = {1, {'po.period' 'T'}, Tspan};
    
    % Run coco
    if ndof == 1
        disp(['Frequency sweep for the forcing value number ' ...
            num2str(idxSweep) ' ...'])
    else
        disp(['Frequency sweep for the set of forcing values number ' ...
            num2str(idxSweep) ' ...'])
    end
    bd  = coco(prob,bd_name, [], cont_args{:});
    
    % Extracts simulation results
    T_vec=coco_bd_col(bd, 'po.period');
    Z0=coco_bd_col(bd, 'X0');
    ph_vec = zeros(ndof,length(T_vec));
    rho_vec = zeros(ndof,length(T_vec));
    for idof = 1:ndof
        Z0_i = Z0(idof,:)+1i*Z0(idof+ndof,:);
        ph_vec(idof,:) = unwrap(angle(Z0_i)); rho_vec(idof,:) = abs(Z0_i);
    end
    stability = sum(abs(coco_bd_col(bd, 'eigs'))>1)==0;
    
    % Compute instantaneous amplitude
    if ndof == 1
        theta_plot = linspace(0,2*pi,51); theta_plot = theta_plot(1:end-1);
        [RR,TT] = meshgrid(rho_vec,theta_plot);
        ZZ = RR.*exp(1i*TT);
        zEval = ZZ(:); zEval = transformationComplexConj(transpose(zEval));
        ampEval = abs(amplitudeFunction( SSM_func( T(zEval) ) ));
        amp_vec = max(reshape(transpose(ampEval),size(RR)),[],1);
    else
        amp_vec = zeros(1,length(T_vec));
        for ii = 1:length(T_vec)
            [~, xyEval] = ode45(@(t,x) fr_coco(t,x,T_vec(ii)), ...
                [0 T_vec(ii)], Z0(:,ii));
            xyEval = transpose(xyEval); zEval = xyEval;
            for idof = 1:ndof
                zEval(idof,:) = xyEval(idof,:)+1i*xyEval(idof+ndof,:);
                zEval(idof+ndof,:) = conj(zEval(idof,:));
            end
            amp_vec(ii) = max(abs(amplitudeFunction( SSM_func( T(zEval) ...
                ))));
        end
    end
    
    % Save and output
    FRC_data.(['F' num2str(idxSweep)]) = struct('Freq',2*pi./T_vec,'Amp',...
        amp_vec,'Nf_Amp',rho_vec,'Nf_Phs',ph_vec,'Stab',stability);
end

end
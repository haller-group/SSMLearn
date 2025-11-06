function FRC_data = continuationFRCpo_mod(IMInfoF, RDInfoF, reducedForcing, ...
    frequencySpan, amplitudeFunction,outdofs)
%   FRC = continuationFRCpo(IMInfoF, RDInfoF, fRed, frequencySpan,amplitudeFunction)
%   Compute the forced response curves on a 2m-dim. SSM for each forcing
%   amplitude in fRed using numerical continuation of periodic orbits
%   within a given frequencySpan. Available only for flows.
%
%   INPUT
%   IMInfoF            struct            Time-periodic manifold.
%   RDInfoF            struct            Time-periodic reduced dynamics.
%   f_red              (1 x nAmp)        Normal form forcing constants.
%   amplitudeFunction  function handle   Map from y to (signed) scalar. Can
%                                        be used to predict forced response
%                                        for a given component or quantity.
%   OUTPUT
%   FRC           struct            Forced response curve with one field F1,
%                                   F2, ... for each forcing amplitude in
%                                   fRed

% Definitions and Inputs
if strcmp(RDInfoF.conjugacyStyle,'normalform') == 1
    N = RDInfoF.conjugateDynamics.map;
    T = RDInfoF.transformation.map;
else
    N = RDInfoF.reducedDynamics.map;
end
SSMFunction = IMInfoF.parametrization.map;
k = length(RDInfoF.eigenvaluesLinPartFlow);
Tspan = sort(2*pi./frequencySpan);
if strcmp(RDInfoF.dynamicsType,'flow') ~= 1
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
    fAmpls = abs(forceVector); fPhs = atan2(imag(forceVector),real(forceVector));
    if strcmp(RDInfoF.conjugacyStyle,'normalform') == 1
        fr_coco = @(t,x,p) transformationReIm(N(t,transformationComplexConj(x(1:ndof,:)+1i*x(ndof+1:end,:)),2*pi./p(1),fAmpls,fPhs));
    else
        fr_coco = @(t,x,p) N(t,x,2*pi./p(1),fAmpls,fPhs);
    end
    
    
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
    prob = coco_add_slot(prob, 'ode_cons_bddat1', @add_slot_Traj, [], 'bddat');
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
    Traj_vec=coco_bd_col(bd, 'Traj');
    if strcmp(RDInfoF.conjugacyStyle,'normalform') == 1
        Z0 = coco_bd_col(bd, 'X0'); 
        theta_vec = atan2(Z0(ndof+1:end,:),Z0(1:ndof,:));
        Z0 = transformationComplexConj( Z0(1:ndof,:)+1i*Z0(ndof+1:end,:) );
        rho_vec = abs(Z0(1:ndof,:)); 
    else
        Z0 = coco_bd_col(bd, 'X0'); rho_vec = 0*Z0; theta_vec = 0*Z0;
    end
    stability = sum(abs(coco_bd_col(bd, 'eigs'))>1)==0;
    
    % Compute instantaneous amplitude
    nSteps = length(T_vec); Omega = 2*pi./T_vec;
    u = zeros(length(outdofs),nSteps); uPhase = zeros(length(outdofs),nSteps);
    for ind = 1:nSteps
        iOmega = Omega(ind);
        normTimeEval = transpose(Traj_vec{ind}.NormalizedTime);
        timeEval = normTimeEval*T_vec(ind);
        cPhi = cos(normTimeEval*2*pi); sPhi = sin(normTimeEval*2*pi);
        if strcmp(RDInfoF.conjugacyStyle,'normalform') == 1
            zEvalTemp = reshape(Traj_vec{ind}.StateVector,k,length(timeEval));
            zEval = transformationComplexConj( zEvalTemp(1:ndof,:)+1i*zEvalTemp(ndof+1:end,:) );
            etaEval = T(timeEval,zEval,iOmega,fAmpls,fPhs);
        else
            etaEval = reshape(Traj_vec{ind}.StateVector,k,length(timeEval));
        end
        y = SSMFunction(timeEval,etaEval,iOmega,fAmpls,fPhs);
        yAmplitudeFunction = amplitudeFunction(y);
        u(:,ind) = max(abs(yAmplitudeFunction), [], 2);
        uPhase(:,ind) = phase( normTimeEval(2) * sum( yAmplitudeFunction.*(cPhi-1i*sPhi ) ,2));
    end
    
    % Save and output
    FRC_data.(['F' num2str(idxSweep)]) = struct('Freq',2*pi./T_vec,'Amp',...
        u,'Nf_Amp',rho_vec,'Nf_Phs',uPhase,'Nf_Phs2',theta_vec,'Stab',stability,'z_IC',Z0);
end

end
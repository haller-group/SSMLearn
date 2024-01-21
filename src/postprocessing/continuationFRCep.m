function varargout = continuationFRCep(IMInfoF, RDInfoF, epsilon, frequencySpan, amplitudeFunction, mFreqs, resonantModes, runid)
%   FRC = continuationFRCfp(IMInfoF, RDInfoF, fRed, frequencySpan,amplitudeFunction)
%   Compute the forced response curves on a 2m-dim. SSM for each forcing
%   amplitude in epsilon using numerical continuation of fixed points 
%   within a given omegaSpan. Available only for flows with polar form.
%
%   INPUT
%   IMInfoF            struct            Time-periodic manifold.
%   RDInfoF            struct            Time-periodic reduced dynamics.
%   epsilon            (1 x nAmp)        forcing amplitude of the forcing 
%                                        vector.
%   frequencySpan      (1 x 2)           frequency interval of interest
%   amplitudeFunction  function handle   Map from y to N signed amplitudes
%                                        of interest. Can be used to 
%                                        predict forced response for
%                                        given components or quantities.
%   mFreqs             (1 x nFreqs)      frequency ratio for internally
%                                        resonant ROMs
%   resonantModes      (1 x 2nFreqs)     resonant reduced coordinates 
%   runid              string            name for coco runs
%   OUTPUT
%   FRC           struct            Forced response curve with one field F1,
%                                   F2, ... for each forcing amplitude in
%                                   fRed


%% This function converts the normal form style reduced dynamics on the SSM
% to fixed point problem to compute the periodic response (FRC) for
% internally resonant systems according to the following article
% 

%% preliminary options setup
coordinates = 'polar';
parName = 'freq';
initialSolver = 'forward';% 'fsolve';
parRange = frequencySpan;
contOptions = cocoOptions();
set(contOptions, 'PtMX', 1000, 'h_max', 0.5); %, 'ItMX',20,'TOL',1e-3);
nPar = 200;
nCycle = 500;
sampStyle = 'cocoBD'; % 'uniform'
m = numel(mFreqs);
nt = 100;
%% Collecting autonomous reduced dynamics coefficients
[lambda, R0, maxOrder] = collect_reduced_dynamics_coeffs(RDInfoF);

check_spectrum_and_internal_resonance(real(lambda),imag(lambda),mFreqs,0.1);

[beta,kappa] = extract_beta_kappa(R0, maxOrder, mFreqs);

%% Collecting non-autonomous reduced dynamics ceofficients
forcingCoeffs = zeros(2*m,1);
forcingCoeffs(1:2:end) = RDInfoF.conjugateDynamics.forcingVectors;

forcingFourier.harmonics = 1;
forcingFourier.coeffs = forcingCoeffs;
[fdata,data_dir,Omega0] = create_reduced_dynamics_data(beta,kappa,lambda,mFreqs,forcingFourier,maxOrder,resonantModes);

%% Performing continuation
FRCdata = struct();
for ii = 1:length(epsilon)
p0 = [Omega0(1); epsilon(ii)];
[FRC, runid] = coco_compute_FRC(runid,fdata, p0, parName,coordinates,nCycle,nPar, parRange,initialSolver,contOptions, sampStyle);
FRC = FRC_polar2complex(FRC,nt,mFreqs);
FRCdata = computePhysicalAmplitude(ii,IMInfoF, RDInfoF, FRCdata, FRC, epsilon(ii), amplitudeFunction);
end
FRC = FRCdata; varargout{1} = FRC; 
fdir = fullfile(data_dir,runid,'SSMep.mat');
save(fdir, 'FRC');
end

function [lambda, R0, maxOrder] = collect_reduced_dynamics_coeffs(RDInfoF)
RDcoeffs = RDInfoF.conjugateDynamics.coefficients;
RDexponents = RDInfoF.conjugateDynamics.exponents;
[RDexponentsCC,RDcoeffsCC] = transform_complex_conjugate(RDexponents, RDcoeffs);

orders = sum(RDexponentsCC,2);
maxOrder = max(orders);
R0 = repmat(struct('coeffs',[],'ind',[]),1,maxOrder);

for j = 1:maxOrder
    [I] = find(orders==j);
    if ~isempty(I)
        R0(j).ind = RDexponentsCC(I,:);
        R0(j).coeffs = RDcoeffsCC(:,I);
    end
end
lambda = diag(R0(1).coeffs);
end

function check_spectrum_and_internal_resonance(lambdaRe,lambdaIm,mFreqs,tol)
% check spectrum and internal resonance

flags1 = abs(lambdaRe(1:2:end-1)-lambdaRe(2:2:end))<1e-6*abs(lambdaRe(1:2:end-1)); % same real parts
flags1 = all(flags1);
flags2 = abs(lambdaIm(1:2:end-1)+lambdaIm(2:2:end))<1e-6*abs(lambdaIm(1:2:end-1)); % opposite imag parts
flags2 = all(flags2);
freqs  = lambdaIm(1:2:end-1);
freqso = freqs - dot(freqs,mFreqs(:))*mFreqs(:)/sum(mFreqs.^2);
flags3 = norm(freqso)<tol*norm(freqs);
assert(flags1, 'Real parts do not follow complex conjugate relation');
assert(flags2, 'Imaginary parts do not follow complex conjugate relation');
assert(flags3, 'Internal resonnace is not detected for given master subspace');
end

function [beta,kappa] = extract_beta_kappa(R_0,order,mFreqs)
% check reduced dynamics (factor out exp(i*Omega*t) part
m  = numel(mFreqs);
Em = eye(m);
beta  = cell(m,1); % coefficients - each cell corresponds to one mode
kappa = cell(m,1); % exponants
for k = 2:order
    R = R_0(k);
    coeffs = R.coeffs;
    ind = R.ind;
    if ~isempty(coeffs)
        for i=1:m
            betai = coeffs(2*i-1,:);
            [~,ki,betai] = find(betai);
            kappai = ind(ki,:);
            % check resonance condition
            l = kappai(:,1:2:end-1);
            j = kappai(:,2:2:end);
            nk = numel(ki);
            rm = repmat(mFreqs(:)',[nk,1]);
            flagi = dot(l-j-repmat(Em(i,:),[nk,1]),rm,2);
            assert(all(flagi==0), 'Reduced dynamics is not consisent with desired IRs');
            % assemble terms
            beta{i}  = [beta{i} betai];
            kappa{i} = [kappa{i}; kappai];
        end
    end
end
end

function [fdata,data_dir,Omega0] = create_reduced_dynamics_data(beta,kappa,lambda,mFreqs,forcingFourier,order,resonantModes)
% CREATE_REDUCED_DYNAMICS_DATA This function create a data structure for
% the vector field of reduced dynamics (leading-order nonautonomous
% approximation)
forcingHarmonics = forcingFourier.harmonics;
forcingCoeffs = forcingFourier.coeffs;
nHarmonics = length(forcingHarmonics);
iNonauto = []; % indices for resonant happens
rNonauto = []; % value of leading order contribution
Omega0 = [];
for k=1:nHarmonics
    idm = find(mFreqs(:)==forcingHarmonics(k)); % idm could be vector if there are two frequencies that are the same
    Omega0 = [Omega0 imag(lambda(2*idm(1)-1))];
    r = forcingCoeffs(2*idm-1,k);
    iNonauto = [iNonauto; idm];
    rNonauto = [rNonauto; r];
end

fdata = struct();
fdata.beta  = beta;
fdata.kappa = kappa;
fdata.lamdRe = real(lambda(1:2:end-1));
fdata.lamdIm = imag(lambda(1:2:end-1));
fdata.mFreqs = mFreqs;
fdata.iNonauto = iNonauto;
fdata.rNonauto = rNonauto;

data_dir = fullfile(pwd,'data');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end
fdata.order = order;
fdata.modes = resonantModes;

end

function [exponentsCC, coeffsCC] = transform_complex_conjugate(exponents, coeffs)

[nMonomials, nVars] = size(exponents);
exponentsCC = sparse(2*nMonomials, nVars);

exponentsCC(1:2:end,1:2:end) = exponents(:,1:nVars/2);
exponentsCC(1:2:end,2:2:end) = exponents(:,(nVars/2+1):end);
exponentsCC(2:2:end,1:2:end) = exponents(:,(nVars/2+1):end);
exponentsCC(2:2:end,2:2:end) = exponents(:,1:nVars/2);

[nEqs, nMonomials] = size(coeffs);

coeffsCC = sparse(2*nEqs, 2*nMonomials);
coeffsCC(1:2:end,1:2:end) = coeffs;
coeffsCC(2:2:end,2:2:end) = conj(coeffs);
end


function [FRC, runid] = coco_compute_FRC(oid, fdata, p0, parName,coordinates,nCycle,nPar, parRange,initialSolver,contOptions,sampStyle)

%% basic setup
switch parName
    case 'freq'
        isomega = true;
    case 'amp'
        isomega = false;
    otherwise
        error('Continuation parameter should be freq or amp');
end

ispolar = strcmp(coordinates, 'polar');
fdata.ispolar = ispolar;
m = numel(fdata.mFreqs);
if ispolar
    odefun = @(z,p) ode_2mDSSM_polar(z,p,fdata);
    z0 = 0.1*ones(2*m,1);
else
    odefun = @(z,p) ode_2mDSSM_cartesian(z,p,fdata);
    z0 = zeros(2*m,1);
end
%% construct initial solution
funcs  = {odefun};
prob = coco_prob();
prob = cocoSet(contOptions, prob);

z0 = get_initial_sol(z0,p0,initialSolver,odefun,nCycle,ispolar);
%% continuation of reduced dynamics w.r.t. parName
% call ep-toolbox
prob = ode_isol2ep(prob, '', funcs{:}, z0, {'om' 'eps'}, p0);
% define monitor functions to state variables
[prob, args1, args2] = monitor_states(prob, ispolar, m);

if strcmp(sampStyle, 'uniform')
    if isomega
        omSamp = linspace(parRange(1),parRange(2), nPar);
        prob   = coco_add_event(prob, 'UZ', 'om', omSamp);
    else
        epSamp = linspace(parRange(1),parRange(2), nPar);
        prob   = coco_add_event(prob, 'UZ', 'eps', epSamp);
    end
end

runid = coco_get_id(oid, 'ep');

fprintf('\n Run=''%s'': Continue equilibria along primary branch.\n', ...
    runid);

if isomega
    cont_args = [{'om'},args1(:)' ,args2(:)',{'eps'}];
else
    cont_args = [{'eps'},args1(:)' ,args2(:)',{'om'}];
end

coco(prob, runid, [], 1, cont_args, parRange);

%% extract results of reduced dynamics at sampled frequencies
FRC = ep_reduced_results(runid,sampStyle,ispolar,isomega,args1,args2);
end

function z0 = get_initial_sol(z0,p0,initialSolver,odefun,nCycle,ispolar)
switch initialSolver
    case 'fsolve'
        % fsolve to approximate equilibrium
        fsolveOptions = optimoptions('fsolve','MaxFunctionEvaluations',100000,...
            'MaxIterations',1000000,'FunctionTolerance', 1e-10,...
            'StepTolerance', 1e-8, 'OptimalityTolerance', 1e-10);
        z0 = fsolve(@(z) odefun(z,p0),z0,fsolveOptions);
    case 'forward'
        % forward simulation to approach equilibirum
        tspan = [0 nCycle*2*pi/p0(1)]; %nCycle
        odefw = @(t,z,p) odefun(z,p);
        opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
        [~,y0] = ode45(@(t,y) odefw(t,y,p0), tspan, z0, opts);
        [~,y] = ode45(@(t,y) odefw(t,y,p0), [0 2*pi/p0(1)], y0(end,:));
        [~, warnId] = lastwarn;
        
        if any(isnan(y(:))) || strcmp(warnId,'MATLAB:ode45:IntegrationTolNotMet')
            warning('Reduced dynamics with IRs in polar form diverges with [0.1 0.1 0.1 0.1]');
        else
            z0 = y(end,:)';
        end
    otherwise
        error('initialSolver must be set to forward or fsolve')
end

if ispolar % regularize initial solution if it is in polar form
    z0(2:2:end) = mod(z0(2:2:end),2*pi); % phase angles in [0,2pi]
    m = numel(z0)/2;
    for k=1:m
        if z0(2*k-1)<0
            z0(2*k-1) = -z0(2*k-1);      % positive amplitudes
            z0(2*k) = z0(2*k)+pi;
        end
    end
end
end

function FRC = FRC_polar2complex(FRC, nt, mFreqs)
Omega = FRC.om;
rho = FRC.rho;
stab = FRC.st;
eps = FRC.ep;
Z = FRC.z;
theta = FRC.th;
SNidx = FRC.SNidx;
HBidx = FRC.HBidx;

numPts = numel(Omega);
nSN = numel(SNidx);
isSN = logical(sparse(SNidx, ones(nSN,1), ones(nSN,1),numPts,1));
nHB = numel(HBidx);
isHB = logical(sparse(HBidx, ones(nHB,1), ones(nHB,1),numPts,1));

FRC   = cell(numPts,1);
phi = linspace(0,2*pi,nt)/min(mFreqs);
m = numel(mFreqs);
for i=1:numPts
    % phase space coordinates
    state = Z(i,:);
    z = transpose(state).*exp(1i*(mFreqs'.*phi));
    p = zeros(2*m,nt);
    p(1:2:end-1,:) = z;
    p(2:2:end,:)   = conj(z);
    
    FRC{i} = struct('rho', rho(i,:), 'stability', stab(i), 'Omega', Omega(i) ,...
        'epsilon', eps(i), 'p', p, 'th', theta(i,:),'isSN',isSN(i),'isHB',isHB(i));
end

FRC = cat(1,FRC{:});
end

% Compute physical amplitude 
function FRCout = computePhysicalAmplitude(iAmp,IMInfo, RDInfo, FRCout, FRC, epsilon, amplitudeFunction)
Omega = [FRC.Omega];
nSteps = length(Omega);
rho = [FRC.rho]; rho = reshape(rho,length(rho)/nSteps,nSteps);
stab = [FRC.stability];
psi = [FRC.th]; psi = reshape(psi,length(psi)/nSteps,nSteps);

% Define intervals for approximation of the physical periodic response
SSMFunction = IMInfo.parametrization.map;
T = RDInfo.transformation.map;
dimSSM = 2*size(rho,1);
z = zeros(dimSSM,nSteps); 
dimAmp = size(amplitudeFunction(SSMFunction(0,...
                       T(0,zeros(dimSSM,1),Omega(1),0,0),Omega(1),0,0)),1);
if dimAmp == 1
    u = zeros(1,nSteps); uPhase = zeros(1,nSteps); 
else
    u = zeros(1,nSteps,dimAmp); uPhase = zeros(1,nSteps,dimAmp);     
end
for ind = 1:nSteps
    iOmega = Omega(ind);
    zEval = [FRC(ind).p]; zEval = [zEval(1:2:end,:); zEval(2:2:end,:)];
    normTimeEval = linspace(0,1,size(zEval,2));
    phiEval =  normTimeEval*2*pi;
    cPhi = cos(phiEval); sPhi = sin(phiEval); 
    timeEval = phiEval/iOmega;
    etaEval = T(timeEval,zEval,iOmega,epsilon,0);
%     y = SSMFunction(etaEval);
    y = SSMFunction(timeEval,etaEval,iOmega,epsilon,0);
    yAmplitudeFunction = amplitudeFunction(y);
    z(:,ind) = zEval(:,1);
    if dimAmp == 1 
    u(ind) = max(abs(yAmplitudeFunction));
    uPhase(ind) = phase( normTimeEval(2) * sum( yAmplitudeFunction.*(cPhi-1i*sPhi ) ));
    else
    u(1,ind,:) = max(abs(yAmplitudeFunction),[],2);   
    uPhase(1,ind,:) = phase( normTimeEval(2) * sum( yAmplitudeFunction.*(cPhi-1i*sPhi ),2 ));
    end
end

FRCout.(['F' num2str(iAmp)]) = struct('Freq',Omega,'Amp',...
        u,'Nf_Amp',rho,'Nf_Phs2',psi,'Nf_Phs',uPhase,'Stab',stab,'z_IC',z);
end




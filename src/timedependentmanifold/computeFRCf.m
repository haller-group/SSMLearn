function FRC = computeFRCf(IMInfoF, RDInfoF, fRed, amplitudeFunction, varargin)
%   FRC = computeFRCf(IMInfo, RDInfo, fRed, amplitudeFunction)
%   Compute the forced response curves on a 2m SSM using coco for each
%   normal form forcing amplitude in fRed.
%
%   INPUT
%   IMInfo             struct            Time-periodic manifold.
%   RDInfo             struct            Time-periodic reduced dynamics.
%   f_red              (1 x nAmp)        Normal form forcing constants.
%   amplitudeFunction  function handle   Map from y to (signed) scalar. Can
%                                        be used to predict forced response
%                                        for a given component or quantity.
%   varargin           scalar            continuations options for coco
%
%   OUTPUT
%   FRC           struct            Forced response curve with one field F1,
%                                   F2, ... for each forcing amplitude in
%                                   fRed

% Initialize functions & checks
SSMFunction = IMInfoF.parametrization.map;
T = RDInfoF.transformation.map;
n = length(RDInfoF.eigenvaluesLinPartFlow);
if strcmp(RDInfoF.conjugacyStyle,'normalform') == 0
    error('Error. For using this function, the reduced dynamics must be known in normal form.')
end
if strcmp(RDInfoF.dynamicsType,'flow') == 0
    error('Error. For using this function, the reduced dynamics must be a flow (continuous time dynamical system).')
end
if mod(n,2)>0
    error('Error. For using this function, all the modes must be oscillatory ones.')
else
    ndof = n / 2;
end

% Setup reduced dynamics in SSMTool-compatible form ...


% Define intervals for approximation of the physical periodic response
Npers = 1;
normTimeEval = linspace(0,1*Npers,250*Npers+1); normTimeEval(end) = [];
normTimeEval = normTimeEval(1:end-1);
phiEval =  normTimeEval*2*pi;
cPhi = cos(phiEval); sPhi = sin(phiEval);

% Iterate for each forcing amplitude value (epsilon)
for iAmp = 1:length(fRed)
    
    % Compute responses via coco
    rho = []; Omega = []; psi = []; stab = [];
    % Here we know the initial conditions of nCont periodic orbits from numerical continuation, namely:
    % rho: ndof x nCont matrix (amplitudes)
    % Omega: 1 x nCont vector (frequencies)
    % psi: ndof x nCont matrix (phase shifts: psi_j = theta_j - Omega t)
    % stab: 1 x nCont vector (logical vector, 1 is stable)
    
    nCont = size(rho,2);
    u = zeros(1,nCont); uPhase = u;
    % Compute physical amplitudes
    for iCont = 1:nCont
        iOmega = Omega(iCont);
        timeEval = phiEval/iOmega;
        thetaEval = phiEval + psi(:,iCont);
        zEval = rho(:,iCont).*exp(1i*thetaEval);
        etaEval = T(timeEval,[zEval; conj(zEval)],iOmega,fRed(iAmp),0);
        y = SSMFunction(timeEval,etaEval,iOmega,fRed(iAmp),0);
        yAmplitudeFunction = amplitudeFunction(y);
        u(iCont) = max(abs(yAmplitudeFunction));
        % Integral phase definition over physical observable
        uPhase(iCont) = angle( normTimeEval(2)* sum( yAmplitudeFunction.*(cPhi-1i*sPhi ) ));
    end
    % Save in the final struct
    FRC.(['F' num2str(iAmp)]) = struct('Freq',Omega,'Amp',...
        u,'Nf_Amp',rho,'Nf_Phs2',psi,'Nf_Phs',uPhase,'Stab',stab);
end
function [IMInfoF,RDInfoF] = forcedSSMROM(IMInfo,RDInfo,varargin)
%   [IMInfoF,RDInfoF] = forcedSSMROM(IMInfo,RDInfo,varargin)
%   Derive a forced m-dim. SSM model. The default is a time-periodic SSM 
%   model, but any finite number of (positive) frequencies can be inserted, 
%   standing hence for quasi-periodic forcing. 
%   The output model does not take of resonances between the linearized
%   frequencies and forcing ones: the reduced dynamics features forcing in
%   every mode and the normal form retains only the +\Omega t rotation,
%   while the normal form coordinate change has the other rotation.
%
%   INPUT
%   IMInfo               struct      Manifold.
%   RDInfo               struct      Reduced dynamics.
%   optional
%   nForcingFrequencies  integer     Number of forcing frequencies 
%   forcingVectors       matrix      ( n_observables x nForcingFrequencies)
%                                    are the forcing vectors in the
%                                    physical space
%   We                   matrix      ( m x n_observables ) projection to
%                                    the tangent space of the SSM at the
%                                    orgin
%   Wo, Lo, Vo           matrices    projection, dynamics and eigenspace of
%                                    the outer directions. They could be
%                                    partial, i.e. for a limited number of
%                                    modes.
%   OUTPUT
%   IMInfoF              struct      Time-periodic manifold.
%   RDInfoF              struct      Time-periodic reduced dynamics.
%
%   For full state space observables with knowledge of additional (outer)
%   modes, the time periodic parametrization also consider their effect on 
%   the time periodic response. Instead, if only SSM-related (inner) modes
%   are available (e.g. generic observables), then forcing is only assumed
%   along modal directions.

% Define inputs
p = inputParser;
addOptional(p, 'nForcingFrequencies', 1, @(x) isnumeric(x) && isscalar(x));
addOptional(p, 'forcingVectors', []);
addOptional(p, 'We', []);
addOptional(p, 'Lo', []);
addOptional(p, 'Vo', []);
addOptional(p, 'Wo', []);
addOptional(p, 'outdof',[]);
parse(p, varargin{:});
outdof = p.Results.outdof;
% Eventual errors
if strcmp(RDInfo.dynamicsType,'flow')~=1
    error('Forced SSM ROMs are only available for flows.'); end
if sum(imag(RDInfo.eigenvaluesLinPartFlow)==0)>0
    error('Forced SSM ROMs are only available for complex conjugate eigenvalues.'); end
ndof = length(RDInfo.eigenvaluesLinPartFlow)/2;
% Adjust forcing vectors
if isempty(p.Results.forcingVectors) == 1
    flagNonModalForcing = 0;
    % Forcing is assumed on each mode of each frequency
    modalForcing = 1i*ones(ndof,p.Results.nForcingFrequencies);
    disp('Forced SSM reduced-order model assumes external forcing only along the tangent (modal) subspace at the origin. ')
else
    flagNonModalForcing = 1; % Assume to observe the full state space and know the linear part
    ROMForcing = p.Results.We*p.Results.forcingVectors;
    if strcmp(RDInfo.conjugacyStyle,'default')==1
        modalForcing = RDInfo.eigenvectorsLinPart \ ROMForcing; 
    else
        modalForcing = RDInfo.inverseTransformation.lintransf * ROMForcing;
    end
    modalForcing = 0.5*modalForcing(1:ndof,:);
    % Setup 
    if isempty(p.Results.Lo) == 1
        error('Outer linearized dynamics needed.');
    else
        if isvector(p.Results.Lo) == 1
            Vo = p.Results.Lo;
            Wo = p.Results.Wo;
            if size(p.Results.Lo,2) == length(p.Results.Lo)
                Lo = transpose(p.Results.Lo); 
            else
                Lo = p.Results.Lo; 
            end
        else
            if isdiag(p.Results.Lo) == 1
                Lo = diag(p.Results.Lo);
                Vo = p.Results.Lo;
                Wo = p.Results.Wo;
            else
                [So,Lo] = eig(p.Results.Lo);
                Lo = diag(Lo);
                Vo = p.Results.Vo*So;
                Wo = So\p.Results.Wo;
            end
        end
    end
    fullNonModalForcing = 0.5*Wo*p.Results.forcingVectors;
end
% Parametrization correction
IMInfoF = IMInfo;
IMInfoF.parametrization.numberForcingFrequencies = p.Results.nForcingFrequencies;
autParam = IMInfoF.parametrization.map;
autParamOut = IMInfoF.parametrization.mapOut;
if flagNonModalForcing == 1
    IMInfoF.parametrization.forcingVectors = p.Results.forcingVectors;
    IMInfoF.parametrization.forcingVectorsModal = IMInfoF.parametrization.tangentSpaceAtOrigin*ROMForcing;
    IMInfoF.parametrization.forcingVectorsNonModal = real(Vo*fullNonModalForcing);
    forcingNonAutParam= @(t,fFreqs,fAmpls,fPhs) -real(Vo * ( ( ((fAmpls.*exp(1i*fPhs)).*fullNonModalForcing)./(Lo-1i*fFreqs) ) * exp(+1i*(transpose(fFreqs).*t)) + ( ((fAmpls.*exp(-1i*fPhs)).*fullNonModalForcing)./(Lo+1i*fFreqs) ) * exp(-1i*(transpose(fFreqs).*t))  ) ) ;
    forcingNonAutParamOut= @(t,fFreqs,fAmpls,fPhs) -real(Vo(outdof,:) * ( ( ((fAmpls.*exp(1i*fPhs)).*fullNonModalForcing)./(Lo-1i*fFreqs) ) * exp(+1i*(transpose(fFreqs).*t)) + ( ((fAmpls.*exp(-1i*fPhs)).*fullNonModalForcing)./(Lo+1i*fFreqs) ) * exp(-1i*(transpose(fFreqs).*t))  ) ) ;
    nonAutParam = @(t,q,fFreqs,fAmpls,fPhs) autParam(q) + forcingNonAutParam(t,fFreqs,fAmpls,fPhs);
    nonAutParamOut = @(t,q,fFreqs,fAmpls,fPhs) autParamOut(q) + forcingNonAutParamOut(t,fFreqs,fAmpls,fPhs);
else
    forcingNonAutParam = @(t,fFreqs,fAmpls,fPhs) 0 ;
    nonAutParam = @(t,q,fFreqs,fAmpls,fPhs) autParam(q);
    nonAutParamOut = @(t,q,fFreqs,fAmpls,fPhs) autParamOut(q);
end
IMInfoF.parametrization.map = nonAutParam;
IMInfoF.parametrization.mapOut = nonAutParamOut;
IMInfoF.parametrization.forcingPart = forcingNonAutParam;

% Dynamics correction
RDInfoF = RDInfo;
RDInfoF.numberForcingFrequencies = p.Results.nForcingFrequencies;
autRedDyn = RDInfoF.reducedDynamics.map;
autConDyn = RDInfoF.conjugateDynamics.map;
autInvTransf = RDInfoF.inverseTransformation.map;
autTransf = RDInfoF.transformation.map;
eigLinROM = RDInfo.eigenvaluesLinPartFlow(1:ndof);
forcingNonAutRedDyn = @(t,fFreqs,fAmpls,fPhs)  real(RDInfoF.eigenvectorsLinPart * transformationComplexConj( ((fAmpls.*exp(1i*fPhs)).*modalForcing) * exp(+1i*(transpose(fFreqs).*t)) + ((fAmpls.*exp(-1i*fPhs)).*modalForcing) * exp(-1i*(transpose(fFreqs).*t))  ) );
forcingVecNonAutRedDyn = real(RDInfoF.eigenvectorsLinPart * transformationComplexConj( modalForcing * 2 ) );
switch RDInfoF.conjugacyStyle
    case 'normalform'
        forcingNonAutConDyn = @(t,fFreqs,fAmpls,fPhs)  transformationComplexConj( ((fAmpls.*exp(1i*fPhs)).*modalForcing) * exp(+1i*(transpose(fFreqs).*t) ) );
        forcingNonAutInvTransf = @(t,fFreqs,fAmpls,fPhs) transformationComplexConj( (((fAmpls.*exp(-1i*fPhs)).*modalForcing)./(eigLinROM+1i*fFreqs) ) * exp(-1i*(transpose(fFreqs).*t) ) );
        forcingNonAutTransf = @(t,fFreqs,fAmpls,fPhs) real(RDInfoF.eigenvectorsLinPart* transformationComplexConj(-(((fAmpls.*exp(-1i*fPhs)).*modalForcing)./(eigLinROM+1i*fFreqs) ) * exp(-1i*(transpose(fFreqs).*t) ) ) );
        forcingVecNonAutConDyn = modalForcing;
    case 'modal'
        forcingNonAutConDyn = @(t,fFreqs,fAmpls,fPhs)  transformationComplexConj( ((fAmpls.*exp(1i*fPhs)).*modalForcing) * exp(+1i*(transpose(fFreqs).*t)) + ((fAmpls.*exp(-1i*fPhs)).*modalForcing) * exp(-1i*(transpose(fFreqs).*t))  );
        forcingNonAutInvTransf = @(t,fFreqs,fAmpls,fPhs) 0;
        forcingNonAutTransf = @(t,fFreqs,fAmpls,fPhs) 0;
        forcingVecNonAutConDyn = modalForcing;
    otherwise
        forcingNonAutConDyn = @(t,fFreqs,fAmpls,fPhs) forcingNonAutRedDyn(t,fFreqs,fAmpls,fPhs);
        forcingNonAutInvTransf = @(t,fFreqs,fAmpls,fPhs) 0;
        forcingNonAutTransf = @(t,fFreqs,fAmpls,fPhs) 0;
        forcingVecNonAutConDyn = forcingVecNonAutRedDyn;
end
% Assignments
RDInfoF.reducedDynamics.forcingVectors = forcingVecNonAutRedDyn;
RDInfoF.reducedDynamics.map = @(t,x,fFreqs,fAmpls,fPhs) real(autRedDyn(x)) + forcingNonAutRedDyn(t,fFreqs,fAmpls,fPhs);
RDInfoF.reducedDynamics.forcingPart = @(t,fFreqs,fAmpls,fPhs) forcingNonAutRedDyn(t,fFreqs,fAmpls,fPhs);
RDInfoF.conjugateDynamics.forcingVectors = forcingVecNonAutConDyn;
RDInfoF.conjugateDynamics.map = @(t,z,fFreqs,fAmpls,fPhs) autConDyn(z) + forcingNonAutConDyn(t,fFreqs,fAmpls,fPhs);
RDInfoF.conjugateDynamics.forcingPart = @(t,fFreqs,fAmpls,fPhs) forcingNonAutConDyn(t,fFreqs,fAmpls,fPhs);
RDInfoF.inverseTransformation.map = @(t,x,fFreqs,fAmpls,fPhs) autInvTransf(x) + forcingNonAutInvTransf(t,fFreqs,fAmpls,fPhs);
RDInfoF.inverseTransformation.forcingPart = @(t,fFreqs,fAmpls,fPhs) forcingNonAutInvTransf(t,fFreqs,fAmpls,fPhs);
RDInfoF.transformation.map = @(t,z,fFreqs,fAmpls,fPhs) autTransf(z) + forcingNonAutTransf(t,fFreqs,fAmpls,fPhs);
RDInfoF.transformation.forcingPart = @(t,fFreqs,fAmpls,fPhs) forcingNonAutTransf(t,fFreqs,fAmpls,fPhs);
if flagNonModalForcing == 0
    IMInfoF.parametrization.forcingVectors = IMInfoF.parametrization.tangentSpaceAtOrigin*forcingVecNonAutRedDyn;
    IMInfoF.parametrization.forcingPart = @(t,fFreqs,fAmpls,fPhs) IMInfoF.parametrization.tangentSpaceAtOrigin*forcingNonAutRedDyn(t,fFreqs,fAmpls,fPhs);    
end
end




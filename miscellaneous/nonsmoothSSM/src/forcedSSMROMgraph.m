function [IMInfoF,RDInfoF] = forcedSSMROMgraph(IMInfo,RDInfo,forcingVectors,WoC,LoC,VoC,W0,V0)
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

% Parametrization correction
IMInfoF = IMInfo;
IMInfoF.parametrization.numberForcingFrequencies = 1;
autParam = IMInfoF.parametrization.map;

fullNonModalForcing = WoC*(V0*W0-eye(size(V0,1)))*forcingVectors;
IMInfoF.parametrization.forcingVectors = fullNonModalForcing;
forcingNonAutParam= @(t,fFreqs,fAmpls,fPhs) 0.5*real(VoC * ( ( ((fAmpls.*exp(1i*fPhs)).*fullNonModalForcing)./(LoC-1i*fFreqs) ) * exp(+1i*(transpose(fFreqs).*t)) + ( ((fAmpls.*exp(-1i*fPhs)).*fullNonModalForcing)./(LoC+1i*fFreqs) ) * exp(-1i*(transpose(fFreqs).*t))  ) ) ;
nonAutParam = @(t,q,fFreqs,fAmpls,fPhs) autParam(q) + forcingNonAutParam(t,fFreqs,fAmpls,fPhs);
IMInfoF.parametrization.map = nonAutParam;
IMInfoF.parametrization.forcingPart = forcingNonAutParam;

% Dynamics correction
RDInfoF = RDInfo;
RDInfoF.numberForcingFrequencies = 1;
autRedDyn = RDInfoF.reducedDynamics.map;
autConDyn = RDInfoF.conjugateDynamics.map;
autInvTransf = RDInfoF.inverseTransformation.map;
autTransf = RDInfoF.transformation.map;
ROMForcing = W0*forcingVectors;
forcingNonAutRedDyn = @(t,fFreqs,fAmpls,fPhs) real(((fAmpls.*exp(1i*fPhs)).*ROMForcing) * exp(+1i*(transpose(fFreqs).*t))) ;
forcingVecNonAutRedDyn = ROMForcing;
forcingNonAutConDyn = @(t,fFreqs,fAmpls,fPhs) forcingNonAutRedDyn(t,fFreqs,fAmpls,fPhs);
forcingNonAutInvTransf = @(t,fFreqs,fAmpls,fPhs) 0;
forcingNonAutTransf = @(t,fFreqs,fAmpls,fPhs) 0;
forcingVecNonAutConDyn = forcingVecNonAutRedDyn;

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
end
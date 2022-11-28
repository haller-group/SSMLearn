function FRC = analyticalFRC(IMInfoF, RDInfoF, fRed, amplitudeFunction, varargin)
%   FRC = analyticalFRCf(IMInfoF, RDInfoF, fRed, amplitudeFunction)
%   Compute the forced response curves on a 2D SSM analytically for each
%   normal form forcing amplitude in fRed.
%
%   INPUT
%   IMInfoF            struct            Time-periodic manifold.
%   RDInfoF            struct            Time-periodic reduced dynamics.
%   fRed               (1 x nAmp)        Normal form forcing constants.
%   amplitudeFunction  function handle   Map from y to (signed) scalar. Can
%                                        be used to predict forced response
%                                        for a given component or quantity.
%   varargin           scalar            optional max amplitude rho at
%                                        which the model is trusted (avoid
%                                        spurious isolas)
%   OUTPUT
%   FRC           struct            Forced response curve with one field F1,
%                                   F2, ... for each forcing amplitude in
%                                   fRed

% The conventions used are
% \dot{z} = (\alpha(\rho) + i*(omega(\rho) - \Omega) ) z + i*f*e^{i\Omega t}
% z = \rho * e^{i\theta t}
% f = fscale * fRed * exp(1i fphase)
% psi = \theta - \Omega t - fphase [negative normal form phase psi is indeed a lag for positive damping].

% Initialize functions
SSMFunction = IMInfoF.parametrization.map;
damp = RDInfoF.conjugateDynamics.damping;
freq = RDInfoF.conjugateDynamics.frequency;
T = RDInfoF.transformation.map;
fscale = abs(RDInfoF.conjugateDynamics.forcingVectors);
zfTemp = RDInfoF.conjugateDynamics.forcingVectors / 1i / fscale ;
fphase = atan2(imag(zfTemp), real(zfTemp));
coeffs = RDInfoF.conjugateDynamics.coefficients;
exponents = RDInfoF.conjugateDynamics.exponents;

% Define useful terms
nPolyTerms = 2 * length(coeffs); orders = sum(exponents, 2);
% Coefficients of the damping & frequency
dampcoeffs = zeros(1,nPolyTerms); dampcoeffs(orders) = real(coeffs);
freqcoeffs = zeros(1,nPolyTerms);freqcoeffs(orders) = imag(coeffs);
% Setup derivatives for the jacobians
if nPolyTerms == 2
    DomegaDrho = @(r) zeros(size(r));
    DaDrho = @(r) dampcoeffs(1)*ones(size(r));
else
    Dfreqcoeffs = freqcoeffs.*[0:(nPolyTerms-1)];
    Dfreqcoeffs = Dfreqcoeffs(3:2:end);
    Dfreqexpnts = orders-2; Dfreqexpnts = Dfreqexpnts(2:end);
    DomegaDrho = @(r) Dfreqcoeffs*(r.^Dfreqexpnts);
    Dacoeffs = dampcoeffs.*[1:nPolyTerms]; Dacoeffs = Dacoeffs(1:2:end);
    Daexpnts = Dfreqexpnts+1;
    DaDrho = @(r) Dacoeffs(1) + Dacoeffs(2:end)*(r.^Daexpnts);
end
% Define intervals for approximation of the physical periodic response
Npers = 1;
normTimeEval = linspace(0,1*Npers,250*Npers+1); normTimeEval(end) = [];
normTimeEval = normTimeEval(1:end-1);
phiEval =  normTimeEval*2*pi;
cPhi = cos(phiEval); sPhi = sin(phiEval);

% Spacing for the rho parametrization of the frequency response
nEvalInt = 301;
chebyNodes = (cos([(nEvalInt-1):-1:0]/(nEvalInt-1)*pi)+1)/2;

% Iterate for each forcing
for iAmp = 1:length(fRed)
    % Compute roots of the sqrt argument ( (f/r)^2 - a^2(r) )
    rhoSol = roots([fliplr(dampcoeffs(1:end-1)),-(fRed(iAmp)*fscale)]);
    rhoSol = sort(abs(rhoSol(imag(rhoSol)==0))); % eliminate spurios
    rhoMin = fminsearch(@(r) (freq(r)-sqrt((fRed(iAmp)*fscale./r).^2-damp(r).^2 )).^2,(fRed(iAmp)*fscale)/abs(coeffs(1)));
    rhoSol = [rhoMin; rhoSol];
    if isempty(varargin) == 0
        rhoMax = varargin{:};
        rhoSol = [rhoSol(1:sum(rhoSol<rhoMax)); rhoMax];
    end
    % Compute the response when the sqrt argument is postivie
    rho = []; Omega = []; psi = []; stab = []; u = []; uPhase = []; z = [];
    for iTip = 1:(length(rhoSol)-1)
        % Chebyshev nodal spacing (denser at the sides)
        rhoInt = rhoSol(iTip) + chebyNodes*(rhoSol(iTip+1)-rhoSol(iTip));
        % Eventually one can use
        % rhoInt = linspace(rhoSol(iTip),rhoSol(iTip+1),nEvalInt);
        % rhoInt = logspace(log10(rhoSol(iTip)),log10(rhoSol(iTip+1)),nEvalInt);
        rhoE = rhoInt(round(nEvalInt/2));
        if ((fRed(iAmp)*fscale/rhoE)^2-damp(rhoE)^2)>=0
            % Evaluate fixed points
            rhoDamp = damp(rhoInt); rhoFreq = freq(rhoInt);
            rhoSqrt = real(sqrt( (fRed(iAmp)*fscale./rhoInt).^2-damp(rhoInt).^2 ));
            rhoPhs = atan(rhoDamp./rhoSqrt) ;
            rho = [rho; rhoInt; rhoInt];
            Omega = [Omega; rhoFreq-rhoSqrt; rhoFreq+rhoSqrt];
            psi = [psi; rhoPhs+pi; -rhoPhs];
            % Evaluate their stability
            rhoDa = DaDrho(rhoInt); rhoDomega = DomegaDrho(rhoInt);
            traceJ = rhoDa + rhoDamp;
            detJP = rhoDa.*rhoDamp - (+rhoInt.*rhoSqrt).*(rhoDomega+(-rhoSqrt)./rhoInt);
            detJM = rhoDa.*rhoDamp - (-rhoInt.*rhoSqrt).*(rhoDomega+(+rhoSqrt)./rhoInt);
            stab = [stab; (((traceJ<0)+(detJM>0))==2); (((traceJ<0)+(detJP>0))==2)];
            % Obtain physical amplitude and phase of the first harmonic
            for iSol = [1 0]
                uTemp = zeros(1,nEvalInt); uPhaseTemp = zeros(1,nEvalInt);
                zTemp = zeros(1,nEvalInt);
                for iRho = 1:nEvalInt
                    iOmega = Omega(end-iSol,iRho);
                    timeEval = phiEval/iOmega;
                    thetaEval = phiEval + fphase + psi(end-iSol,iRho);
                    zEval = rho(end-iSol,iRho)*exp(1i*thetaEval);
                    etaEval = T(timeEval,[zEval; conj(zEval)],iOmega,fRed(iAmp),0);
                    y = SSMFunction(timeEval,etaEval,iOmega,fRed(iAmp),0);
                    yAmplitudeFunction = amplitudeFunction(y);
                    zTemp(iRho) = zEval(1);
                    uTemp(iRho) = max(abs(yAmplitudeFunction));
                    zfTemp = normTimeEval(2)* sum( yAmplitudeFunction.*(cPhi-1i*sPhi ) );
                    uPhaseTemp(iRho) = atan2(imag(zfTemp), real(zfTemp));
                end
                z = [z; zTemp]; u = [u; uTemp]; uPhase = [uPhase; uPhaseTemp];
            end
        end
    end
    
    FRC.(['F' num2str(iAmp)]) = struct('Freq',Omega,'Amp',...
        u,'Nf_Amp',rho,'Nf_Phs2',psi,'Nf_Phs',uPhase,'Stab',stab,'z_IC',z);
end
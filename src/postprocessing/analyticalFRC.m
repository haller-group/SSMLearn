function FRC = analyticalFRC(IMInfo, RDInfo, fRed, amplitudeFunction,varargin)
%   FRC = analyticalFRC(IMInfo, RDInfo, fRed, amplitudeFunction)
%   Compute the forced response curves on a 2D SSM analytically for each 
%   normal form forcing amplitude in fRed. 
%   
%   INPUT
%   IMInfo             struct            Manifold.
%   RDInfo             struct            Reduced dynamics.
%   f_red              (1 x nAmp)        Normal form forcing constants.
%   amplitudeFunction  function handle   Map from y to scalar. Can be
%                                        used to predict forced response 
%                                        for a given component or quantity.
%   varargin           scalar            optional max amplitude at which
%                                        the model is trusted (avoid 
%                                        spurious isolas)
%   OUTPUT
%   FRC           struct            Forced response curve with one field F1, 
%                                   F2, ... for each forcing amplitude in
%                                   fRed

SSMFunction = IMInfo.parametrization.map;
damp = RDInfo.conjugateDynamics.damping;
freq = RDInfo.conjugateDynamics.frequency;
T = RDInfo.transformation.map;
coeffs = RDInfo.conjugateDynamics.coefficients;
exponents = RDInfo.conjugateDynamics.exponents;

nPolyTerms = 2 * length(coeffs);
dampcoeffs = zeros(1,nPolyTerms);
dampcoeffs(nPolyTerms-sum(exponents, 2)) = real(coeffs);
% freqcoeffs = zeros(1,nPolyTerms);
% freqcoeffs(1+nPolyTerms-sum(exponents, 2)) = imag(coeffs);

for iAmp = 1:length(fRed)
    rhoTip = roots([dampcoeffs(1:end-1),-fRed(iAmp)]);
    rhoTip = abs(rhoTip(imag(rhoTip)==0));
    rhoTip = [min(rhoTip)*0.003; sort(rhoTip)];
    rho = [];
    for iTip = 1:2:length(rhoTip)
        rho = [rho; linspace(rhoTip(iTip), rhoTip(iTip+1), 1000)];
    end
    rho = [rho, -fliplr(rho)];
    Omega = real(freq(rho) + -1./rho.*sqrt(fRed(iAmp)^2-(rho.*damp(rho)).^2));
    rho = abs(rho);
    
    u = zeros(size(rho));
    eitheta = exp(1i*linspace(-pi,pi,51)); eitheta(end) = [];
    for iPart = 1:size(rho,1)
        for iRho = 1:size(rho,2)
            y = SSMFunction(T([rho(iPart,iRho)*eitheta; rho(iPart,iRho)*conj(eitheta)]));
            u(iPart,iRho) = max(abs(amplitudeFunction(y)));
        end
    end
    
    psi = acos((Omega - freq(rho)).*rho./fRed(iAmp));
    
    eps = 1e-10;
    stab = zeros(size(rho));
    for iPart = 1:size(rho,1)
        rhoPart = rho(iPart,:);
        dadrho = ((rhoPart+eps).*damp(rhoPart+eps) - (rhoPart-eps).*damp(rhoPart-eps)) / eps * 0.5;
        dbdrho = (freq(rhoPart+eps) - freq(rhoPart-eps)) / eps * 0.5;
        J = zeros(2,2,size(rho,2));
        J(1,1,:) = dadrho;
        J(2,1,:) = dbdrho - (Omega(iPart,:) - freq(rhoPart))./rhoPart;
        J(1,2,:) = (Omega(iPart,:) - freq(rhoPart)).*rhoPart;
        J(2,2,:) = (rhoPart.*damp(rhoPart))./rhoPart;
        for iRho = 1:length(rhoPart)
            eigenvalues = eig(det(J(:,:,iRho)));
            stab(iPart,iRho) = all(eigenvalues > 0);
        end
    end
    if isempty(varargin) == 0
       ampMax = varargin{:};
       pos = u>ampMax;
       Omega(pos) = []; u(pos) = []; rho(pos) = []; 
       psi(pos) = []; stab(pos) = [];
    end
    FRC.(['F' num2str(iAmp)]) = struct('Freq',Omega,'Amp',...
                u,'Nf_Amp',rho,'Nf_Phs',psi,'Stab',stab);
end
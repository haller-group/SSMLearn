function [fRed, OmegaFm, OmegaFp, uF, psiF] = miniComputeFRC(yCal, Omega, BBCInfo, Hmap, iHmap, Tmap, iTmap, amplitudeFunction)

%% calibration
ndof = size(iHmap(yCal), 1)*0.5;
nAmp = size(yCal,2);
fRed = zeros(nAmp,ndof);
for iAmp = 1:nAmp
    zCal = iTmap(iHmap(yCal(:,iAmp)));
    rhoCal = abs(zCal(1:ndof,:));
    fRed(iAmp,:) = sqrt(rhoCal.^2.*(BBCInfo.freq(rhoCal)-Omega(iAmp)).^2 ...
        + rhoCal.^2.*(BBCInfo.damp(rhoCal)).^2);
end

%% FRC computation
OmegaFm = @(rho, fR) BBCInfo.freq(rho) - 1./rho.*sqrt(fR.^2-(rho.*BBCInfo.damp(rho)).^2);
OmegaFp = @(rho, fR) BBCInfo.freq(rho) + 1./rho.*sqrt(fR.^2-(rho.*BBCInfo.damp(rho)).^2);

includesConj = size(zCal, 1)-1;
uF = @(rho) getamp(rho, Hmap, Tmap, amplitudeFunction, includesConj);

psiF = @(rho, fR, Omega) acos((Omega - BBCInfo.freq(rho)).*rho./fR);
end
function u = getamp(rho, Hmap, Tmap, amplitudeFunction, includesConj)
    if includesConj
        eitheta = exp([1i;-1i].*linspace(-pi,pi,151)); eitheta(:,end) = [];
    else
        eitheta = exp(1i*linspace(-pi,pi,151)); eitheta(:,end) = [];
    end
    u = zeros(1,size(rho,2));
    for ie = 1:size(eitheta,2)
        u(ie,:) = amplitudeFunction(Hmap(Tmap(rho.*eitheta(:,ie))));
    end
    u = max(abs(u), [], 1);
end
    
%     eps = 1e-10;
%     stab = zeros(size(rho));
%     for iPart = 1:size(rho,1)
%         rhoPart = rho(iPart,:);
%         dadrho = ((rhoPart+eps).*damp(rhoPart+eps) - (rhoPart-eps).*damp(rhoPart-eps)) / eps * 0.5;
%         dbdrho = (freq(rhoPart+eps) - freq(rhoPart-eps)) / eps * 0.5;
%         J = zeros(2,2,size(rho,2));
%         J(1,1,:) = dadrho;
%         J(2,1,:) = dbdrho - (Omega(iPart,:) - freq(rhoPart))./rhoPart;
%         J(1,2,:) = (Omega(iPart,:) - freq(rhoPart)).*rhoPart;
%         J(2,2,:) = (rhoPart.*damp(rhoPart))./rhoPart;
%         for iRho = 1:length(rhoPart)
%             eigenvalues = eig(det(J(:,:,iRho)));
%             stab(iPart,iRho) = all(eigenvalues > 0);
%         end
%     end
%     
%     FRC.(['F' num2str(iAmp)]) = struct('Freq',Omega,'Amp',...
%                 u,'Nf_Amp',rho,'Nf_Phs',psi,'Stab',stab);

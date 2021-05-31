function FRC = computeFRC(f_red, damp, freq, SSMFunction, T, yObservable)

options = optimoptions('fsolve', 'MaxFunctionEvaluations', 10000);
for iAmp = 1:length(f_red)
    rhoTip = abs(fsolve(@(rho) 1e5*(f_red(iAmp)-(rho*damp(rho))), f_red(iAmp), options));
    rho = logspace(log10(rhoTip*0.003), log10(rhoTip), 1000);
    rho = [rho, -fliplr(rho)];
    Omega = real(freq(rho) + -1./rho.*sqrt(f_red(iAmp)^2-(rho.*damp(rho)).^2));
    rho = abs(rho);
    eitheta = exp(1i*linspace(-pi,pi,51)); eitheta(end) = [];
    for iRho = 1:length(rho)
        y = SSMFunction(T([rho(iRho)*eitheta;rho(iRho)*conj(eitheta)]));
        u(iRho) = max(abs(yObservable(y)));
    end
    
    psi = acos((Omega - freq(rho)).*rho./f_red(iAmp));
    
    eps = 1e-10;
    dadrho = ((rho+eps).*damp(rho+eps) - (rho-eps).*damp(rho-eps)) / eps * 0.5;
    dbdrho = (freq(rho+eps) - freq(rho-eps)) / eps * 0.5;
    J(1,1,:) = dadrho;
    J(2,1,:) = dbdrho - (Omega - freq(rho))./rho;
    J(1,2,:) = (Omega - freq(rho)).*rho;
    J(2,2,:) = (rho.*damp(rho))./rho;
    stab = zeros(size(rho));
    for iRho = 1:length(rho)
        eigenvalues = eig(det(J(:,:,iRho)));
        stab(iRho) = all(eigenvalues > 0);
    end
    
    FRC.(['F' num2str(iAmp)]) = struct('Freq',Omega,'Amp',...
                u,'Nf_Amp',rho,'Nf_Phs',psi,'Stab',stab);
end
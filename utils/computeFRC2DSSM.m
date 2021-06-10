function FRC = computeFRC(f_red, damp, freq, SSMFunction, T, yObservable, N_info)

plen = 2*length(N_info.coeff);
dampcoeffs = zeros(1,plen); freqcoeffs = zeros(1,plen);
dampcoeffs(plen-sum(N_info.exponents,2)) = real(N_info.coeff);
freqcoeffs(1+plen-sum(N_info.exponents,2)) = imag(N_info.coeff);

for iAmp = 1:length(f_red)
    
    rhoTip = roots([dampcoeffs(1:end-1),-f_red(iAmp)]);
    rhoTip = abs(rhoTip(imag(rhoTip)==0));
    rhoTip = [min(rhoTip)*0.003; sort(rhoTip)];
    rho = [];
    for iTip = 1:2:length(rhoTip)
        rho = [rho, linspace(rhoTip(iTip), rhoTip(iTip+1), 1000)];
    end
    rho = [rho, -fliplr(rho)];
    J = zeros(2,2,length(rho)); u = 0*rho;
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
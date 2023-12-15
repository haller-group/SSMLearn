function FRC_data = SSMToolFRCcustom(M, C, K, fnl, fext, fcoeffs, outdof, omegaRange, orders,varargin)
% Use SSMTool to extract forced response curves from a dynamical system. 
% fext is a vector containing the periodic forcing
% outdof is the index of the plotted degree of freedom
% orders is a scalar/vector of the degree of Taylor expansions to be plotted

[DS, S, ~] = getSSM(M, C, K, fnl, 1:2, orders);

epsilon = 1;
kappas = [-1; 1];

set(S.FRCOptions, 'method', 'continuation ep', 'z0', 1e-4*[1; 1]) % 'level set'
set(S.Options, 'reltol', 1, 'IRtol', 0.02, 'notation', 'multiindex', 'contribNonAuto', true)
set(S.FRCOptions, 'nt', 2^7, 'nRho', 200, 'nPar', 200, 'nPsi', 200, 'rhoScale', 2)
set(S.FRCOptions, 'outdof', outdof,'nCycle',300)
if isempty(varargin) == 0 % Increase continuation steps
    contSteps = varargin{:};
    set(S.contOptions,'PtMX',contSteps); set(S.contOptions,'ItMX',100);
end

FRC_data = struct();
for ii = 1:length(fcoeffs)
    coeffs = [fext fext]/2*fcoeffs(ii);
    DS.add_forcing(coeffs, kappas, epsilon);
    FRC = S.extract_FRC('freq', omegaRange, orders); 
    FRC_data.(['F' num2str(ii)]) = struct('Freq',[FRC.Omega],'Amp',...
   [FRC.Aout],'Nf_Amp',[FRC.rho],'Nf_Phs',[FRC.th],'Stab',[FRC.stability]);
end

set(S.FRCOptions, 'method', 'level set') % 'level set'
for ii = 2:length(fcoeffs)
    coeffs = [fext fext]/2*fcoeffs(ii);
    DS.add_forcing(coeffs, kappas, epsilon);
    FRC = S.extract_FRC('freq', omegaRange, orders);
    FreqO = FRC_data.(['F' num2str(ii)]).Freq; 
    Nf_ampO = FRC_data.(['F' num2str(ii)]).Nf_Amp; 
    Nf_PhsO = FRC_data.(['F' num2str(ii)]).Nf_Phs;
    AmpO = FRC_data.(['F' num2str(ii)]).Amp; 
    StabO = FRC_data.(['F' num2str(ii)]).Stab;
    Freq = [FRC.Omega]; Nf_amp = [FRC.rho]; Nf_Phs = [FRC.psi];
    Amp = [FRC.Aout]; Stab = [FRC.stability];
    if ii > 1
    idxUns = find(Stab==0); 
    [unsMinAmp,pos] = min(Amp(idxUns)); unsMinFreq = 660;%Freq(idxUns(pos)); 
    idxStab2 = find((Freq>unsMinFreq).*(Amp<unsMinAmp));
    newPos = [idxUns idxStab2];
    newPos = [setdiff(1:length(Stab),newPos) newPos];
    newPos = [idxStab2];
    else
       newPos = 1:length(Stab);
    end
    FRC_data.(['F' num2str(ii)]) = struct('Freq',[FreqO FreqO(end) Freq(newPos)],'Amp',...
    [AmpO AmpO(end) Amp(newPos)],'Nf_Amp',[Nf_ampO Nf_ampO(end) Nf_amp(newPos)],'Nf_Phs',[Nf_PhsO Nf_PhsO(end) Nf_Phs(newPos)],'Stab',[StabO Stab(newPos(1)) Stab(newPos)]);
end

end

function FRC_data = SSMToolFRC(M, C, K, fnl, fext, fcoeffs, outdof, omegaRange, orders,varargin)
% Use SSMTool to extract forced response curves from a dynamical system. 
% fext is a vector containing the periodic forcing
% outdof is the index of the plotted degree of freedom
% orders is a scalar/vector of the degree of Taylor expansions to be plotted

SSMOrder = 9; 
[DS, S, ~] = getSSM(M, C, K, fnl, 1:2, SSMOrder);

epsilon = 1;
kappas = [-1; 1];

set(S.FRCOptions, 'method', 'continuation ep', 'z0', 1e-4*[1; 1]) % 'level set' 
set(S.Options, 'reltol', 1, 'IRtol', 0.02, 'notation', 'multiindex', 'contribNonAuto', true)
set(S.FRCOptions, 'nt', 2^7, 'nRho', 200, 'nPar', 200, 'nPsi', 200, 'rhoScale', 2)
set(S.FRCOptions, 'outdof', outdof,'nCycle',300)
if isempty(varargin) == 0 % Increase continuation steps
    contSteps = varargin{:};
    set(S.contOptions,'PtMX',contSteps);
end

FRC_data = struct();
for ii = 1:length(fcoeffs)
    coeffs = [fext fext]/2*fcoeffs(ii);
    DS.add_forcing(coeffs, kappas, epsilon);
    FRC = S.extract_FRC('freq', omegaRange, orders);
    FRC_data.(['F' num2str(ii)]) = struct('Freq',[FRC.Omega],'Amp',...
   [FRC.Aout],'Nf_Amp',[FRC.rho],'Nf_Phs',[FRC.th],'Stab',[FRC.stability]);
end
end

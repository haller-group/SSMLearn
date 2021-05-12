function FRC_data = getFRC_full(M, C, K, fnl, f_vec, f_coeffs, ...
     outdof, omegaRange, orders)
% Use SSMTool to extract forced response curves from a dynamical system. 
% fext is a vector containing the periodic forcing
% outdof is the index of the plotted degree of freedom
% orders is a scalar/vector of the degree of Taylor expansions to be plotted

% Setup
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
masterModes = [1,2]; 
S.choose_E(masterModes);

epsilon = 1;
kappas = [-1; 1];
set(S.Options, 'reltol', 1,'IRtol',0.02,'notation', 'multiindex','contribNonAuto',true)
set(S.FRCOptions, 'nt', 2^7, 'nRho', 200, 'nPar', 200, 'nPsi', 100, 'rhoScale', 2 )
set(S.FRCOptions, 'method','continuation ep', 'z0', 1e-4*[1; 1]) 
set(S.FRCOptions, 'outdof',outdof)
set(S.contOptions, 'PtMX', 500, 'h_min', 1e-3)

FRC_data = struct();
for ii = 1:length(f_coeffs)
    coeffs = [f_vec f_vec]/2*f_coeffs(ii);
    DS.add_forcing(coeffs, kappas, epsilon);
    FRC = S.extract_FRC('freq', omegaRange, orders);
    FRC_data.(['F' num2str(ii)]) = struct('Freq',[FRC.Omega],'Amp',...
   [FRC.Aout],'Nf_Amp',[FRC.rho],'Nf_Phs',[FRC.th],'Stab',[FRC.stability]);
end
end

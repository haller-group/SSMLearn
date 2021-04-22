function [FRC] = getFRC(M, C, K, fnl, fext, outdof, orders)
% Use SSMTool to extract forced response curves from a dynamical system. 
% fext is a vector containing the periodic forcing
% outdof is the index of the plotted degree of freedom
% orders is a scalar/vector of the degree of Taylor expansions to be plotted

SSMOrder = 17;
[DS, S, ~] = getSSM(M, C, K, fnl, 2, SSMOrder);

epsilon = 1;
kappas = [-1; 1];
coeffs = [fext fext]/2;
DS.add_forcing(coeffs, kappas, epsilon);

set(S.Options, 'reltol', 1, 'IRtol', 0.02, 'notation', 'multiindex', 'contribNonAuto', true)
set(S.FRCOptions, 'nt', 2^7, 'nRho', 200, 'nPar', 200, 'nPsi', 100, 'rhoScale', 2)
set(S.FRCOptions, 'method', 'continuation ep', 'z0', 1e-4*[1; 1]) % 'level set' 
set(S.FRCOptions, 'outdof', outdof)

omega0 = imag(S.E.spectrum(1));
omegaRange = omega0*[0.9 1.5];

FRC = S.extract_FRC('freq', omegaRange, orders);

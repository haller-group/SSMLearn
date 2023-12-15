function [SSMparametrization, SSMdynamics, DS, S] = get2DSSMprops(M, C, K, fnl, mfdOrder)
% [IC, mfd, DS, S] = getSSMIC(M, C, K, fnl, nTraj, ICRadius, SSMDim, seed)
% [IC, mfd, DS, S] = getSSMIC(M, C, K, fnl, nTraj, ICRadius, SSMDim, seed, order, E)
% Computes an SSM and returns points on it at a given distance from the origin.
% The distance metric is in the parametrization coordinates on the manifold.
%
% INPUT
% M         (nxn)            Mass matrix.
% C         (nxn)            Damping matrix.
% K         (nxn)            Stiffness matrix.
% fnl       {nx(2n)^2, nx(2n)^3, ...}  Nonlinearities tensor.
% nTraj     int              Number of initial conditions to generate.
% ICRadius  real >=0         Distance in parametrization coordinates from
%                               the origin to initial conditions.
% SSMDim    int >0           Manifold dimension.
% seed      real or rand     Set to a value for reproducible results, set
%                               to rand for random IC spacing
% order     int >0           Optional manifold order, default 9.
% E         int array        Optional eigenvalue indices, default 1:SSMDim.
%                               Overrides SSMDim when passed.
%
% OUTPUT
% IC        (nxnTraj)        Initial conditions matrix.
% mfd       {1xorder}        Coefficients and indices of SSM expansion.
% DS        DynamicalSystem  Obtained from SSMTool.
% S         SSM              Obtained from SSMTool.


E = 1:2;
DS = DynamicalSystem();
set(DS, 'M', M, 'C', C, 'K', K, 'fnl', fnl);
set(DS.Options, 'Emax', 5, 'Nmax', 10, 'notation', 'multiindex')
S = SSM(DS);
set(S.Options, 'reltol', 0.1, 'notation', 'multiindex')
S.choose_E(E);
[mfd,dyn] = S.compute_whisker(mfdOrder);

SSMparametrization = @(z) real(getManifoldPoint(mfd, z));
gamma = compute_gamma(dyn);
lambda = S.E.spectrum(1);
SSMdynamics = @(z) getManifoldVF(abs(z), gamma, lambda);
end

function z = getManifoldVF(rho,gamma,lambda)
a = rho * real(lambda);
b = rho * imag(lambda);

for j = 1:length(gamma)
    a = a + real(gamma(j))* rho.^(2*j+1);
    b = b + imag(gamma(j)) * (rho.^(2*j+1));
end
z = [a+1i*b; a-1i*b];
end

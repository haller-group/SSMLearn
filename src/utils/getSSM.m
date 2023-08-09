function [DS, S, mfd] = getSSM(M, C, K, fnl, E, varargin)
% [DS, S, mfd] = getSSM(M, C, K, fnl, E, order)
% Generates a DynamicalSystem and an SSM for given tensors.
%   
% INPUT
% M         (nxn)            Mass matrix.
% C         (nxn)            Damping matrix.
% K         (nxn)            Stiffness matrix.
% fnl       {nx(2n)^2, nx(2n)^3, ...}  Nonlinearities tensor.
% E         int array        Eigenvalue indices of the tangent space, e.g.
%                               1:SSMDim or [3,4,7,8]
% order     int >0           Optional manifold order, default 9.
%
% OUTPUT
% DS        DynamicalSystem  Obtained from SSMTool.
% S         SSM              Obtained from SSMTool.
% mfd       {1xorder}        Coefficients and indices of SSM expansion.

mfdOrder = 9;
if ~isempty(varargin)
    mfdOrder = varargin{:};
end

DS = DynamicalSystem();
set(DS, 'M', M, 'C', C, 'K', K, 'fnl', fnl);
set(DS.Options, 'Emax', 5, 'Nmax', 10, 'notation', 'tensor')
S = SSM(DS);
set(S.Options, 'reltol', 0.1, 'notation', 'tensor')
S.choose_E(E);
mfd = S.compute_whisker(mfdOrder);
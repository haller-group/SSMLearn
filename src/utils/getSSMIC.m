function [IC, mfd, DS, S] = getSSMIC(M, C, K, fnl, nTraj, ICRadius, SSMDim, seed, varargin)
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


mfdOrder = 9;
E = [1:SSMDim];
if ~isempty(varargin)
    mfdOrder = varargin{1};
    if numel(varargin) > 1
        E = varargin{2};
    end
end

[DS, S, mfd] = getSSM(M, C, K, fnl, E, mfdOrder);

N = size(M,1);
IC = zeros(2*N, nTraj);
z = ICRadius * pickPointsOnHypersphere(nTraj, length(E), seed);
z = repelem(z(1:2:end,:) + 1j*z(2:2:end,:),2,1);
z(2:2:end) = conj(z(2:2:end));
for iTraj = 1:nTraj
    IC(:,iTraj) = real(getManifoldPoint(mfd, z(:,iTraj)));
end
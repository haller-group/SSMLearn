function [F, lambda, V, G, DG] = functionFromTensors(M, C, K, fnl, varargin)
% functionFromTensors(M, C, K, fnl)
% functionFromTensors(M, C, K, fnl, fExt, Omega)
% Creates an anonymous function for the evolution of a mechanical system.
%
% INPUT
% M      (nxn)            Mass matrix.
% C      (nxn)            Damping matrix.
% K      (nxn)            Stiffness matrix.
% fnl    {nx2nx2n, nx2nx2nx2n, ...}   Cell array containing nonlinear
%                         polynomial coefficent tensors.
% fExt   (nx1)            Optional: Forcing tensor.
% Omega  scalar>=0        Optional: Forcing frequency.
%
% OUTPUT
% F      function handle     F(t, q) = \dot{q}
% lambda (2nx1)              eigenvalues of the linearized system at 0
% V      (2nx2n)             eigenvectors of the linearized system at 0
%
% EXAMPLES
% Get a function handle for the evolution of the system
% m*qddot1 + c*qdot1 + k*q1 = 0
% m*qddot2 + c*qdot2 + k*q2 = 0:
%   functionFromTensors(m*eye(2), c*eye(2), k*eye(2), {})
%
% Get a function handle for the evolution of the system
% m*qddot1 + k*q1 + b*q1dot*q2^2 = 0
% m*qddot2 + k*q2 + a*q1^2 = 0:
%   fnl = {sptensor(zeros(2, 2*2, 2*2)), sptensor(zeros(2, 2*2, 2*2, 2*2))}
%   fnl{1}(2,1,1) = a;        % q1^2
%   fnl{2}(1,2+1,2,2) = b;    % q1dot*q2^2
%   functionFromTensors(m*eye(2), zeros(2), k*eye(2), fnl)
%
% Mechanical tensors should be passed as compatible with SSMTool. See
%   DynamicalSystem for more info.
% Integrating F is considerably faster than the tensor-formulated system.

forced = 0;
if ~isempty(varargin)
    forced = 1;
    fExt = varargin{1};
    Omega = varargin{2};
end

n = size(M,1);
q = sym('q', [2*n 1]);
qdot = sym(zeros(n,1));
for iT = 1:length(fnl)
    Fi = fnl{iT};
    idx = find(Fi);
    for i = 1:size(idx,1)
        qdot(idx(i,1)) = qdot(idx(i,1)) + Fi(idx(i,:))*prod(q(idx(i,2:end)));
    end
end
dqdot = sym(zeros(n,2*n));
for ii = 1:2*n
   dqdot(:,ii) = diff(qdot,q(ii));
end
g = matlabFunction(qdot, 'vars', {q});
Dg = matlabFunction(dqdot, 'vars', {q});

Minv = inv(M);
A = [zeros(n), eye(n);
    -Minv*K,  -Minv*C];
G = @(x) [zeros(n,1);
         -Minv*g(x)];
DG = @(x) [zeros(n,2*n);
         -Minv*Dg(x)];     
if forced
    H = @(t,x) [zeros(n,1);
        -Minv*fExt*cos(Omega*t)];
    F = @(t,x) A*x + G(x) + H(t,x);
else
    F = @(t,x) A*x + G(x);
end

[V, D] = eigSorted(full(A));
lambda = diag(D);
end
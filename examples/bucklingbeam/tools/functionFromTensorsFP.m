function [F, lambda, V] = functionFromTensorsFP(M, C, K, fnl, fExt, unstable_fp, varargin)
% functionFromTensorsFP(M, C, K, fnl, loadvector, unstable_fp);
% functionFromTensorsFP(M, C, K, fnl, loadvector, unstable_fp, A_t);
% Move the system to one fixed point, and creates an anonymous function for
% the evolution of a mechanical system.
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

recal = 1;
if ~isempty(varargin)
    recal = 0; A_t = varargin{1};
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
g = matlabFunction(qdot, 'vars', {q});

if recal
    A_t = zeros(n,n);
    for i = 1:n
        for j = 1:n
            temp = diff(qdot(i), q(j));
            A_t(i,j) = subs(temp, q , unstable_fp);
        end
    end
end

Minv = inv(M);
A = [zeros(n), eye(n);
    -Minv*K,  -Minv*C];
A_tilde = [zeros(n), eye(n);
    -Minv*(K+A_t),  -Minv*C];

G = @(x) [zeros(n,1);
         -Minv*g(x + unstable_fp)];
H = @(t,x) [zeros(n,1);
        -Minv*( fExt + K*unstable_fp(1:n) )];
F = @(t,x) A*x + G(x) + H(t,x);


[V, D] = eigSorted(full(A_tilde));
lambda = diag(D);
end
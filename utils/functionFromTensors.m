function [F, lambda, V] = functionFromTensors(M, C, K, fnl, varargin)

forced = 0;
if ~isempty(varargin)
    forced = 1;
    fext = varargin{1};
    omega = varargin{2};
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

Minv = inv(M);
A = [zeros(n), eye(n);
    -Minv*K,     -Minv*C];
G = @(x) [zeros(n,1);
         -Minv*g(x)];
     
if forced
    H = @(t,x) [zeros(n,1);
        -Minv*fext*cos(omega*t)];
    F = @(t,x) A*x + G(x) + H(t,x);
else
    F = @(t,x) A*x + G(x);
end

[V,D] = eig(full(A));
[lambda,pos] = sort(diag(D)); 
V = V(:,pos);
end
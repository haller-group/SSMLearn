function [F, lambda] = tensorfunction2(M, C, K, fnl)

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
f = matlabFunction(qdot, 'vars', {q});

Minv = inv(M);
A = [zeros(n), eye(n);
    -Minv*K,     -Minv*C];
G = @(x) [zeros(n,1);
         -Minv*f(x)];
F = @(t,x) A*x + G(x);

lambda = sort(eig(full(A)));
function [F, lambda] = tensorfunction(M, C, K, fnl)

n = size(M,1);

q = sym('q', [2*n 1]);
qdot = sym(zeros(n,1));
t = sym('t');
for iT = 1:length(fnl)
    Fi = fnl{iT};
    idx = find(Fi);
    for i = 1:size(idx,1)
        qdot(idx(i,1)) = qdot(idx(i,1)) + Fi(idx(i,:))*prod(q(idx(i,2:end)));
    end
end

A = [zeros(n), eye(n);
    -M\K,     -M\C];
G = [zeros(n,1);
    -M\qdot];
F = matlabFunction(A*q + G, 'vars', {t q});
lambda = sort(eig(full(A)));
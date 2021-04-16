function f = tensorfunction2(fnl)

n = size(fnl{1},1);
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
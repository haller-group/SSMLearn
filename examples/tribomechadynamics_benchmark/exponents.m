function exps = exponents(d, k)
    B = repmat({0:max(k)}, 1, d);
    A = combvec(B{:}).';
    exps = A(ismember(sum(A, 2), k),:);
    [~,ind] = sort(sum(exps, 2));
    exps = exps(ind,:);
end
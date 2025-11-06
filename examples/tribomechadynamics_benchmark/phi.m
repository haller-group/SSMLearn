function u = phi(xi, r) % return monomials
    x = reshape(xi, 1, size(xi, 1), []);
    exps = exponents(size(xi, 1), r);
    u = reshape(prod(x.^exps, 2), size(exps, 1), []);
end
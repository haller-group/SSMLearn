function p = getManifoldPoint(mfd, z)
    % Computes coordinates on a manifold given reduced coordinates z
    %
    % INPUT
    % mfd   (1 x md) cell array computed from compute_whisker in
    %       SSMTool-2.0
    % z     (2 x 1) reduced coordinates in complex form
    
    n = size(mfd{1}.coeffs,1);
    p = zeros(n,1);
    for i=1:length(mfd)
        p = p + mfd{i}.coeffs*phi(z,i);
    end

end


function S=Sigma(d,k)
    B = repmat({0:k},1,d);
    A = combvec(B{:})';
    S = A(sum(A,2)==k,:);
end


function u = phi(xi,r)
    % order 2 and up
    u = [];
    get_prods = @(x,rr) prod(x'.^Sigma(size(x,1),rr), 2);
    for l = r
        u = [u; get_prods(xi,l)];
    end
end
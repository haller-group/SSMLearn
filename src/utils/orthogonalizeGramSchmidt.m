function Vort = orthogonalizeGramSchmidt(V)
% Given the space col(V), the function outputs an orthonormal basis for
% this space in the columns of Vort.

Vort = V; Vort(:,1) = V(:,1)/norm(V(:,1)); [~,k] = size(V);
if k>1
    for ii = 2:k
        Vort(:,ii) = Vort(:,ii) - Vort(:,1:(ii-1)) * ( ...
            transpose(Vort(:,1:(ii-1))) * Vort(:,ii) );
        Vort(:,ii) = Vort(:,ii)/norm(Vort(:,ii));
    end
end
end
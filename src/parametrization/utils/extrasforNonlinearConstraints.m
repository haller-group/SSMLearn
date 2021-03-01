function idx_info = extrasforNonlinearConstraints(n,k,phi_dim)
% Helper for IMparametrization

% Computations of indexes for the computation of the nonlinear constraints

% For the computation of V'*V=I
idx_info = cell(1,4);
idx_VV = []; for ii = 1:k; idx_VV = [idx_VV [ii:k]+(ii-1)*k ]; end
idx_info{1} = idx_VV;
% For the computation of the constraint derivatives
% Derivatives for V'*V = I
Der_index = transpose(multivariante_exponents(k,2));
[Ipos,Jpos] = ind2sub(size(Der_index),find(Der_index));
[~,pos] = sort(Ipos); Jpos = Jpos(pos);
rows_v = repmat(transpose(1:n),k,k)+repmat([0:k-1]*n,n*k,1);
rows_Dceq_VV = rows_v(:);
cols_v = repmat(transpose(Jpos),n,1);
cols_Dceq_VV = cols_v(:);
idx_info{2} = reshape(ones(k)+diag(ones(1,k)),1,k*k);
% Derivatives for V'*H = 0
% Derivatives w.r.t. V
rows_Dceq_VH_v = repmat(transpose(1:n*k),phi_dim,1);
cols_h = repmat(1:k*phi_dim,n,1);
cols_Dceq_VH_v = cols_h(:) + k*(k+1)/2;
% Derivatives w.r.t. H
rows_h = repmat(n*k+[1:n],phi_dim,1)+transpose([0:phi_dim-1]*n);
rows_Dceq_VH_h = reshape(repmat(transpose(rows_h),k,1),n*k*phi_dim,1);
cols_Dceq_VH_h = cols_Dceq_VH_v;
rows_Dceq = [rows_Dceq_VV; rows_Dceq_VH_v; rows_Dceq_VH_h];
cols_Dceq = [cols_Dceq_VV; cols_Dceq_VH_v; cols_Dceq_VH_h];
idx_info{3} = rows_Dceq;
idx_info{4} = cols_Dceq;
end


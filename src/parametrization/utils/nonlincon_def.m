function [c,ceq,Dc,Dceq] = nonlincon_def(z,n,k,phi_dim,idx_info)
% Helper for IMparametrization

% Definition of nonlinear equality constraints and their gradients
% These are the vectorized version of V'*V = I and V'*H = 0

V = reshape(z(1:n*k),n,k);
H = reshape(z(1+n*k:end),n,phi_dim);
% We have no inequality constraints
c = []; Dc = [];
% Equality constraints
Cv = transpose(V)*V-eye(k,k); Cv = Cv(:); Ch = transpose(V)*H;
ceq = [Cv([idx_info{1}]) ; Ch(:)];
if nargout>2
    vals_v = repmat(V,1,k).*[idx_info{2}];
    vals_Dceq_VV = vals_v(:);
    vals_Dceq_VH_v = reshape(repmat(H,k,1),n*k*phi_dim,1);
    vals_Dceq_VH_h = repmat(V(:),phi_dim,1);
    vals_Dceq = [vals_Dceq_VV; vals_Dceq_VH_v; vals_Dceq_VH_h];
    Dceq = sparse([idx_info{3}],[idx_info{4}],vals_Dceq,length(z),...
        length(ceq));
end
end
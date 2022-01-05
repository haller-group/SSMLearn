function [phi,Exp_mat,D_phi_info] = multivariatePolynomialSelection(k,MinDeg,MaxDeg,idx)
% Definition of the polynomial map phi as a Matlab function. The map phi has
% monomials in k variables from order MinDeg to order MaxDeg for the row
% indexes of the exponents matrix idx. These exponents are collected in the
% Exp_mat_full, while D_phi_info is a cell array containing the information
% (function & indexes) for computing the derivatives of phi.

Exp_mat_full = []; 
for ii = MinDeg : MaxDeg; Exp_mat_full = [Exp_mat_full; multivariateExponents(k,ii)]; end
u = sym('u',[1 k]);
Exp_mat = Exp_mat_full(idx,:);
phi_sym = prod(u.^Exp_mat,2);
phi = matlabFunction(phi_sym,'Vars', {transpose(u)} );
if nargout > 2 % Define derivatives
    D_phi_info = cell(k,2);
    for ii = 1:k
        Di_phi_sym = diff(phi_sym,['u' num2str(ii)]);
        ind_nz = find(Di_phi_sym); ind_z = find(Di_phi_sym==0);
        D_phi_info{ii,1} = matlabFunction(Di_phi_sym(ind_nz),'Vars', {transpose(u)} );
        [~,ind_sort] = sort([ind_nz; ind_z]);
        D_phi_info{ii,2} = ind_sort;
    end
end
end

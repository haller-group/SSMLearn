function [phi,Exp_mat,D_phi_info] = multivariateFractionalPolynomial(k,MinDeg,MaxDeg, varargin)
% Definition of the polynomial map phi as a Matlab function. The map phi has
% monomials in k variables from order MinDeg (which can be only 1 or 2) to
% order MaxDeg. The exponents are collected in the matrix Exp_mat, while
% D_phi_info is a cell array containing the information (coefficients &
% indexes) for computing the derivatives of phi based on phi.

Exp_mat = []; 
for ii = 1 : MaxDeg; Exp_mat = [Exp_mat; multivariateExponents(k,ii)]; end

% fractional terms
if nargin==4
    if k>1
        disp("Fractional terms for multivariate expressions not implemented.");
    else 
        if ~isempty(cell2mat(varargin(1))) % default to regular polynomials if empty array is supplied
            Exp_mat = cell2mat(varargin(1));
        end
    end
end


u = sym('u',[1 k]);
phi_sym = prod(u.^Exp_mat,2);
if nargout > 2 % Define derivatives
    D_phi_info = cell(k,2);
    for ii = 1:k
        li = zeros(1,k); li(ii) = 1;
        Di_phi_sym = Exp_mat(:,ii).*prod(u.^(Exp_mat-li),2); % Alternative to diff(phi_sym,['u' num2str(ii)]);
        if MinDeg == 1
        ind_nz = find(Di_phi_sym); ind_1 = ind_nz(1);
        ind_var = ind_nz(2:end); ind_0 = find(Di_phi_sym==0);
        else
        ind_var = find(Di_phi_sym(k+1:end)); ind_1 = [];
        ind_0 = find(Di_phi_sym(k+1:end)==0);
        end
        D_phi_info{ii,1} = Exp_mat(ind_var,ii);
        [~,D_phi_info{ii,2}] = sort([ind_var; ind_1; ind_0]);
        % sum([Di_phi_sym(k*(MinDeg-1)+ind_var)-Exp_mat(k*(MinDeg-1)+ind_var,ii).*phi_sym(1:length(ind_var))])
    end
end
if MinDeg == 1
    phi = matlabFunction(phi_sym,'Vars', {transpose(u)} );    
else
    phi = matlabFunction(phi_sym(k+1:end),'Vars', {transpose(u)} );
    Exp_mat = Exp_mat(k+1:end,:); % Exclude linear terms
end
end

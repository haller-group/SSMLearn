function [Aine,bine,Aequ,bequ] = alignmentLinearConstraint(z_dim,V_ref)
% Helper for IMparametrization

% Linear equality and inequality constraints that ensure an alignment that
% makes the optimization routine to have a unique solution: given a refe-
% rence representation V_ref, we aim to find V such that is "aligned" to
% V_ref: is k = 1, then we set V_ref'*V > 0. If k = 2, we set V_ref_1'*V_1
% >0, V_ref_2'*V_1 = 0 and V_ref_2'*V_2 > 0. For general k, we need to set
% diagonal of the matrix V_ref'*V greater than zero, while all the entries
% below its diagonal equal to zero. Finally, these constraints have to be
% represented in terms of matrices and vectors such that
% Aine * z <= bine and Aequ * z = bequ, where z is the optimization vector


[n,k] = size(V_ref);
% Inequality
rows_ine = repmat(1:k,n,1); rows_ine = rows_ine(:);
cols_ine = transpose(repmat(1:n,k,1))+[0:k-1]*n; cols_ine = cols_ine(:);
vals_ine = -V_ref(:);
Aine = sparse(rows_ine,cols_ine,vals_ine,k,z_dim);
bine = sparse(k,1);
% Equality
rows_equ = zeros(n*k*(k-1)/2,1);
cols_equ = zeros(n*k*(k-1)/2,1);
vals_equ = zeros(n*k*(k-1)/2,1);
crows = 0; ctot = 0;
for ii = 1:k-1
    vals_equ_i = V_ref(:,ii+1:end);
    ind_i = ctot+1:ctot+length(vals_equ_i(:));
    rows_equ_i = transpose(crows+[1:k-ii]);
    rows_equ_i = transpose(repmat(rows_equ_i,1,n));
    cols_equ_i = 1+(ii-1)*n:ii*n;
    cols_equ_i = transpose(repmat(cols_equ_i,k-ii,1));
    rows_equ(ind_i) = rows_equ_i(:);
    cols_equ(ind_i) = cols_equ_i(:);
    vals_equ(ind_i) = vals_equ_i(:);
    ctot = ind_i(end); crows = crows + k - ii;
end
Aequ = sparse(rows_equ,cols_equ,vals_equ,k*(k-1)/2,z_dim);
bequ = sparse(k*(k-1)/2,1);
end
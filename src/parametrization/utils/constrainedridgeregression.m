function H = constrainedridgeregression(X,Y,XL,lI,V)
% Weighted Ridge Regression with Equality Constraints on the weights H
%
% Solves the least square fit for the model
%                   Y = X*H + e
% Each error is penalized according the vector L and the weights are
% penalized according with a coefficient lambda, defined above. 
% An eventual cross validation loop can find the best parameter set.
% Additionally, the weights satisfy the relation H*V = 0 (equality
% constraint). The code is intended for size(Y,2) to be large, hence
% sparse matrix are used.

% Normalization of features
[~,n] = size(X);
m = size(Y,2); l = n*size(V,2);

% Construction of weight constraints
A_w = []; I = [1:n*m]'; J = repmat([1:n]',m,1);
for ii = 1:size(V,2)
    V_tmp = repmat(transpose(V(:,ii)),n,1);
    A_w = horzcat(A_w,sparse(I,J,V_tmp(:)));
end

% Model Estimation
A11b = sparse(XL'*X+lI); ACell = repmat({A11b}, 1, m); 
A11 = blkdiag(ACell{:});
A = vertcat(horzcat(A11,A_w),horzcat(transpose(A_w),sparse(l,l)));
B = XL'*Y; b = sparse([B(:);   zeros(size(A_w,2),1)]);
sol = A\b;
% Final coefficient reshaping
H = full(transpose(reshape(sol(1:n*m,:),n,m)));
end

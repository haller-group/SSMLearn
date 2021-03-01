function [H,l_opt,Err] = ridgeregression(X,Y,L,idx_folds,l_vals)
% Implementation of the classic weighted ridge regression with Nfolds-
% crossvalidation for the regularizing parameters lambda_vec,c1_vec,c2_vec.

% Check on dimensions
[Ndata,Nfeatures] = size(X);
if Nfeatures > Ndata
    trans_option  = 1;
    X = transpose(X); Y = transpose(Y); L = transpose(L);
end
[~,Nfeatures] = size(X);

% Normalization & definition of regularizers
X_normal = max(abs(X),[],1); X = X./X_normal;
XL = L.*X;
if isempty(idx_folds) == 1
    l_opt = l_vals; Err = [];
else
    % Cross Validation
    Err = zeros(length(l_vals),length(idx_folds));
    for ii = 1:length(l_vals)
        lI = diag(l_vals(ii)*ones(1,Nfeatures));
        Err(ii,:) = ridgeregression_CVloop(X,Y,XL,lI,idx_folds);
    end
    [~,pos] = min(mean(Err,2)); l_opt = l_vals(pos);
end
H =  diag(X_normal.^(-1)) * ...
               (( XL'*X + diag(l_opt*ones(1, Nfeatures)) )\( XL'*Y ));
if trans_option == 1; H = transpose(H); end
end

%---------------------------Subfunctions---------------------------------

function Err = ridgeregression_CVloop(X,Y,XL,lI,idx_folds)
Err = zeros(1,length(idx_folds));
for ii = 1:length(idx_folds)
    ind_train = [idx_folds{:}];
    ind_train([idx_folds{ii}])=[];
    ind_test = [idx_folds{ii}];
    X_train = X(ind_train,:); Y_train = Y(ind_train,:);
    XL_train = XL(ind_train,:);
    X_test = X(ind_test,:); Y_test = Y(ind_test,:);
    sol = ( XL_train'*X_train + lI )\( XL_train'*Y_train );
    Err_i = Y_test-X_test*sol;
    Err(ii) = mean(sqrt(sum(Err_i.*conj(Err_i))));
end
end


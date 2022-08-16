function [W_r,phi,Dxphi,Dpphi,Expmat,l_opt,Err] = RidgeRegressionConstrainedParametric(t,X,P,Y,XPo,Yo,options)
% [W_r,phi,Dxphi,Expmat,l_opt,Err] = ...
%    RidgeRegressionConstrainedParametric(t,X,P,Y,Xo,Po,Yo,infoLin,options)
% Coinstrained Weighted Ridge Regression for parameter dependend maps. This
% functionfinds a map y = f(x,p) where f is expressed as a polynomial in x
% and p of maximal order options.Rs_PolyOrd in x and maximal order
% options.Rp_PolyOrd in p. It is assumed that 0 = f(0,0) and by setting
% option.origin_fixed to true we impose that 0 = f(0,p).
% We can optionally impose that the function takes known values
% yo = f(xo,po) and also constrain its linear part Ae = Dxf(xe,pe).

% Kill warning on singularity
warning('off','MATLAB:nearlySingularMatrix')
% Phase space dimension & Error Weghting
k = size(X,1); l = size(P,1); normF = sqrt(size(X,2));
L2 = (1+options.c1*exp(-options.c2*t)).^(-2);
options.L2 = L2; % = setfield(options,'L2',L2);
% Compute the exponents of the Taylor expansion and the polynomial map
ExpmatState = []; % State variable
for iOrd = 1 : options.Rs_PolyOrd
    ExpmatState = [ExpmatState; multivariateExponents(k,iOrd)];
end
if l>0
    ExpmatParam = zeros(1,l); % Parameters variable
    for iOrd = 1 : options.Rp_PolyOrd
        ExpmatParam = [ExpmatParam; multivariateExponents(l,iOrd)];
    end
    Expmat = []; % State + param. variables
    for iRow = 1:size(ExpmatState,1)
        Expmat = [Expmat; repmat(ExpmatState(iRow,:),size(ExpmatParam,1),1) ExpmatParam];
    end
    if options.origin_fixed == 0
        % Only skip the constant term so that R(0,0) = 0
        Expmat = [zeros(size(ExpmatParam,1)-1,k) ExpmatParam(2:end,1); Expmat];
    end
else
    Expmat = ExpmatState;
end
u = sym('u',[1 (k+l)]);
phi_sym = prod(u.^Expmat,2);
phi = matlabFunction(phi_sym,'Vars', {transpose(u)} );
Dxphi_sym = diff(phi_sym,['u' num2str(1)]);
if k > 1
    for iVar = 2:k
        Dxphi_sym = [Dxphi_sym diff(phi_sym,['u' num2str(iVar)])];
    end
end
Dxphi = matlabFunction(Dxphi_sym,'Vars', {transpose(u)} );
if l > 0
    Dpphi_sym = diff(phi_sym,['u' num2str(1+k)]);
    if l > 1
        for iVar = 2:l
            Dpphi_sym = [Dpphi_sym diff(phi_sym,['u' num2str(iVar+k)])];
        end
    end
    Dpphi = matlabFunction(Dpphi_sym,'Vars', {transpose(u)} );
else
    Dpphi = 0;
end
% Not constrained or constrained
linearParts = options.lin_part;
if isempty(XPo)+isempty(fieldnames(linearParts)) == 2
    [W_r,l_opt,Err] = ridgeRegression(phi([X; P])/normF,Y/normF,...
        options.L2,options.idx_folds,options.l_vals);
else
    % Fixed point constraints
    PhiC = phi(XPo);%reshape(repmat(phi(XPo),k,1),size(Expmat,1),size(Yo,2)*k);
    C = Yo;
    % Linear part constraints
    if isempty(fieldnames(linearParts)) == 0
        for iCon = 1:length(linearParts)
            if strcmp(options.type,'dynamics') == 1
                Wlin = linearParts(iCon).reducedDynamics;
            else
                Wlin = linearParts(iCon).parametrization;
            end
            if isfield(linearParts,'parameter') == 1
                p_i = linearParts(iCon).parameter;
            else
                p_i = [];
            end
            if isfield(linearParts,'reducedState') == 1
                x_i = linearParts(iCon).reducedState;
            else
                x_i = zeros(k,1);
            end
            C = [C Wlin]; PhiC = [PhiC Dxphi([x_i; p_i])];
        end
    end
    [W_r,l_opt,Err] = cRidgeRegression(phi([X; P])/normF,Y/normF,...
        options.L2,options.idx_folds,options.l_vals,PhiC,C);
end
% Re-activate warning on singularity
warning('on','MATLAB:nearlySingularMatrix')
end

%---------------------------Subfunctions---------------------------------

function [H,l_opt,Err] = cRidgeRegression(X,Y,L,idx_folds,l_vals,Xc,C)
% Implementation of the classic weighted ridge regression with Nfolds-
% crossvalidation for the regularizing parameters lambda_vec,c1_vec,c2_vec.

% Check on dimensions
[Ndata,Nfeatures] = size(X);
if Nfeatures > Ndata
    trans_option  = 1;
    X = transpose(X); Y = transpose(Y); L = transpose(L);
end
[~,Nfeatures] = size(X);
if size(Xc,2)~= Nfeatures
    Xc = transpose(Xc);
end
if size(C,1)~= size(Xc,1)
    C = transpose(C);
end
% Normalization & definition of regularizers
X_normal = max(abs(X),[],1); X = X./X_normal; XcO = Xc; Xc = Xc./X_normal;
XL = L.*X;
if isempty(idx_folds) == 1
    l_opt = l_vals; Err = [];
else
    % Cross Validation
    Err = zeros(length(l_vals),length(idx_folds));
    for ii = 1:length(l_vals)
        lI = diag(l_vals(ii)*ones(1,Nfeatures));
        Err(ii,:) = cridgeregression_CVloop(X,Y,XL,lI,idx_folds,Xc,C);
    end
    [~,pos] = min(mean(Err,2)); l_opt = l_vals(pos);
end
Amat = [(transpose(XL)*X + diag(l_opt*ones(1, Nfeatures))) ...
    transpose(Xc); Xc zeros(size(Xc,1))];
bvec = [transpose(XL)*Y; C];
sol = Amat\bvec; sol = sol(1:Nfeatures,:);
H =  diag(X_normal.^(-1)) * sol;
if trans_option == 1; H = transpose(H); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Err = cridgeregression_CVloop(X,Y,XL,lI,idx_folds,Xc,C)
Err = zeros(1,length(idx_folds)); ind_all = 1:size(X,1);
for ii = 1:length(idx_folds)
    ind_train = ind_all;
    ind_train([idx_folds{ii}])=[];
    ind_test = [idx_folds{ii}];
    X_train = X(ind_train,:); Y_train = Y(ind_train,:);
    XL_train = XL(ind_train,:);
    X_test = X(ind_test,:); Y_test = Y(ind_test,:);
    Amat = [(transpose(XL_train)*X_train + lI) transpose(Xc); ...
        Xc zeros(size(Xc,1))];
    bvec = [transpose(XL_train)*Y_train; C];
    sol = Amat\bvec; sol = sol(1:size(X_train,2),:);
    Err_i = Y_test-X_test*sol;
    Err(ii) = mean(sqrt(sum(Err_i.*conj(Err_i),2)));
end
end
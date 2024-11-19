function [V_e,IMParam,IM_param_info] = IMGeometryGraphT0(X,k,M,varargin)
% Construct the parametrization for an invariant manifold related to a
% fixed point (assumed to be the origin) as a graph the coordinates
% y = V_e'*x:
%
%               x = IM_para(y) = V_e*y + H phi(y)
%
% V_e is an orthonormal representation of the eigenspace to which the
% manifold is tangent at the origin. Moreover, phi(y) is a k-variate
% polynomial mapping from degree 2 to degree M and its coefficients H are
% orthogonal to V_e, i.e., V_e'*H = 0.
% A classic ridge penality can be inserted in the optimization of H.
% Specifically, reg_info{1} is the regularizer value, while reg_info{2}
% contains the normalization coefficients for the weights in H.
% Eventual error weighting can be performed in a slow manifold fashion: the
% cell err_info contains the time instances in err_info{1} and a vector
% c in err_info{2}. Each error is multiplied by a  coefficient
% (1+c(1)*exp(-c(2)*t)).^(-1) which is useful for imperfect measurements
% converging to the slow manifold which has to be identified.
% The best hyperparameter set may be obtained via crossvali-
% dation routines to be implemented at higher level.
%
% The code obtains the parametrization either if V_e is known or not. If 
% V_e is known, then the algorithm solves a ridge regression problem in
% closed form with constraints on the coefficients H. If V is unknown,
% then V_e and H are the result of a constrained optimization procedure.
% This latter is solved using the fmincon Matlab function (see the
% documentation for further details)
%
% REQUIRED INPUT
%    X    - matrix of dimension n x N, where n is the number of features
%           and N that of the number of data-points or a cell array of
%           dimension (N_traj,2) where the first column contains time
%           instances (1 x mi each) and the second column the trajectories
%           (n x mi each)
%    k    - invariant manifold dimension
%    M    - polynomial degree of phi for the parametrization
%
% OPTIONAL INPUT
%  varargin = options list: 'field1', value1, 'field2', value2, ... . The
%             options fields and their default values are:
%   'V_e' - real representation for the tangent space at the origin
%    'l'  - nonlinear coefficients regularization, default 0
%    'c1' - error coefficient for slow manifolds weighting
%           (1+c1*exp(-c2*t)).^(-1), default 0
%    'c2' - error coefficient for slow manifolds weighting
%           (1+c1*exp(-c2*t)).^(-1), default 0
%    't'  - time instances at which X data points are known, default 1.
%           Is set automatically if X is a cell array.
%    'V_e0' - optimizer initial condition for V, default []. Unless specified
%           , the algorithm uses the weighted k-dim. PCA of the data
%    'H0' - optimizer initial condition for H, default []. Unless
%           specified, the algorithm sets H0 from a constrained ridge 
%           regression using V as V_e0.
%    'Display' - default 'iter'
%    'OptimalityTolerance' - default 1e-5
%    'MaxIter' - default 100
%    'MaxFunctionEvaluations' - default 300
%    'SpecifyObjectiveGradient' - default true
%    'SpecifyConstraintGradient' - default true
%    'CheckGradients' - default false
%     the last seven options are for Matlab function fmincon.
%     For more information, check out its documentation.
%
% OUTPUT
%   V_e         - an orthonormal representation of the subspace to which
%                 the manifold is tangent at the origin
% IM_para       - function that parametrize the manifold as function of the
%                 the coordinates q
% OPTIONAL OUTPUT
% IM_param_info - contains information on the parametrization function
%
% Developed by Mattia Cenedese. Updated June 2021.

%--------------------------Main function-----------------------------------

% Default options and custom ones
optsParam = IMParametrizationOptions;
if rem(length(varargin),2) > 0 && length(varargin) > 1
    error('Error on input arguments. Missing or extra arguments.')
end
if nargin > 4
    for ii = 1:length(varargin)/2
        optsParam.(varargin{2*ii-1}) = varargin{2*ii};
    end
end
if iscell(X)==1
    Xcell = X; X = []; t = [];
    for ii = 1:size(Xcell,1)
        X = [X [Xcell{ii,2}]]; t = [t [Xcell{ii,1}]];
    end
    optsParam.t = t;
end
L = (1+optsParam.c1*exp(-optsParam.c2*optsParam.t)).^(-1);
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
V_e = optsParam.V_e;
if isempty(optsParam.V_e) == 0 % KNOWN V
    % Check if V is an orthogonal representation, if V'*V-I does not
    % satisfy a tolerance, we compute an orthogonal representation V
    if sum(sum(abs(transpose(V_e)*V_e-eye(k))))>k^2*1e-5
        V_e = orthogonalizeGramSchmidt(V_e);
    end
    if M < 2
        IMParam=@(q) V_e*q; H = []; Exp_mat = []; phi = @(q) [];
    else
        % Construct phi and the projection coordinates q
        [phi,Exp_mat] = multivariatePolynomial(k,2,M);
        Q = transpose(V_e) * X; Phi = phi(Q);
        idxs = sum(V_e==1,2);
        if sum(idxs) == k
            idxs = idxs==0;
            % Perform simple regression
            [H1,~,~] = ridgeregression(Phi,X(idxs,:),L,[],optsParam.l);
            H = zeros(size(X,1),size(Phi,1)); H(idxs,:) = H1;
        else
            % Perform constrained regression
            lI = (optsParam.l).*(max(abs(Phi),[],2).^(-2));
            H = constrainedRidgeRegression(transpose(Phi),...
                transpose(X-V_e*Q),transpose(Phi.*L),optsParam.l,V_e);
        end
        IMParam=@(q) V_e*q + H*phi(q);
    end
else % UNKNOWN V
    if M < 2 % Classic (weighted) PCA if polynomial degree = 1
        [~,~,V_e]=svds(transpose(X.*L),k);
        IMParam = @(q) V_e * q; H = []; Exp_mat = []; phi = @(q) [];
    else
        % Optimization over V and H.
        % The optimization vector is defined as z = [V(:); H(:)]
        [n,~] = size(X); XL = X.*L;
        
        % Construct phi and its derivatives. phi_dim is the number of
        % multivariate monomials in phi
        [phi,Exp_mat,D_phi_info] = multivariatePolynomial(k,2,M);
        phi_dim = size(Exp_mat,1);
        
        % Initial guess: classic (weighted) PCA with nonlinearities compu-
        % from this PCA estimate, unless specified by the user
        if isempty(optsParam.V_e0)==1
            [~,~,V_e0]=svds(transpose(XL),k);
        else
            V_e0 = optsParam.V_e0;
        end
        Q0 = transpose(V_e0) * X; Phi_0 = phi(Q0);
        if optsParam.l~= 0; optsParam.Phi_n = max(abs(Phi_0),[],2);
            lI = (optsParam.l).*(max(abs(Phi_0),[],2).^(-2));
        else
            lI = zeros(phi_dim,1);
        end
        if isempty(optsParam.H0)==1
            H0 = constrainedRidgeRegression(transpose(Phi_0),...
                transpose(X-V_e0*Q0),transpose(Phi_0.*L),lI,V_e0);
            H0 = H0(:);
        else
            H0 = optsParam.H0;
        end
        H0 = H0(:); z_0 = [reshape(V_e0,n*k,1); H0(1:n*phi_dim)];
        % Define function to minimize
        fun = @(z) fMinimize(z,X,n,k,phi,phi_dim,D_phi_info,L,...
            repmat(transpose(lI),n,1));
        
        % Define constraints
        % Linear constraints: ensure alignment with the initial guess for V
        [Aine,bine,Aequ,bequ] = alignmentLinearConstraint(length(z_0),V_e0);
        % Nonlinear equality constraints
        idx_info = extrasNonlinearConstraints(n,k,phi_dim);
        nonlincon = @(z) defineNonlinConstraints(z,n,k,phi_dim,idx_info);
        
        % Define options for the optimization algorithm
        opt_options = optimoptions('fmincon',...
        'Display',optsParam.Display,...
        'OptimalityTolerance',optsParam.OptimalityTolerance,...
        'MaxIter',optsParam.MaxIter,...
        'MaxFunctionEvaluations',optsParam.MaxFunctionEvaluations,...
        'SpecifyObjectiveGradient',optsParam.SpecifyObjectiveGradient,...
        'SpecifyConstraintGradient',optsParam.SpecifyConstraintGradient,...
         'CheckGradients',optsParam.CheckGradients);
        
        % Minimization process
        z_opt = fmincon(fun,z_0,Aine,bine,Aequ,bequ,[],[],nonlincon,...
            opt_options);
        % Final output
        V_e = reshape(z_opt(1:n*k),n,k);
        H = reshape(z_opt(1+n*k:end),n,phi_dim);
        IMParam = @(q) V_e * q + H * phi(q);
    end
end
IM_param_info = struct('map',IMParam,'polynomialOrder',M,...
    'dimension', k, 'tangentSpaceAtOrigin',V_e,...
    'nonlinearCoefficients',H,'phi',phi,...
    'exponents',Exp_mat,'l',optsParam.l,...
    'c1',optsParam.c1,'c2',optsParam.c2);
end

% Default options

function optsParam = IMParametrizationOptions
optsParam = struct('V_e',[],'l', 0,'c1',0,'c2',0,'t',1,...
    'V_e0',[],...
    'H0',[],...
    'Display','iter',...
    'OptimalityTolerance',1e-5,...
    'MaxIter',100,...
    'MaxFunctionEvaluations',300,...
    'SpecifyObjectiveGradient',true,...
    'SpecifyConstraintGradient',true,...
    'CheckGradients',false);
end

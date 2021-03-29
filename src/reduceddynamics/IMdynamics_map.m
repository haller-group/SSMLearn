function [R,iT,N,T,Maps_info] = IMdynamics_map(X_traj,varargin)
% Identification of the reduced dynamics in k coordinates, i.e. the map
%
%                        x_{k+1} = R(x_k)
%
% via a weighted ridge regression. R(x) = W_r * phi(x) where phi is a
% k-variate polynomial from order 1 to order M. Cross-validation can be
% performed on random folds or on the trajectories.
% Upon request, the dynamics is returned via a coordinate change, i.e.
%
%                         R = iT o N o T
%
% where iT, T and N depend on the selected style.
% If the style is selected as modal, then the coordinate change is a linear
% map that transforms the linear part of R into the diagonal matrix of its
% eigenvalues. 
% If the style is selected as normalform, then the functions seeks from
% data the maps iT, N and T such that the dynamics N is in normal form.
% This option is only available for purely oscillatory dynamics (i.e., the
% eigenvalues of the linear part of R are complex conjugated only). The
% normal form is detected by evaluating the small denominators that would
% occur when computing analytically the normal form from the knowledge of
% the vector field. These small denominators occur when the real part of
% the eigenvalues of the linear part of R is small or resonant. With small,
% it is intended to be smaller then a user-defined tolerance. To overcome
% tolerance-based issues, the user can enforce the detection of the
% coefficients of the normal form of the equivalent center manifold reduced
% order model, specifying eventual resonances among the frequencies.
%
% INPUT 
% X_traj - cell array of dimension (N_traj,2) where the first column 
%          contains time instances (1 x mi each) and the second column the 
%          trajectories (k x mi each). Sampling time is assumed to be
%          constant
%  varargin = polynomial order of R
%     or
%  varargin = options list: 'field1', value1, 'field2', value2, ... . The
%             options fields and their default values are:
%     'c1' - error coefficient for slow manifolds weighting
%           (1+c1*exp(-c2*t)).^(-1), default 0
%     'c2' - error coefficient for slow manifolds weighting
%           (1+c1*exp(-c2*t)).^(-1), default 0
% 'l_vals' - regularizer values for the ridge regression, default 0
% 'n_folds'- number of folds for the cross validation, default 0
% 'fold_style' - either 'default' or 'traj'. The former set random folds
%               while the latter exclude 1 trajectory at time for the cross
%               validation
% 'style' - default, modal or normalform
% 'nf_style' - 'center_mfld' or 'actual eigs'
% 'tol_nf' - parameter for the tolerance in the detection of small
%            denominators in the normal form. The tolreance is set to be
%            max(abs(Eig(R)))*(1-exp(-tol_nf*max(abs(real(log(eig(R))))))). 
%            Default value for tol_nf is 10
% 'frequencies_norm' - expected frequencies ratios (e.g. [1 2]) for 1:2
%                      resonance with the first and the second frequencies,
%                      by the default the code uses the actual values of
%                      of the frequencies. This option only works with the
%                      normal form style center manifold
% 'IC_nf' - initial condition for the optimization in the normal form.
%           0 (default): zero initial condition;
%           1: initial estimate based on the coefficients of R
%           2: normally distributed with the variance of case 1
% 'rescale' - rescale for the modal coordinates. 
%           0: no rescale
%           1 (default): the maximum amplitude is 0.5 (ratios kept)
%           2: the maximum amplitude of all coordinates is 0.5   
% 'Display' - default 'iter'
% 'OptimalityTolerance' - default 1e-4 times the number of datapoints
% 'MaxIter' - default 100
% 'MaxFunctionEvaluations' - default 300
% 'SpecifyObjectiveGradient' - default true
%     the last five options are for Matlab function fminunc.
%     For more information, check out its documentation.

if rem(length(varargin),2) > 0 && length(varargin) > 1
    error('Error on input arguments. Missing or extra arguments.')
end

% Reshape of trajectories into matrices
t   = []; % time values
X   = []; % coordinates at time k
X_1 = []; % coordinates at time k + 1
ind_traj = cell(1,size(X_traj,1)); idx_end = 0;
for ii = 1:size(X_traj,1)
    t_i = X_traj{ii,1}; X_i = X_traj{ii,2};
    t = [t t_i(1:end-1)]; X = [X X_i(:,1:end-1)]; X_1 = [X_1 X_i(:,2:end)];
    ind_traj{ii} = idx_end+[1:length(t_i)]; idx_end = length(t);
end
options = IMdynamics_options(nargin,varargin,ind_traj,size(X,2));
% Phase space dimension & Error Weghting
k = size(X_i,1); L2 = (1+options.c1*exp(-options.c2*t)).^(-2);
options = setfield(options,'L2',L2);

% Construct phi and its derivatives
disp('Estimation of the reduced dynamics... ')
[phi,Expmat] = multivariate_polynomial(k,1,options.R_PolyOrd);
[W_r,l_opt,Err] = ridgeregression(phi(X),X_1,options.L2,...
                                         options.idx_folds,options.l_vals);
R = @(x) W_r*phi(x);
R_info = assemble_struct(R,W_r,phi,Expmat,l_opt,Err);
options.l = l_opt;
fprintf('\b Done. \n')
% Find the change of coordinates desired
switch options.style
    case 'modal'
        % Linear transformation
        [V,~,~] = eig_sorted(W_r(:,1:k),1); iT = @(x) V\x; T = @(y) V*y;
        iT_info = assemble_struct(iT,inv(V),@(x) x,eye(k),[],[]);
        T_info = assemble_struct(T,V,@(y) y,eye(k),[],[]);
        % Nonlinear modal dynamics coefficients
        V_M = multivariate_polynomial_lintransf(V,k,options.R_PolyOrd);
        W_n = V\W_r*V_M; N = @(y) V\(W_r*phi(V*y));
        N_info = assemble_struct(N,W_n,phi,Expmat,[],[]);
    case 'normalform'
        disp('Estimation of the reduced dynamics in normal form...')
        Dt = t(2)-t(1);
        [V,D,d_cont] = eig_sorted(W_r(:,1:k),Dt);
        n_real_eig = sum(imag(d_cont)==0);
        if n_real_eig>0
            disp('Normal form not available. Returning modal style.')
            % Linear transformation
            iT = @(x) V\x; T = @(y) V*y;
            iT_info = assemble_struct(iT,inv(V),@(x) x,eye(k),[],[]);
            T_info = assemble_struct(T,V,@(y) y,eye(k),[],[]);
            % Nonlinear modal dynamics coefficients
            V_M = multivariate_polynomial_lintransf(V,k,options.R_PolyOrd);
            W_n = V\W_r*V_M; N = @(y) V\(W_r*phi(V*y));
            N_info = assemble_struct(N,W_n,phi,Expmat,[],[]);
        else
            if options.rescale == 1
               v_rescale = max(abs(V\X),[],2); 
               V = 2*V*diag(max(v_rescale(1:k/2))*ones(1,k));
            end
            if options.rescale == 2
               v_rescale = max(abs(V\X),[],2); 
               V = 2*V*diag(v_rescale);
            end
            % Initialize the normal form
            Maps_info_opt=initialize_nf_map(V,D,d_cont,W_r,X,X_1,Dt,...
                                                                  options);
            % Get normal form mappings T, N and T^{-1}
            Maps = dynamics_coordchange_nf(Maps_info_opt,options);
            % Final output
            T_info_opt = Maps.T; N_info_opt = Maps.N;
            iT_info_opt = Maps.iT;
            iT  = iT_info_opt.Map; N = N_info_opt.Map; T = T_info_opt.Map;
            iT_info = assemble_struct(iT,iT_info_opt.coeff,...
                              iT_info_opt.phi,iT_info_opt.Exponents,[],[]);
            N_info = assemble_struct(N,N_info_opt.coeff,N_info_opt.phi,...
                                               N_info_opt.Exponents,[],[]);
            T_info = assemble_struct(T,T_info_opt.coeff,T_info_opt.phi,...
                                               T_info_opt.Exponents,[],[]);
            iT_info.lintransf = inv(V); T_info.lintransf = V;
            % Display the obtained normal form
            fprintf('\n')
            disp(['The data-driven normal form dynamics reads:'])
            fprintf('\n')
            table_nf = disp_normalform(N_info_opt.coeff,...
                                                  N_info_opt.Exponents,Dt);
            disp(table_nf)
            if k == 2
            disp(['Notation: z is a complex number; z` is the ' ...
                  'complex conjugated of z; z^k is the k-th power of z.'])
            else
            disp(['Notation: each z_j is a complex number; z`_j is the '...
              'complex conjugated of z_j; z^k_j is the k-th power of z_j.'])    
            end
        end
    otherwise
        T = @(x) x; iT=@(y) y; N =@(y) R(y);
        T_info = assemble_struct(@(x) x,eye(k),@(x) x,eye(k),[],[]);
        iT_info = T_info; N_info = R_info;
end
Maps_info = struct('R',R_info,'iT',iT_info,'N',N_info,'T',T_info);
end

%---------------------------Subfunctions-----------------------------------

function str_out = assemble_struct(fun,W,phi,Emat,l,Err)
str_out = struct('Map',fun,'coeff',W,'phi',phi,'exponents',Emat,...
    'l_opt',l,'CV_error',Err);
end

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function Maps_info_opt=initialize_nf_map(V,D,d_cont,W_r,X,X_1,Dt,options)
% Preparation function for the estimate of the normal form maps. Based on
% the optimization properties, the functions seeks the coefficients of the
% normal form dynamics and sets to zero those coefficients for the
% transformation T^{-1}. Their indexes are stored in the output struct.
% The error at time instant k for the successive optimization is
%
% Err_k = Y_1-D*Y+W_it_nl*phi_it(Y_1)-W_n*phi_n(Y+W_it_nl*phi_it(Y))
% 
% and this function also precomputes the difference Y_1-D*Y and the
% transformations phi_it(Y_1) and phi_it(Y) for a more efficient
% optimization. The overall process consider complex numbers, and the
% conjugated are ignored.

% Modal transformation
k = size(W_r,1); ndof = k/2;
% Transformation of coordinates
Y = V\X; Y_1 = V\X_1;
Y_red = Y(1:ndof,:); Y_1_red = Y_1(1:ndof,:); 
d_red = diag(D); d_red=d_red(1:ndof);
V_M = multivariate_polynomial_lintransf(V,k,options.R_PolyOrd);
W_modal = V\W_r*V_M;
% Initialize nonlinear maps
[phi_it,Expmat_it] = multivariate_polynomial(k,2,options.iT_PolyOrd);
[phi_n,Expmat_n] = multivariate_polynomial(k,2,options.N_PolyOrd);
Phi_iT_Y = phi_it(Y); Phi_iT_Y_1 = phi_it(Y_1);
% Get the terms of the normal form
if strcmp(options.nf_style,'center_mfld')
    tol_nf = 1e-8;
    if isempty(options.frequencies_norm) == 1
        d_nf = exp(1i*imag(d_cont)*Dt);
    else
        d_nf = transpose(exp([+1i*options.frequencies_norm ...
                        -1i*options.frequencies_norm]*imag(d_cont(1))*Dt));
    end
else
    d_nf = diag(D); tol_nf = max(abs(d_nf))*...
                        (1-exp(-Dt*options.tol_nf*max(abs(real(d_cont)))));
end
lidx_n = find(abs(d_nf-transpose(phi_n(d_nf)))<tol_nf);
if options.R_PolyOrd>options.N_PolyOrd
    W_n_0 = W_modal(:,k+[1:size(Expmat_n,1)]);
else
    W_n_0 = [W_modal(:,k+1:end) zeros(k,size(Expmat_n,1)+k-size(W_r,2))]; 
end
if options.R_PolyOrd>options.iT_PolyOrd
    W_it_0 = -W_modal(:,k+[1:size(Expmat_it,1)]);
else
    W_it_0 = -[W_modal(:,k+1:end) zeros(k,size(Expmat_it,1)+k-size(W_r,2))]; 
end
W_it_0=W_it_0./(diag(D)-transpose(phi_it(diag(D))));
lidx_elim_it = lidx_n;
if options.iT_PolyOrd<options.N_PolyOrd
    lidx_elim_it(lidx_elim_it>numel(W_it_0)) = [];
end
lidx_it  = transpose(1:numel(W_it_0)); 
lidx_it(lidx_elim_it)  = [];
% Set the indexes for the coefficients of T^{-1} and N
[idx_it(:,1),idx_it(:,2)] = ind2sub(size(W_it_0),lidx_it); 
idx_it(idx_it(:,1)>ndof,:) = []; % Eliminate cc rows
W_it_0_up = W_it_0(1:ndof,:);
lidx_it_up = sub2ind(size(W_it_0_up),idx_it(:,1),idx_it(:,2)); 
% Eliminate useless exponents for N
[idx_n(:,1),  idx_n(:,2)] = ind2sub(size(W_n_0),lidx_n);  
idx_n(idx_n(:,1)>ndof,:) = []; % Eliminate cc rows
W_n_0_up = W_n_0(1:ndof,:);
lidx_n_up = sub2ind(size(W_n_0_up),idx_n(:,1),idx_n(:,2)); 
W_n_0_up(setdiff(1:numel(W_n_0_up),lidx_n_up)) = 0;
IDX_expnts = ones(size(W_n_0_up));
IDX_expnts(setdiff(1:numel(W_n_0_up),lidx_n_up)) = 0;
idx_expnts = find(sum(IDX_expnts,1));
[phi_n,Expmat_n,D_phi_n_info] = multivariate_polynomial_sel(k,2,...
                                             options.N_PolyOrd,idx_expnts);
W_n_0_up = W_n_0_up(:,idx_expnts); IDX_expnts = IDX_expnts(:,idx_expnts);
lidx_n_up = find(IDX_expnts);
[idx_n(:,1),  idx_n(:,2)] = ind2sub(size(W_n_0_up),lidx_n_up);  
% Initial condition for the optimization
if ndof == 1
    IC_opt_complex = transpose([W_it_0_up(lidx_it_up) W_n_0_up(lidx_n_up)]);   
else
    IC_opt_complex = [W_it_0_up(lidx_it_up); W_n_0_up(lidx_n_up)];
end
IC_opt = [real(IC_opt_complex); imag(IC_opt_complex)];
% Maps info
N_info_opt  = struct('phi',phi_n,'Exponents',Expmat_n,'idx',idx_n,... 'D_phi_info',D_phi_n_info
                                  'lidx',lidx_n_up);
iT_info_opt = struct('phi',phi_it,'Exponents',Expmat_it,'idx',idx_it,...
                                                        'lidx',lidx_it_up);
D_phi_n = cell2struct(D_phi_n_info,{'Derivative','Indexes'},2); 

% Final Output
Y_1_DY_red = Y_1_red - d_red.*Y_red;
Maps_info_opt = struct('IC_opt',IC_opt,'V',V,'d_r',d_red,'Yk_r',Y_red,...
                       'Yk_1_r',Y_1_red,'Yk_1_DYk_r',Y_1_DY_red,...
                       'Phi_iT_Yk',Phi_iT_Y,'Phi_iT_Yk_1',Phi_iT_Y_1,...
                       'iT',iT_info_opt,'N',N_info_opt,...
                       'D_phi_n_info',D_phi_n);
end

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Default options

function options = IMdynamics_options(nargin_o,varargin_o,idx_traj,Ndata)
options = struct('Style','default','R_PolyOrd', 1,'iT_PolyOrd',1,...
    'N_PolyOrd',1,'T_PolyOrd',1,'c1',0,...
    'c2',0,...
    'L2',[],...
    'n_folds',0,...
    'l_vals',0,...
    'idx_folds',[],'fold_style',[],...
    'style','default','nf_style','center_mfld','frequencies_norm',[],...
    'tol_nf',1e1,...
    'IC_nf',0,...
    'rescale',1,...
    'Display','iter',...
    'OptimalityTolerance',10^(-4-floor(log10(Ndata))),...
    'MaxIter',3000,...
    'MaxFunctionEvaluations',10000,...
    'SpecifyObjectiveGradient',true);
% Default case
if nargin_o == 2; options.R_PolyOrd = varargin_o{:};
    options.N_PolyOrd = varargin_o{:}; end
% Custom options
if nargin_o > 2
    for ii = 1:length(varargin_o)/2
        options = setfield(options,varargin_o{2*ii-1},...
            varargin_o{2*ii});
    end
    % Some default options for polynomial degree
    if strcmp(options.style,'normalform')==1 && ...
                options.iT_PolyOrd*options.N_PolyOrd*options.T_PolyOrd == 1
        options.N_PolyOrd = options.R_PolyOrd;
    end
    if strcmp(options.style,'normalform')==1 && options.N_PolyOrd > 1
       if options.T_PolyOrd*options.T_PolyOrd == 1
           options.T_PolyOrd = options.N_PolyOrd;
           options.iT_PolyOrd = options.N_PolyOrd;
       else
           PolyOrdM = max([options.T_PolyOrd options.iT_PolyOrd]);
           options.T_PolyOrd = PolyOrdM;
           options.iT_PolyOrd = PolyOrdM;
       end
    end
    % Fold indexes
    if options.n_folds > 1
        if strcmp(options.fold_style,'traj') == 1
            options = setfield(options,'n_folds',length(idx_traj));
            idx_folds = idx_traj;
        else
            idx_folds = cell(options.n_folds,1);
            ind_perm = randperm(Ndata);
            fold_size = floor(Ndata/options.n_folds);
            for ii = 1:options.n_folds-1
                idx_folds{ii} = ind_perm(1+(ii-1)*fold_size:ii*fold_size);
            end
            idx_folds{ii} = ind_perm(1+(ii-1)*fold_size:length(ind_perm));
        end
        options = setfield(options,'idx_folds',idx_folds);
    end
end
end


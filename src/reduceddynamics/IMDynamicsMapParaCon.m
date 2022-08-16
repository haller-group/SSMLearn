function RDInfo = IMDynamicsMapParaCon(xpData,varargin)
% RDInfo = IMDynamicsParaCon(etaData)
% Identification of the reduced dynamics in k coordinates and l parameters,
% i.e. the map
%
%                        x_{j+1} = R(x_j,p)
%
% via a weighted ridge regression. R(x,p) = W_r * phi(x,p) where phi is
% a (k+l)-variate polynomial featuring orders from 1 to order Mk in x and
% orders from 0 to order Ml in mu. Cross-validation can be
% performed on random folds or on the trajectories. It is assumed that
% 0 = R(0,0).
% The presence of non-trivial fixed points can be enforced. If we have that
% x_o = R(x_o,p_o), then we enforce the regression weights to be
% identified to satisfy x_o = W_r * phi(x_o,p_o) for each of the
% prescribed fixed points. With a specific flag, we can enforced the fact
% that the origin keeps being a fixed point for all parameter values, i.e.,
% 0 = R(0,p). The linear part of the dynamics at some locations can be also
% enforced, i.e. A_o  = DR(x_o,p_o), in an analogous manner.
%
% OUTPUTS
% RDInfo - struct containing the information of all these mapings
%
% INPUTS
% etaData - cell array of dimension (N_traj,3) where the first column
%          contains time instances (1 x mi each), the second column the
%          trajectories (k x mi each), the third column contain the
%          specificied parameters mk. Sampling time is assumed to be
%          constant
%  varargin = polynomial order of R in eta
%     or
%  varargin = options list: 'field1', value1, 'field2', value2, ... . The
%             options fields and their default values are:
%
%     'c1' - error coefficient for slow manifolds weighting
%           (1+c1*exp(-c2*t)).^(-1), default 0
%     'c2' - error coefficient for slow manifolds weighting
%           (1+c1*exp(-c2*t)).^(-1), default 0
% 'l_vals' - regularizer values for the ridge regression, default 0
% 'n_folds'- number of folds for the cross validation, default 0
% 'fold_style' - either 'default' or 'traj'. The former set random folds
%               while the latter exclude 1 trajectory at time for the cross
%               validation
% 'fixed_points' - default struct(), specify fixed points with an
%                  appropriate struct field:
%                  fixedPoints = struct('reducedState',[],'parameter',[]);
% 'origin_fixed' - default false, specify if 0 = R(0,p) for any parameter
%                  value
% 'lin_part' - default struct(), specify fixed points with an
%              appropriate struct field: linearParts = struct(
%              'reducedDynamics',[],'reducedState',[],'parameter',[]);

if rem(length(varargin),2) > 0 && length(varargin) > 1
    error('Error on input arguments. Missing or extra arguments.')
end

% Reshape of trajectories into matrices
t   = []; % time values
X   = []; % coordinates at time k
X_1 = []; % coordinates at time k + 1
P   = []; % parameters at time k
ind_traj = cell(1,size(xpData,1)); idx_end = 0;
for ii = 1:size(xpData,1)
    t_i = xpData{ii,1}; X_i = xpData{ii,2}; 
    t = [t t_i(1:end-1)]; X = [X X_i(:,1:end-1)]; X_1 = [X_1 X_i(:,2:end)];
    if size(xpData,2)>2 
        P_i = xpData{ii,3}; P = [P repmat(P_i,1,length(t_i)-1)]; end
    ind_traj{ii} = idx_end+[1:length(t_i)]; idx_end = length(t);
end
Dt = t(2)-t(1);
options = IMdynamics_options(nargin,varargin,ind_traj,size(X,2));
% Phase space dimension & Error Weghting
k = size(X,1); l = size(P,1);
% Add fixed points
fixedPoints = options.fixed_points;
if isempty(fieldnames(fixedPoints)) == 0
    nFixedPoints = length(fixedPoints);
    Xo = reshape([fixedPoints(:).reducedState],k,nFixedPoints);
    if isfield(fixedPoints,'parameter') == 1
        Po = reshape([fixedPoints(:).parameter],l,nFixedPoints);
    else
        Po = [];
    end
else
    Xo = []; Po = [];
end
% Regression
disp('Estimation of the reduced dynamics... ')
[W_r,phi,Dxphi,Dpphi,Expmat,l_opt,Err] = ...
    RidgeRegressionConstrainedParametric(t,X,P,X_1,[Xo; Po],Xo,options);
if l > 0
    R = @(x,p) W_r*phi([x; repmat(p,1,size(x,2))]);
    DxR = @(x,p) W_r*Dxphi([x; p]); DpR = @(x,p) W_r*Dpphi([x; p]); 
    R_info = assembleStruct(R,W_r,phi,Expmat,l_opt,Err);
    R_info.jacobianState = DxR; R_info.jacobianParameter = DpR;
else
    R = @(x) W_r*phi(x);
    DxR = @(x) W_r*Dxphi(x);
    R_info = assembleStruct(R,W_r,phi,Expmat,l_opt,Err);
    R_info.jacobianState = DxR;
end
options.l = l_opt;
fprintf('\b Done. \n')
T = @(x) x; iT=@(y) y; N = @(y,p) R(y,p);
T_info = assembleStruct(@(x) x,eye(k),@(x) x,eye(k));
iT_info = T_info; N_info = R_info; V = []; D = []; d_cont  = [];
RDInfo = struct('reducedDynamics',R_info,'inverseTransformation',...
    iT_info,'conjugateDynamics',N_info,'transformation',T_info,...
    'conjugacyStyle',options.style,'dynamicsType','map',...
    'mapTimeStep',Dt,'eigenvaluesLinPartMap',diag(D),...
    'eigenvaluesLinPartFlow',d_cont,'eigenvectorsLinPart',V);
end

%---------------------------Subfunctions-----------------------------------

function str_out = assembleStruct(fun,W,phi,Emat,varargin)
PolyOrder = sum(Emat(end,:));
if isempty(varargin) == 0
    str_out = struct('map',fun,'coefficients',W,'polynomialOrder',...
        PolyOrder,'phi',phi,'exponents',Emat,'l_opt',varargin{1},...
        'CV_error',varargin{2});
else
    str_out = struct('map',fun,'coefficients',W,'polynomialOrder',...
        PolyOrder,'phi',phi,'exponents',Emat);
end
end

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Default options

function options = IMdynamics_options(nargin_o,varargin_o,idx_traj,Ndata)
options = struct('style','default','Rs_PolyOrd', 1,'Rp_PolyOrd',0,...
    'c1',0,'c2',0,'fixed_points',struct(),'origin_fixed',false,'lin_part',struct(),...
    'L2',[],'n_folds',0,'l_vals',0,'idx_folds',[],'fold_style',[],'type','dynamics');
% Default case
if nargin_o == 2; options.R_PolyOrd = varargin_o{:}; end
% Custom options
if nargin_o > 2
    for ii = 1:length(varargin_o)/2
        options = setfield(options,varargin_o{2*ii-1},...
            varargin_o{2*ii});
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
            ii = ii+1;
            idx_folds{ii} = ind_perm(1+(ii-1)*fold_size:length(ind_perm));
        end
        options = setfield(options,'idx_folds',idx_folds);
    end
end
end


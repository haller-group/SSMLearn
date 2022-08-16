function IMInfo = IMGeometryParaCon(yData,etaData,varargin)
% [RDInfo,R,iT,N,T] = IMDynamicsParaCon(etaData)
% Identification of the parametrization of a manifold in k dynamical 
% reduced coordinates and l parameters, i.e. the map
%
%                          y = V(x,p)
%
% via a weighted ridge regression. V(x,p) = W_v * phi(x,p) where phi is
% a (k+l)-variate polynomial featuring orders from 1 to order Mk in x and
% orders from 0 to order Ml in mu. Cross-validation can be
% performed on random folds or on the trajectories. It is assumed that
% 0 = V(0,0).
% The presence of non-trivial fixed points can be enforced. If we have that
% y_o = V(x_o,p_o), then we enforce the regression weights to be
% identified to satisfy y_o = W_v * phi(x_o,p_o) for each of the
% prescribed fixed points. With a specific flag, we can enforced the fact
% that the origin keeps being a fixed point for all parameter values, i.e.,
% 0 = V(0,p). The linear part of the parametrization at some locations can 
% be also enforced, i.e. A_o  = DV(x_o,p_o), in an analogous manner.
%
% OUTPUTS
% IMInfo - struct containing the information of all these mapings
%
% INPUTS
% yData - cell array of dimension (N_traj,3) where the first column
%          contains time instances (1 x mi each), the second column the
%          trajectories (n x mi each) in the observable coordines, the third column contain the
%          specificied parameters mk. 
% etaData - cell array of dimension (N_traj,3) where the first column
%          contains time instances (1 x mi each), the second column the
%          trajectories (k x mi each), the third column contain the
%          specificied parameters mk. 
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
%                  fixedPoints = struct('observableState',[],...
%                  'reducedState',[],'parameter',[]);
% 'origin_fixed' - default false, specify if 0 = R(0,p) for any parameter
%                  value
% 'lin_part' - default struct(), specify fixed points with an
%              appropriate struct field: linearParts = struct(
%              'parametrization',[],'reducedState',[],'parameter',[]);

if rem(length(varargin),2) > 0 && length(varargin) > 1
    error('Error on input arguments. Missing or extra arguments.')
end

% Reshape of trajectories into matrices
t   = []; % time values
X   = []; % reduced coordinates at time k
Y = []; % observable coordinates at time k
P   = []; % parameters
ind_traj = cell(1,size(etaData,1)); idx_end = 0;
for ii = 1:size(etaData,1)
    n_i = length(etaData{ii,1});
    t = [t etaData{ii,1}]; X = [X etaData{ii,2}]; Y = [Y yData{ii,2}];
    if size(etaData,2)>2 
        P_i = etaData{ii,3}; P = [P repmat(P_i,1,n_i)]; end
    ind_traj{ii} = idx_end+[1:n_i]; idx_end = length(t);
end
options = IMGeometry_options(nargin,varargin,ind_traj,size(X,2));
% Phase space dimension & Error Weghting
k = size(X,1); l = size(P,1); n = size(Y,1);
% Add fixed points
fixedPoints = options.fixed_points;
if isempty(fieldnames(fixedPoints)) == 0
    nFixedPoints = length(fixedPoints);
    Xo = reshape([fixedPoints(:).reducedState],k,nFixedPoints);
    Yo = reshape([fixedPoints(:).observableState],n,nFixedPoints);
    if isfield(fixedPoints,'parameter') == 1
        Po = reshape([fixedPoints(:).parameter],l,nFixedPoints);
    else
        Po = [];
    end
else
    Xo = []; Po = []; Yo = [];
end
% Regression
disp('Estimation of the parametrization... ')
[W_v,phi,Dxphi,Dpphi,Expmat,l_opt,Err] = ...
    RidgeRegressionConstrainedParametric(t,X,P,Y,[Xo; Po],Yo,options);
V = @(x,p) W_v*phi([x; repmat(p,1,size(x,2))]); 
DxV = @(x,p) W_v*Dxphi([x; p]); 
fprintf('\b Done. \n')
paramInfo = struct('map',V,'coefficients',W_v,...
    'polynomialOrderState',options.Rs_PolyOrd,'polynomialOrderParam',...
     options.Rp_PolyOrd,'dimensionState', k,'dimensionParam', l, ...
    'phi',phi,'exponents',Expmat,'l_opt',l_opt,'CV_error',Err,...
    'derivativeState', DxV);
if l > 0
    DpV = @(x,p) W_v*Dpphi([x; p]); V_info.jacobianParameter = DpV;
end
IMInfo = struct('chart',struct(),'parametrization',paramInfo);
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

function options = IMGeometry_options(nargin_o,varargin_o,idx_traj,Ndata)
options = struct('style','default','Rs_PolyOrd', 1,'Rp_PolyOrd',0,...
    'c1',0,'c2',0,'fixed_points',struct(),'origin_fixed',false,'lin_part',struct(),...
    'L2',[],'n_folds',0,'l_vals',0,'idx_folds',[],'fold_style',[]);
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


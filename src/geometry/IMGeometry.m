function [IMInfo, IMChart, IMParam] = IMGeometry(yData, SSMDim, M, varargin)
% [IMInfo, IMChart, IMParam] = IMGeometry(yData, SSMDim, M)
% Description of the k-dim. invariant manifold geometry, in terms of
% coordinate chart and the parametrization. The default or natural method
% assumes that the manifold can be seen as a graph of a function for
% coordinates being the projection to the tangent space at the origin. This
% space and the resulting parametrization are computed with the function
% IMGeometryGraphT0. In this case, the chart is the projection to the
% tangent space at the origin and the parametrization is estimated as a
% polynomial model of order M. See its help for further descriptions.
% Alternatively, user-defined charts or coordinates can be used and the
% parametrization is estimated with a polynomial model of order M from the
% user-defined coordinates.
%
% REQUIRED INPUT
%    yData - matrix of dimension n x N, where n is the number of features
%           and N that of the number of data-points, or a cell array of
%           dimension (N_traj,2) where the first column contains time
%           instances (1 x mi each) and the second column the trajectories
%           (n x mi each)
%    SSMDim - invariant manifold dimension
%    M    - polynomial degree for the parametrization
%
% OPTIONAL INPUT
%  varargin = options list: 'field1', value1, 'field2', value2, ... . The
%             options fields and their default values are:
% 'style' - parametrization style. Default is 'natural'
% 'chart' - user-defined coordinate chart
% 'reducedCoordinates' - user-defined reduced coordinates. It must have the
%                        same dimensions as X.
%    'l'  - coefficients regularization, default 0
%    'c1' - error coefficient for slow manifolds weighting
%           (1+c1*exp(-c2*t)).^(-1), default 0
%    'c2' - error coefficient for slow manifolds weighting
%           (1+c1*exp(-c2*t)).^(-1), default 0
%    't'  - time instances at which X data points are known, default 1.
%           Is set automatically if X is a cell array.
% Other optional inputs are those of IMGeometryGraphT0.
%
% OUTPUT
% IMInfo - struct with information on the chart and on the parametrization
% IMChart - function that return reduced coordinate from a manifold point
% IMParam - function that parametrize the manifold as function of the
%           the reduced coordinates

% Default options and custom ones
optsGeomtery = struct('l', 0,'c1',0,'c2',0,'t',1,...
    'style','natural','chart',[],'reducedCoordinates',[],...
    'Ve',[],'outdof',[]);
if rem(length(varargin),2) > 0 && length(varargin) > 1
    error('Error on input arguments. Missing or extra arguments.')
end
if nargin > 4
    for ii = 1:length(varargin)/2
        optsGeomtery.(varargin{2*ii-1}) = varargin{2*ii};
    end
end
if isempty(optsGeomtery.chart) == 0 || ...
        isempty(optsGeomtery.reducedCoordinates) == 0
    optsGeomtery.style = 'custom';
end

% Compute manifold geometry
switch optsGeomtery.style
    % Describe the manifold as a graph over the tangent space at the origin
    case 'natural'
        [V,IMParam,paramInfo] = IMGeometryGraphT0(yData,SSMDim,M,varargin{:});
        IMChart = @(x) transpose(V)*x;
        chartInfo = struct('map',IMChart,'polynomialOrder',1);
        paramInfo.mapOut = @(q) q;

        % Describe the manifold with user-defined coordinates or their chart
    otherwise
        % Store chart
        Q = optsGeomtery.reducedCoordinates;
        IMChart = optsGeomtery.chart;
        chartInfo = struct('map',IMChart);
        
        % Compute reduced coordinates
        if iscell(yData)==1
            X_cell = yData; yData = []; t = [];
            for ii = 1:size(X_cell,1)
                yData = [yData [X_cell{ii,2}]]; t = [t [X_cell{ii,1}]];
            end
            optsGeomtery.t = t;
        end
        L = (1+optsGeomtery.c1*exp(-optsGeomtery.c2*optsGeomtery.t)).^(-1);
        if iscell(Q)==1
            Q_cell = Q; Q = [];
            for ii = 1:size(Q_cell,1)
                Q = [Q [Q_cell{ii,2}]];
            end
        end
        if isempty(Q) == 1 || size(Q,2)~=size(yData,2)
            Q = IMChart(yData);
        end
        
        % Compute parametrization
        if isempty(optsGeomtery.Ve) == 1
            [phi,Exp_mat] = multivariatePolynomial(SSMDim,1,M);
            Phi = phi(Q);
            [H,~,~] = ridgeRegression(Phi,yData,L,[],optsGeomtery.l);
        else
            linPart = optsGeomtery.Ve*Q;
            [phi,Exp_mat] = multivariatePolynomial(SSMDim,1,M); 
            Phi = phi(Q);
            [H,~,~] = ridgeRegression(Phi(size(Q,1)+1:end,:),...
                yData-linPart,L,[],optsGeomtery.l);
            H = [optsGeomtery.Ve H];
        end
        IMParam = @(q) H * phi(q);
        paramInfo = struct('map',IMParam,'polynomialOrder',M,...
            'dimension', SSMDim, 'tangentSpaceAtOrigin',H(:,1:SSMDim),...
            'nonlinearCoefficients',H(:,SSMDim+1:end),'phi',phi,...
            'exponents',Exp_mat,'l',optsGeomtery.l,...
            'c1',optsGeomtery.c1,'c2',optsGeomtery.c2);
        % Outdof parametrization to save memory
        if ~isempty(optsGeomtery.outdof)
            HOut = H(optsGeomtery.outdof,:);
            IMParamOut = @(q) HOut * phi(q);
            paramInfo.mapOut = IMParamOut;
            
        end

end
IMInfo = struct('chart',chartInfo,'parametrization',paramInfo);
end

function IC = getStaticResponseIC(K, M, F, loadvector, plotMesh, PlotFieldonDefMesh,varargin)
% getStaticResponseIC(K, M, F, loadvector, plotMesh, PlotFieldonDefMesh)
% getStaticResponseIC(K, M, F, loadvector, plotMesh, PlotFieldonDefMesh,w0)
% calculate static response with/without initial guess w0

if ~isempty(varargin)
    w0 = varargin{1};
else
    w0 = -K\loadvector(:,1); % linear initial guess
end
nTraj = size(loadvector, 2);
n = size(K, 1);
IC = zeros(2*n,nTraj);
options = optimoptions('fsolve', 'MaxFunctionEvaluations', 1000*n, ...
    'MaxIterations',1000*n, 'Display', 'off');
for iLoad = 1:nTraj
    f_eq = @(w) ( [zeros(n) M]*F(0,[w; zeros(n,1)]) + loadvector(:,iLoad) );
    [w0, ~, exitflag, output] = fsolve(f_eq, w0, options);    
    if exitflag <= 0
        disp(output);
        error('Warning: No solution found for loading configuration')
    end
    IC(:,iLoad) = [w0; zeros(n,1)];
end
if plotMesh
    figure; 
    PlotFieldonDefMesh(w0,200)
end
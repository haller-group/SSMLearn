function IC = getStaticResponse(K, M, F, loadvector, plotMesh, PlotFieldonDefMesh)

nTraj = size(loadvector, 2);
n = size(K, 1);
w0 = -K\loadvector(:,1); % linear initial guess
IC = zeros(2*n,nTraj);
options = optimoptions('fsolve', 'MaxFunctionEvaluations', 1000*n, 'Display', 'off');
for iLoad = 1:nTraj
    f_eq = @(w)([zeros(n) M]*F(0,[w; zeros(n,1)]) + loadvector(:,iLoad));
    [w0, ~, exitflag, output] = fsolve(f_eq, w0, options);
    if exitflag <= 0
        disp(output);
        error('Warning: No solution found for loading configuration')
    end
    IC(:,iLoad) = [w0; zeros(n,1)];
end
if plotMesh
    figure; PlotFieldonDefMesh(w0,200)
end
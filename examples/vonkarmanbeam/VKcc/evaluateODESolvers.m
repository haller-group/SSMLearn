%% Example setup

odetols = [1e-5 1e-6]
odesolvers = {@ode15s, @ode23tb}
nElements = 12
[M, C, K, fnl, f_vec, outdof, PlotFieldonDefMesh] = build_model(nElements);
n = size(M,1);    % mechanical dofs (axial def, transverse def, angle)
[F, lambda] = functionFromTensors(M, C, K, fnl);
%% Generation of Synthetic Data

loads = [2];
nTraj = size(loads, 2);
%% 

loadvector = loads.*f_vec;
IC = getStaticResponse(K, M, F, loadvector, 1, PlotFieldonDefMesh);

observable = @(x) x(outdof,:);
tEnd = 30;
nSamp = fix(50 * tEnd * abs(imag(lambda(1))) / (2*pi));
dt = tEnd/(nSamp-1);

ii = 1;
for odetol = odetols
    for iSol = 1:length(odesolvers)
        odesolver = odesolvers{iSol}
        ind = (ii-1)*nTraj+1:ii*nTraj;
        tols(ind) = odetol;
        tols(ind(1))
        tic
        xData(ind,:) = integrateTrajectories(F, observable, tEnd, nSamp, nTraj,...
            IC, 'odetol', odetol, 'odesolver', odesolver);
        toc
        for iTraj = ind
            fnames{iTraj} = functions(odesolver).function;
        end
        ii = ii + 1;
    end
end
%%
figure
for ii = 1:size(xData, 1)
    plot(xData{ii,1}, xData{ii,2}, 'DisplayName', ...
        [num2str(nElements), ' elements ', fnames{ii}, ' tol=' num2str(tols(ii))])
    hold on
end
legend
xlabel('time [s]')
ylabel('u [m]')
set(gca,'fontname','times')
set(gca,'fontsize',18)
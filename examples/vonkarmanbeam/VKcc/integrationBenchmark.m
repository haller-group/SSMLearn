clearvars
close all
clc
%% Point load

loads = 2;
nTraj = size(loads, 2);
elementList = 2:2:14;
odesolvers = {@ode15s, @ode23tb};
integrationTimes = zeros(size(elementList));

for iSol = 1:length(odesolvers)
fnames{iSol} = functions(odesolvers{iSol}).function;
for iEl = 1:length(elementList)
    nElements = elementList(iEl);
    [M, C, K, fnl, f_vec, outdof] = build_model(nElements);
    n = size(M,1);    % mechanical dofs (axial def, transverse def, angle)
    [F, lambda] = functionFromTensors(M, C, K, fnl);
    loadvector = loads.*f_vec;
    IC = getStaticResponse(K, M, F, loadvector, 0, 0);
    close all
    
    observable = @(x) x(outdof,:);
    tEnd = 30;
    nSamp = fix(50 * tEnd * abs(imag(lambda(1))) / (2*pi));
    dt = tEnd/(nSamp-1);
    tic
    integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC, 'odetol', 1e-3, 'odesolver', odesolvers{iSol});
    integrationTimes(iSol, iEl) = toc
end
end
%%
bar(elementList, integrationTimes,0.4)
set(gca,'YScale','log')
ylim([0.75*min(min(integrationTimes)), max(max(integrationTimes))*1.5])
xlabel('# elements')
ylabel('time [s]')
grid on
set(gca,'fontname','times')
set(gca,'fontsize',18)
legend(fnames)
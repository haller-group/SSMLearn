clearvars
close all
clc
%% Point load

loads = 1.75;
nTraj = size(loads, 2);
elementList = 2:2:14;
integrationTimes = zeros(size(elementList));

for iter = 1:length(elementList)
    nElements = elementList(iter);
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
    integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC);
    integrationTimes(iter) = toc
end

bar(elementList, integrationTimes,0.4,'y')
set(gca,'YScale','log')
ylim([0.75*min(integrationTimes), max(integrationTimes)*1.5])
xlabel('# elements')
ylabel('time [s]')
grid on
set(gca,'fontname','times')
set(gca,'fontsize',18)
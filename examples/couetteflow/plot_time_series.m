clearvars
close all
clc
load dataRe146
% the cell array aData contains a compressed version of the full flow
nTraj = 8;
nTrain = 1;
nTest = 7;
shape = [32, 35, 32];
%% Measure the flow in a single location: 
probe_index = [16, 14, 16];

customFigure('latex', true);
hold on;
for iTraj = 1 : nTrain
    full_flow = aData{iTraj,2}'*pcaComponents + pcaMean;
    % The U, V, W components of a flow field which have been stacked
    % on top of each other.
    U = full_flow(:, 1:prod(shape));
    U = reshape(U, [size(full_flow,1), shape(1), shape(2), shape(3)]);
    plot(U(:, probe_index(1), probe_index(2), probe_index(3)));
end
xlabel('$t$', 'Interpreter','latex');
ylabel('$U_{probe}$', 'Interpreter','latex');
title('Training trajectory');

customFigure('latex', true);
hold on;
for iTraj = 1 : nTest
    testIndex = indTest(iTraj);
    full_flow = aData{iTraj,2}'*pcaComponents + pcaMean;
    % The U, V, W components of a flow field which have been stacked
    % on top of each other.
    U = full_flow(:, 1:prod(shape));
    U = reshape(U, [size(full_flow,1), shape(1), shape(2), shape(3)]);
    plot(U(:, probe_index(1), probe_index(2), probe_index(3)), 'LineWidth', 0.9);
end
xlabel('$t$', 'Interpreter','latex');
ylabel('$U_{probe}$', 'Interpreter','latex');
title('Test trajectory')

%% The energy input rate and energy dissipation rate of the flow is also recorded
% this is contained in the cell array xData


% Visualize the training trajectories
customFigure('subPlot',[2 1],'latex', true);
subplot(211);
title('Training trajectories')

plot(xData{1,1}, xData{1,2}(1,:), '-', 'Linewidth', 0.9);
xlabel('$t$', 'Interpreter','latex');
ylabel('$I$', 'Interpreter','latex');
 
subplot(212);
plot(xData{1,2}(1,:), xData{1,2}(2,:));
xlabel('$I$', 'Interpreter','latex');
ylabel('$D$', 'Interpreter','latex');
% Visualize the testing trajectories   
customFigure('latex', true); colororder(parula(nTest))
title('Test trajectories') 
for iTraj = 1:nTest
    testIndex = indTest(iTraj);
    plot(xData{testIndex,1}, xData{testIndex,2}(1,:), '-', 'Linewidth', 1.2);
end
xlabel('$t$', 'Interpreter','latex');
ylabel('$I$', 'Interpreter','latex');
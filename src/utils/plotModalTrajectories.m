function plotModalTrajectories(yData, M, C, K, modes)
% plotModalTrajectories(yData, M, C, K, modes)
% Compute and plot trajectories in modal coordinates for a known mechanical 
% system. 
% 
% INPUT
% yData    {nTraj x 2}   Cell array of trajectories to plot. First column 
%                         contains time, second column contains state.
% M        n x n         Mass matrix.
% C        n x n         Damping matrix.
% K        n x n         Stiffness matrix.
% modes    2 x 1         Mode numbers to be plotted: 
%                          xaxis: real(modes(1)), 
%                          yaxis: imag(modes(1))
%                          zaxis: real(modes(2))

n = size(M,1);
A = [zeros(n), eye(n);
    -M\K,     -M\C];

customFigure();
[W, D] = eigSorted(full(A));
nTraj = size(yData, 1);
for iTraj = 1:nTraj
    etaData{iTraj,1} = yData{iTraj,1};
    etaData{iTraj,2} = W\yData{iTraj,2};
    plot3(real(etaData{iTraj,2}(modes(1),:)), imag(etaData{iTraj,2}(modes(1),:)), real(etaData{iTraj,2}(modes(2),:)))
end
xlabel(['Re(\eta_', num2str(modes(1)), ')'])
ylabel(['Im(\eta_', num2str(modes(1)), ')'])
zlabel(['Re(\eta_', num2str(modes(2)), ')'])
title('Trajectories in modal coordinates')
view(50,30)
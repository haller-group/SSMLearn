function plotModalTrajectories(yData, M, C, K, modes)

n = size(M,1);
A = [zeros(n), eye(n);
    -M\K,     -M\C];

figure
[W, D] = eig_sorted(full(A));
nTraj = size(yData, 1);
for iTraj = 1:nTraj
    etaData{iTraj,1} = yData{iTraj,1};
    etaData{iTraj,2} = W\yData{iTraj,2};
    plot3(real(etaData{iTraj,2}(modes(1),:)), imag(etaData{iTraj,2}(modes(1),:)), real(etaData{iTraj,2}(modes(2),:)))
    hold on
end
xlabel(['Re(\eta_', num2str(modes(1)), ')'])
ylabel(['Im(\eta_', num2str(modes(1)), ')'])
zlabel(['Re(\eta_', num2str(modes(2)), ')'])
title('Trajectories in modal coordinates')
set(gca,'fontname', 'times')
set(gca,'fontsize', 18)
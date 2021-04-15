function plotModalTrajectories(xData, A, modes)

figure
[W, D] = eig_sorted(full(A));
nTraj = size(xData, 1);
for iTraj = 1:nTraj
    etaData{iTraj,1} = xData{iTraj,1};
    etaData{iTraj,2} = W\xData{iTraj,2};
    plot3(real(etaData{iTraj,2}(modes(1),:)), imag(etaData{iTraj,2}(modes(1),:)), real(etaData{iTraj,2}(modes(2),:)))
    hold on
end
xlabel(['Re(\eta_', num2str(modes(1)), ')'])
ylabel(['Im(\eta_', num2str(modes(1)), ')'])
zlabel(['Re(\eta_', num2str(modes(2)), ')'])
title('Trajectories in modal coordinates')
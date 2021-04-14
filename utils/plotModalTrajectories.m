function plotModalTrajectories(xData, DS)

[W, D] = eig_sorted(full(DS.BinvA));
nTraj = size(xData, 1);
for iTraj = 1:nTraj
    etaData{iTraj,1} = xData{iTraj,1};
    etaData{iTraj,2} = W\xData{iTraj,2};
    plot3(real(etaData{iTraj,2}(1,:)), imag(etaData{iTraj,2}(1,:)), real(etaData{iTraj,2}(2,:)))
    hold on
end
xlabel('Re(\eta_1)')
ylabel('Im(\eta_1)')
zlabel('Re(\eta_2)')
title('Trajectories in modal coordinates')
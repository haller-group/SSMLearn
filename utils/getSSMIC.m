function [IC, mfd, DS, S] = getSSMIC(M, C, K, fnl, nTraj, ICRadius, SSMDim, seed)

[DS, S, mfd] = getSSM(M, C, K, fnl, SSMDim);

N = size(M,1);
IC = zeros(2*N, nTraj);
z = ICRadius * pickPointsOnHypersphere(nTraj, SSMDim, seed);
z = repelem(z(1:2:end,:) + 1j*z(2:2:end,:),2,1);
z(2:2:end) = conj(z(2:2:end));
for iTraj = 1:nTraj
    IC(:,iTraj) = real(getManifoldPoint(mfd, z(:,iTraj)));
end
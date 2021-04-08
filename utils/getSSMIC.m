function [IC, mfd, DS, S] = getSSMIC(M, C, K, fnl, nTraj, ICRadius, SSMDim, seed)

N = size(M,1);
DS = DynamicalSystem();
set(DS, 'M', M, 'C', C, 'K', K, 'fnl', fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
S.choose_E([1:SSMDim]);
mfd = S.compute_whisker(9);

IC = zeros(2*N, nTraj);
z = ICRadius * pickPointsOnHypersphere(nTraj, SSMDim, seed);
z = repelem(z(1:2:end,:) + 1j*z(2:2:end,:),2,1);
z(2:2:end) = conj(z(2:2:end));
for iTraj = 1:nTraj
    IC(:,iTraj) = real(getManifoldPoint(mfd, z(:,iTraj)));
end
function [IC, mfd, DS, S] = getSSMIC(M, C, K, fnl, nTraj, ICRadius)

N = size(M,1);
DS = DynamicalSystem();
set(DS, 'M', M, 'C', C, 'K', K, 'fnl', fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
S.choose_E([1,2]);
mfd = S.compute_whisker(3);

IC = zeros(2*N, nTraj);
for iTraj = 1:nTraj
    z = ICRadius*(-1)^(2*iTraj/nTraj);
    IC(:,iTraj) = real(getManifoldPoint(mfd, [z; conj(z)]));
end
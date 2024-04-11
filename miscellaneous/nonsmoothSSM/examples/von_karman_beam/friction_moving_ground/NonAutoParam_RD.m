function [v,r] = NonAutoParam_RD(n,V0,W0,A,forcingVectors_na)

% Cosine forcing

MAT_plus = @(fFreqs) (eye(2*n) - V0 * W0) * A - 1i * fFreqs*eye(2*n);
MAT_minus = @(fFreqs) (eye(2*n) - V0 * W0) * A + 1i * fFreqs*eye(2*n);


F_1_plus = @(fAmpls) 1/2 * fAmpls * forcingVectors_na;
F_1_minus = @(fAmpls) 1/2 * fAmpls * forcingVectors_na;

v_1_plus = @(fFreqs,fAmpls) MAT_plus(fFreqs)\(V0*W0 - eye(2*n)) * F_1_plus(fAmpls); 
v_1_minus =  @(fFreqs,fAmpls) MAT_minus(fFreqs)\(V0*W0 - eye(2*n)) * F_1_minus(fAmpls); 

v = @(t,fFreqs,fAmpls) real(v_1_plus(fFreqs,fAmpls) * exp(1i * fFreqs*t) + v_1_minus(fFreqs,fAmpls) * exp(-1i * fFreqs * t));
r = @(t,fFreqs,fAmpls) real(W0 * A * v(t,fFreqs,fAmpls) + W0 * (F_1_plus(fAmpls) * exp(1i*fFreqs * t) + F_1_minus(fAmpls) * exp(-1i*fFreqs * t)));

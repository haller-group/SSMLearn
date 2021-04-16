function [F, M, C, K, fnl, lambda] = oscillator(N, mass, stiffness, damping)
% Example dynamical system
% Oscillator chain
% Returns the anonynous evolution function and initial conditions

M = mass*eye(N);
K = stiffness*(2*eye(N) - diag(ones(1,N-1),-1) - diag(ones(1,N-1),1));
C = K; C(1,1) = 4; C(end,end) = 3; C = damping*C;

% nonlinearities
F2 = sptensor(zeros(N, 2*N, 2*N));
F2(1,N+1,N+1) = 0.33;    % q1dot^2
F3 = sptensor(zeros(N, 2*N, 2*N, 2*N));
F3(1,1,1,1) = 2;         % q1^3
F3(1,1,1,N+1) = 0.3;     % q1^2*q1dot^2 
F3(1,N+1,N+1,N+1) = 0.5; % q1dot^3
fnl = {F2, F3};

[F, lambda] = tensorfunction2(M, C, K, fnl);
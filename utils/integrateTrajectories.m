function xSim = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC, varargin)
% Computes trajectories of F starting at initial conditions IC.
% Anonymous function observable can be passed as @(x) x to save the full
% state space. To save only the first coordinate, pass @(x) x(1,:)
%
% OUTPUT
% xSim      (nTraj x 2) cell array of (1 x nSamp) time sequencies and 
%               (dim x nSamp) trajectories
opts = odeset('AbsTol',1e-7);
xSim = cell(nTraj, 2);
for iTraj = 1:nTraj
    fprintf('simulating trajectory %d of %d...\n', iTraj, nTraj)
    [t, x] = ode23tb(F, linspace(0, tEnd, nSamp), IC(:, iTraj), opts);
    xSim{iTraj,1} = t';
    xSim{iTraj,2} = observable(x');
end

function xSim = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC, varargin)
% Computes trajectories of F starting at initial conditions IC.
% Anonymous function observable can be passed as @(x) x to save the full
% state space. To save only the first coordinate, pass @(x) x(1,:)
%
% OUTPUT
% xSim      (nTraj x 2) cell array of (1 x nSamp) time sequencies and 
%               (dim x nSamp) trajectories

xSim = cell(nTraj, 2);
for iTraj = 1:nTraj
    fprintf('simulating trajectory %d of %d...\n', iTraj, nTraj)
    [t, x] = ode45(F, linspace(0, tEnd, nSamp), IC(:, iTraj));
    xSim{iTraj,1} = t';
    xSim{iTraj,2} = observable(x');
end

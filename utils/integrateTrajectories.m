function ySim = integrateTrajectories(F, observable, tEnd, nSamp, nTraj, IC, varargin)
% Computes trajectories of F starting at initial conditions IC.
% Anonymous function observable can be passed as @(x) x to save the full
% state space. To save only the first coordinate, pass @(x) x(1,:)
%
% OUTPUT
% ySim      (nTraj x 2) cell array of (1 x nSamp) time sequencies and 
%               (dim x nSamp) trajectories

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addParameter(p, 'odetol', 1e-4, validScalarPosNum);
addParameter(p, 'odesolver', @ode15s);
parse(p, varargin{:});

opts = odeset('RelTol', p.Results.odetol);
ySim = cell(nTraj, 2);
for iTraj = 1:nTraj
    fprintf('simulating trajectory %d of %d...\n', iTraj, nTraj)
    [t, x] = p.Results.odesolver(F, linspace(0, tEnd, nSamp), IC(:, iTraj), opts);
    ySim{iTraj,1} = t.';
    ySim{iTraj,2} = observable(x.');
end

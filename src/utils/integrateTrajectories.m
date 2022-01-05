function xData = integrateTrajectories(F, tspan, IC, nSamples, varargin)
% xData = integrateTrajectories(F, tspan, IC, nSamples)
% xData = integrateTrajectories(F, tspan, IC, nSamples, observable)
% Computes trajectories of F starting at initial conditions IC.
% The function handle observable can be passed as @(x) x to save the full
% state space. To save only the first coordinate, pass @(x) x(1,:)
%
% INPUT
% F          function handle	\dot(x) = F(t, x)
% tspan      (2 x 1) or scalar	timespan, or if scalar, endtime
% IC         (n x nTraj)    	one column per trajectory
% nSamples   int                number of evenly sampled datapoints
% observable function handle    optional, default identity. To observe the
%                               mth component, pass @(x) x(m,:);
%
% OUTPUT
% xData      {nTraj x 2}        cell array of (1 x nSamp) time sequencies 
%                               and (dim x nSamp) trajectories

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addOptional(p, 'observable', @(x) x(:,:));
addParameter(p, 'odetol', 1e-4, validScalarPosNum);
addParameter(p, 'odesolver', @ode15s);
parse(p, varargin{:});

if length(tspan) == 1; tspan = [0, tspan]; end
nTraj = size(IC, 2);
opts = odeset('RelTol', p.Results.odetol);
xData = cell(nTraj, 2);
for iTraj = 1:nTraj
    fprintf('simulating trajectory %d of %d...\n', iTraj, nTraj)
    [t, x] = p.Results.odesolver(F, linspace(tspan(1), tspan(2), nSamples), IC(:, iTraj), opts);
    xData{iTraj,1} = t.';
    xData{iTraj,2} = p.Results.observable(x.');
end

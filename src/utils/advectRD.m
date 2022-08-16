function [etaRec, zRec, zData] = advectRD(RDInfo, etaData)
% [yRec, etaRec, zRec, zData] = advect(IMInfo, RDInfo, etaData)
% Advect the initial conditions in etaData by the dynamics in RDInfo.
%   The reconstructed trajectories can be compared to the original ones
% with computeTrajectoryErrors to evaluate the model precision
%   
% INPUT
% RDInfo     struct         Reduced dynamics in map or flow form.
% etaData    {nTraj x 2}    Trajectories. Only the initial condition
%                           and times are used for reconstruction.
% OUTPUT
% etaRec     {nTraj x 2}    Advected trajectories on the manifold.
% zRec       {nTraj x 2}    Advected trajectories in the conjugate
%                           dynamics.
% zData      {nTraj x 2}    Original yData transformed to conjugate
%                           dynamics space.

invT = RDInfo.inverseTransformation.map;
T = RDInfo.transformation.map;
N = RDInfo.conjugateDynamics.map;
zData = transformTrajectories(invT, etaData);
if strcmp(RDInfo.dynamicsType, 'flow')
    zRec = integrateFlows(N, zData);
elseif strcmp(RDInfo.dynamicsType, 'map')
    zRec = iterateMaps(N, zData);
end
etaRec = transformTrajectories(T, zRec);
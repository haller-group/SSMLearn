function yRec = liftTrajectories(IMInfo, etaData)
% yRec = liftTrajectories(IMInfo, etaData)
% 
% Use the parametrization of a manifold to transform trajectories in
% reduced coordinates back to the observable space.
%   
% INPUT
% IMInfo     struct        from IMGeometry function. 
% etaData    {nTraj x 2}   Cell array of trajectories. First column 
%                          contains time, second column contains state.
% OUTPUT
% yRec       {nTraj x 2}   Trajectories in the observable space.
%
% Equivalent to:
%       transformTrajectories(IMInfo.parametrization.map, etaData);

yRec = transformTrajectories(IMInfo.parametrization.map, etaData);
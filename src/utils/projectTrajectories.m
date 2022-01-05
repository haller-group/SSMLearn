function etaData = projectTrajectories(IMInfo, yData)
% etaData = projectTrajectories(IMInfo, yData)
% 
% Use the chart of a manifold to transform trajectories in the observable
% space to reduced coordinates on the manifold.
%   
% INPUT
% IMInfo     struct        from IMGeometry function. 
% yData      {nTraj x 2}   Cell array of trajectories. First column 
%                          contains time, second column contains state.
% OUTPUT
% etaData    {nTraj x 2}   Trajectories in reduced coordinates.
%
% Equivalent to:
%       transformTrajectories(IMInfo.chart.map, yData);

etaData = transformTrajectories(IMInfo.chart.map, yData);
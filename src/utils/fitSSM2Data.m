function [IMInfo,etaData] = fitSSM2Data(yData,SSMOrder,Ve,We,indTrain,indTest,outdof)
%fitSSM2Data This function fits the trajectory data (yData) over an SSM 
% constructed around the given modal subspace (Ve, We) upto specfied order (SSMOrder)
% The fitting is performed over the training trajectory (indTrain) and the 
% fitting error is computed for the testing trajectory (indTest)

% Get projection or modal coordinates 
SSMDim = size(Ve,2);
etaData = yData;
nTraj = size(yData,1);
for iTraj = 1:nTraj
    etaData{iTraj,2} = We*yData{iTraj,2};    
end

plotReducedCoordinates(etaData);
legend({'Test set trajectory', 'Training set trajectory'})
if SSMDim>2
   view(3) 
end

% Compute nonlinear part of the parametrization
IMInfo = IMGeometry(yData(indTrain,:), SSMDim,SSMOrder,...
         'reducedCoordinates',etaData(indTrain,:),'Ve',Ve,'outdof',outdof); 
IMInfo.chart.map = @(x) We*x;                          

% Parametrization error on test trajectory
normedTrajDist = computeTrajectoryErrors(liftTrajectories(IMInfo,...
    etaData), yData);
staticNMTE = mean(normedTrajDist(indTest))*100; % in percentage

disp(['fitting error = ' num2str(staticNMTE) '%'])

if ~isempty(outdof) && SSMDim<=2
idxPlot = [outdof]; % 3D Plot: eta_1, eta_2 and idxPlot coordinate
plotSSMandTrajectories(IMInfo, idxPlot, yData(indTest,:), etaData(indTest,:))
view(-100,20); legend('off')
end
end


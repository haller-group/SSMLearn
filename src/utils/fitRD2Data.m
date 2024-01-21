function [RDInfo,yRec, etaRec, zRec] = fitRD2Data(etaData,Ae,ROMOrder,IMInfo,yData,indTest,indTrain,outdof,freqNorm)
% fitRD2Data This function fits a reduced order model with a given linear part (Ae) to the provided data
% (etaData) with upto user-specified polynomial order (ROMOrder). Then, ROM
% trajectories are simulated and lifted onto the given invariant manifold (IMInfo) error against
% the full system data is computed.


RDInfo = IMDynamicsMech(etaData(indTrain,:), ...
    'R_PolyOrd', 1,'N_PolyOrd', ROMOrder, 'style', 'normalform','R_coeff',Ae,'rescale',1,'frequencies_norm',freqNorm,'MaxIter',5e3);

% We transform the truncated initial condition of our test trajectory according to the obtained change of coordinates, and integrate our reduced order evolution rule to predict the development of the trajectory. 
[yRec, etaRec, zRec] = advect(IMInfo, RDInfo, yData);

% Evaluation of reduced dynamics
% The error NMTE is computed as the average distance of the predicted trajectory to the measured one in the full state space.
normedTrajDist = computeTrajectoryErrors(yRec, yData);
NMTE = mean(normedTrajDist(indTest))*100;
disp(['Normalized mean trajectory error = ' num2str(NMTE) '%'])

% We plot the true test set trajectory in the reduced coordinates and compare it to the prediction. 
plotReducedCoordinates(etaData(indTest,:), etaRec(indTest,:))
legend({'Test set (truncated)', 'Prediction'})
if size(Ae,1)==2
    % Plot SSM with trajectories in the normal form reduced coordinates
    plotSSMandTrajectories(IMInfo, outdof, yData(indTest,:), ...
        zRec(indTest,:), 'NFT', RDInfo.transformation.map)
    view(-100,20); legend('off')
else
    view(3)
end
end


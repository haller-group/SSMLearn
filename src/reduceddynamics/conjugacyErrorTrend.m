function errorInfo = conjugacyErrorTrend(etaData,indTrain,indTest,polyOrders)
% Evaluate the conjugacy error on training and test trajectory for a set
% of polynomial orders.

% The routines fits the normal form on the training data for each polynomial
% order and computes the squared conjugacy error, showing a picture with
% its trend varying the polynomial order
errorInfo = zeros(length(indTrain)+length(indTest),length(polyOrders));
for iOrder = 1:length(polyOrders)
    % Set order
ROMOrder = polyOrders(iOrder);
disp(['Evaluating ROM at order : ' num2str(ROMOrder) ' ...'])
% Compute normal form on training data
RDInfo = IMDynamicsFlow(etaData(indTrain,:), ...
    'R_PolyOrd', ROMOrder, 'style', 'normalform','fig_disp_nf',0,...
    'fig_disp_nfp',-1,'Display','off');
zData = transformTrajectories(RDInfo.inverseTransformation.map, etaData);
% Evaluate error
errorInfo(:,iOrder) = conjugacyError(zData,RDInfo);
%fprintf('\b Done. \n')
end
errorTrain = mean(errorInfo(indTrain,:),1);
errorTest = mean(errorInfo(indTest,:),1);
customFigure();
plot(polyOrders,errorTrain,'Linewidth',2,'Color',[0 0 0])
plot(polyOrders,errorTest,'Linewidth',2,'Color',[0.65 0.65 0.65])
xlabel('normal form polynomial order')
ylabel('conjugacy error')
legend('Training data','Test data','Location','NE')
set(gca,'yscale','log')
end

function error = conjugacyError(zData,RDInfo)
% Evaluate the conjugacy error for a given normal form model
error = zeros(size(zData,1),1);
% Reshape of trajectories into matrices and time derivative
N = RDInfo.conjugateDynamics.map;
for ii = 1:size(zData,1)
    t_in = zData{ii,1}; Z_in = zData{ii,2};
    [dZidt,Zi,~] = finiteTimeDifference(Z_in,t_in,3);
    deltai = dZidt-N(Zi);
    error(ii) = sum(sum((deltai.*conj(deltai))));
end

end
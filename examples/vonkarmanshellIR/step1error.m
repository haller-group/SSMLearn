% One step prediction error comparison

linSyst = @(y) Ae*y;
nMean = 1;
tMax = 10;
indTrajs = indTest;
errors = cell(length(indTrajs),3);
for ii = 1:length(indTrajs)
    iTraj = indTrajs(ii);
    idxMax = sum(xDataTrunc{iTraj,1}<tMax)+1;
    tCurr = xDataTrunc{iTraj,1}(1:idxMax);   
    xCurr = xDataTrunc{iTraj,2}(:,1:idxMax); 
    yCurr = IMInfo.chart.map(xCurr);
    normy = max(vecnorm(yCurr));
    stepping = 1;
    idxStep = 1:stepping:length(tCurr);
    error=zeros(2,length(idxStep)-1);
    for jj = 1:length(idxStep)-1
    
    xDataCurr = {tCurr(idxStep(jj):idxStep(jj+1)), ...
                 xCurr(:,idxStep(jj):idxStep(jj+1))};
    yDataCurr = {tCurr(idxStep(jj):idxStep(jj+1)), ...
                   yCurr(:,idxStep(jj):idxStep(jj+1))};

    [~, yCurrRecNL, ~] = advect(IMInfo, RDInfo, xDataCurr);

    yCurrRecL= integrateFlows(linSyst, yDataCurr);

    error(1,jj) = norm(yCurrRecNL{1,2}(:,end)-yDataCurr{1,2}(:,end));
    error(2,jj) = norm(yCurrRecL{1,2}(:,end)-yDataCurr{1,2}(:,end));
    end
    errors{ii,1} = tCurr(1:end-1);
    errors{ii,2} = error;
    errors{ii,3} = length(tCurr)-1;
end
%% 
customFigure;
colors = colororder; colSSMT = 5; colSSML = 7; colFOM = 1;
[~,refTraj] = min([errors{:,3}]);
tCurr = errors{refTraj,1};
errorNL = zeros(1,length(errors{refTraj,1}));
errorL = errorNL;
for ii = 1:length(indTrajs)
    error = errors{ii,2};
%     if nMean>1
%         errorL = errorL + movmean(error(2,:)/normy*100,nMean)/length(indTrajs);
%         errorNL = errorNL + movmean(error(1,:)/normy*100,nMean)/length(indTrajs);
%     else
        errorL = errorL + interp1(errors{ii,1},error(2,:),tCurr)/normy*100/length(indTrajs);
        errorNL = errorNL + interp1(errors{ii,1},error(1,:),tCurr)/normy*100/length(indTrajs);
%     end
end
plot(tCurr,errorL,'k','Linewidth',2)
plot(tCurr,errorNL,'color',colors(colSSML,:),'Linewidth',2)
plot(tCurr,movmean(errorL,100),'w--','Linewidth',2)
plot(tCurr,movmean(errorNL,100),'w--','Linewidth',2)

xlabel('$t \, [$s$]$','Interpreter','latex');
ylabel('PE$_1^{\mathrm{test}}$ [\%]','Interpreter','latex');
set(gca,'YScale','log')
xlim([0 8])
legend('Linear','SSMLearn')
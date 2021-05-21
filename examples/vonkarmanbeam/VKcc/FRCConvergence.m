clearvars

ROMOrdersFull = [7,9,11,13]
for iOrder = 1:length(ROMOrdersFull)
    FRC_full{iOrder} = load(['FRC_full_',num2str(ROMOrdersFull(iOrder)),'.mat']).FRC_full;
end
ROMOrdersData = [9]
for iOrder = 1:length(ROMOrdersData)
    FRC_data{iOrder} = load(['FRC_data_',num2str(ROMOrdersData(iOrder)),'.mat']).FRC_data;
    if iOrder == 1; load(['backbone_data',num2str(ROMOrdersData(iOrder)),'.mat'], 'frq', 'amp'); end
end
load FRC_NI.mat
pp = 146;
FRC_NI.omega = FRC_NI.omega(1:pp);
FRC_NI.amp = FRC_NI.amp(1:pp);
%%

figure; hold on; colors = colororder;
plot(frq, amp, 'k', 'DisplayName', ['Backbone curve O(', num2str(ROMOrdersData(1)), ')']);
for iOrder = 1:numel(FRC_full)
    plotFRC(FRC_full{iOrder}, colors(iOrder,:), ['SSMTool O(', num2str(ROMOrdersFull(iOrder)), ')'])
end
for iOrder = 1:numel(FRC_data)
    plotFRC(FRC_data{iOrder}, [0,0,0], ['SSMLearn O(', num2str(ROMOrdersData(iOrder)), ')'])
end

scatter1 = scatter(FRC_NI.omega,FRC_NI.amp,48,'b','MarkerFaceColor','c','DisplayName','Numerical integration');
scatter1.MarkerFaceAlpha = 0.6;
scatter1.MarkerEdgeAlpha = 1.0;
xlim([30,42])
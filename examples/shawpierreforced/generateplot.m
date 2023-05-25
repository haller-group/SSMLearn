

customFigure('latex', true);
set(groot,'defaultLegendInterpreter','latex');


surf(Q, P, Z3, 'FaceColor','#4B7935', 'edgecolor','none','FaceAlpha', 0.6, 'HandleVisibility', 'off')
for iTraj = indTest
    if iTraj == indTest(1)
    plot3(yData{iTraj,2}(1,:), yData{iTraj,2}(2,:), yData{iTraj,2}(3,:), '.', 'markersize', 6.5, 'HandleVisibility', 'off');

    else
    plot3(yData{iTraj,2}(1,:), yData{iTraj,2}(2,:), yData{iTraj,2}(3,:), '.', 'markersize', 7, 'HandleVisibility', 'off');
   
    end

end

plot3(yData{1,2}(1,end), yData{1,2}(2,end), yData{1,2}(3,end), 'r.', 'markersize', 30, 'Displayname', 'Upper orbit');
plot3(yData{3,2}(1,end), yData{3,2}(2,end), yData{3,2}(3,end), 'b.', 'markersize', 30, 'Displayname', 'Lower orbit');
plot3(0,0,0, 'k.', 'markersize', 30, 'Displayname', 'Saddle');

legend('location', 'east')

xlabel('$q_1$');
ylabel('$p_1$');
zlabel('$q_2$');
set(gcf, 'units', 'inches', 'position', [0.001, 0.001, 6.8, 4.5]);
view([-47 15]);

exportgraphics(gcf, 'ssmplot_shawpierre.pdf', 'Resolution', 1600);%,'ContentType','vector');



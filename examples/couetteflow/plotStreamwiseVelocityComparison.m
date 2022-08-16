function plotStreamwiseVelocityComparison(vxFull,pcaCompsShifted,ssmCompsShifted,pcaComponents,pcaMean)
% Compare streamwise velocity fields of the Couette flow

% recover the full flow from the PCA compressed data
full = (pcaCompsShifted)' * pcaComponents + pcaMean;
% the shape of the original data
shape = [32, 35, 32];
[U, ~, ~] = unravelField(full', shape);
% create the x and z axes
xx = linspace(0, 5*pi/2, 32);
zz = linspace(0, 4*pi/3, 32);
tiLay = tiledlayout(1,3, 'TileSpacing','loose');

% Full flow
ax1 = nexttile;
cMap = purpleorange(128);
colormap(cMap);
xlabel('x');
ylabel('z');
contourf(xx, zz, squeeze(vxFull(:,17, :)));

% Pca compressed velocity field
ax2 = nexttile;
colormap(cMap);
contourf(xx, zz, squeeze(U(:,17, :)));

% Velocity field recovered from the reduced order model
ax3 = nexttile;
colormap(cMap);
full = (ssmCompsShifted)' * pcaComponents + pcaMean; 
[U, ~, ~] = unravelField(full', shape);
contourf(xx, zz, squeeze(U(:,17, :)));colorbar;

%title(tiLay, 'Streamwise velocity');
title(ax1, 'Original velocity');
title(ax2, 'PCA compression');
title(ax3, 'SSM-based prediction')
for a=[ax1, ax2, ax3]
    xlabel(a, '$x$', 'Interpreter','latex');
    ylabel(a, '$z$', 'Interpreter','latex');
    set(a,'defaulttextinterpreter','latex')
    set(a, 'fontname', 'Latin Modern Roman')
    set(a, 'fontsize', 14)   
end

end
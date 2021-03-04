function plotSSM(radius, plotInds, SSMFunction)
% Draws the shape of an SSM in the space of components plotInds.

rho_plot = linspace(0,radius,50);
theta_plot = linspace(0,2*pi,51);
[rho_mesh,theta_mesh] = meshgrid(rho_plot,theta_plot);

U1 = rho_mesh.*cos(theta_mesh); U2 = rho_mesh.*sin(theta_mesh);
x_SSM = SSMFunction(transpose([U1(:) U2(:)]));
q1 = x_SSM(plotInds(1),:); q2 = x_SSM(plotInds(2),:); q3 = x_SSM(plotInds(3),:);

Q1 = reshape(q1,size(rho_mesh,1),size(rho_mesh,2));
Q2 = reshape(q2,size(rho_mesh,1),size(rho_mesh,2));
Q3 = reshape(q3,size(rho_mesh,1),size(rho_mesh,2));

h = surf(Q1,Q2,Q3,rho_mesh);
colormap spring
h.EdgeColor = 'none';
h.FaceAlpha = .7;

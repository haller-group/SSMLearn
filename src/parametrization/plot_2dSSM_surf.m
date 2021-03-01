function h = plot_2dSSM_surf(c_plot,Q,IM_para,out_domani_incr,...
                                                     grid_evals,color_surf)
% Plot a 2 dim. SSM as a surface in the three cordinates of c_plot. 
% If c_plot is a single value, then z = x_SSM(c_plot), while x and y are
% the parametrizing coordinates.

% Generate radial mesh in standard coordinates
rho_plot = linspace(0,1+out_domani_incr/100,grid_evals+1);
theta_plot = linspace(0,2*pi,grid_evals+1);
[rho_mesh,theta_mesh] = meshgrid(rho_plot,theta_plot);
U1_std = rho_mesh.*cos(theta_mesh); U2_std = rho_mesh.*sin(theta_mesh);
U_std = transpose([U1_std(:) U2_std(:)]);

% Construct transformation to actual coordinates and map the standard mesh
% to actual coordinates
[B,~] = red_coordinates_standardization(Q);
U = B(U_std);
x_SSM = IM_para(U); c_numb = length(c_plot);

% Different plotting options
switch c_numb
case 1
    c1 = U(1,:); c2 = U(2,:); c3 = x_SSM(c_plot(1),:);
case 2
    c1 = U_std(1,:); c2 = U_std(2,:); c3 = x_SSM(c_plot(1),:);
otherwise 
c1 = x_SSM(c_plot(1),:); c2 = x_SSM(c_plot(2),:); c3 = x_SSM(c_plot(3),:);
end
C1 = reshape(c1,size(rho_mesh,1),size(rho_mesh,2));
C2 = reshape(c2,size(rho_mesh,1),size(rho_mesh,2));
C3 = reshape(c3,size(rho_mesh,1),size(rho_mesh,2));
h = surf(C1,C2,C3,rho_mesh);
h.EdgeColor = 'none';
h.FaceAlpha = .5;
if color_surf == 0
    colormap cool
else
    h.FaceColor = color_surf;
end
view(3)
end
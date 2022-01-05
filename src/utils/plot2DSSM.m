function h = plot2DSSM(cPlot, Q, IMParam, outDomainIncr, ...
                                                   gridEvals, colorSurf)

% Plot a 2 dim. SSM as a surface in the three cordinates of cPlot. 
% If cPlot is a single value, then z = x_SSM(c_plot), while x and y are
% the parametrizing coordinates. Eventually, cPlot could be also set as a
% projection matrix to a 3D system, i.e., cPlot has size 3 x length(x_SSM).

% Generate radial mesh in standard coordinates
rho_plot = linspace(0, 1+outDomainIncr/100, gridEvals+1);
theta_plot = linspace(0, 2*pi, gridEvals+1);
[rho_mesh, theta_mesh] = meshgrid(rho_plot, theta_plot);
U1_std = rho_mesh.*cos(theta_mesh);
U2_std = rho_mesh.*sin(theta_mesh);
U_std = transpose([U1_std(:) U2_std(:)]);

% Construct transformation to actual coordinates and map the standard mesh
% to actual coordinates
[B,~] = rcoordinatesStandardization(Q);
U = B(U_std);
x_SSM = IMParam(U);

% Different plotting options
[rowsC,colsC] = size(cPlot);
if rowsC~=1 && colsC~=1
    c1 = cPlot(1,:)*x_SSM;
    c2 = cPlot(2,:)*x_SSM;
    c3 = cPlot(3,:)*x_SSM;
else
    c_numb = length(cPlot);
switch c_numb
    case 1
        c1 = U(1,:);
        c2 = U(2,:);
        c3 = x_SSM(cPlot(1),:);
    case 2
        c1 = U_std(1,:);
        c2 = U_std(2,:);
        c3 = x_SSM(cPlot(1),:);
    case 3
        c1 = x_SSM(cPlot(1),:);
        c2 = x_SSM(cPlot(2),:);
        c3 = x_SSM(cPlot(3),:);
    otherwise
        c1 = U(1,:);
        c2 = U(2,:);
        c3 = cPlot*x_SSM;    
end
end
C1 = reshape(c1,size(rho_mesh,1), size(rho_mesh,2));
C2 = reshape(c2,size(rho_mesh,1), size(rho_mesh,2));
C3 = reshape(c3,size(rho_mesh,1), size(rho_mesh,2));
h = surf(C1, C2, C3, rho_mesh);
h.EdgeColor = 'none';
h.FaceAlpha = .5;
if sum(colorSurf == [-1 -1 -1]) == 3
    colormap cool
else
    h.FaceColor = colorSurf;
end
view(3)
end
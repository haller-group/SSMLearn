function plotSSMandTrajectories(IMInfo, plotInds, yData, etaData, varargin)
% plotSSMandTrajectories(IMInfo, plotInds, yData, etaData, 'Margin', m)
%
% Draws the shape of a 1D or 2D SSM along with the trajectories in yData,
% in the space of components plotInds.
%
% INPUT
% IMInfo     struct         Manifold from IMGeometry function.
% plotInds   int or (3x1)   Either a vector of three components of the
%            or scalar      observable space to plot, or a single component
%            function or    to plot on the z axis. If a scalar value is
%            3 scalar       passed the x and y axes are aligned with the
%            functions      reduced coordinate directions. Analogous if
%            yPlot = g(y)   scalar functions g are given as input. The
%                           function should take yData{iTraj,2}, whose
%                           dimensions are n x m, as input and return as
%                           output a 1 x m or 3 x m matrices, where m are
%                           the time instances.
% yData      {nTraj x 2}    Cell array of trajectories. First column
%                           contains time, second column contains state.
% etaData    {nTraj x 2}    Cell array of trajectories. First column
%                           contains time, second column contains reduced
%                           states. Could also be zData.
% Margin     real, percent  Optional. Default 20. Controls the margin
%                           between the manifold edge and the maximum
%                           extent of trajectories
% ColorTraj  (nAmp x 3)     Optional RGB colors, one row per curve, or pass
%                           a single row vector for all curves.
% ColorSSM      (1 x 3)     Optional RGB color for the SSM, ...
%                           default [-1 -1 -1]='cool'
% NFT        function       Normal form transformation h, i.e., z = h(y)

% EXAMPLES
% Plot a 1D or 2D manifold in the first three components of the observable
% space:
%    plotSSMWithTrajectories(IMInfo, [1,2,3], yData, etaData)
% Plot a 2D manifold with y_1 on the z axis and the reduced coordinates
% eta_1 and eta_2 on the x and y axes, with a margin of 50 %:
%    plotSSMWithTrajectories(IMInfo, 1, yData, etaData, 'Margin', 50)

customFigure(); colors = colororder;
SSMDim = IMInfo.parametrization.dimension;
SSMFunction = IMInfo.parametrization.map;
if size(yData,1) ~= size(etaData,1)
    error('The trajectories in yData are different from those in etaData.')
end
p = inputParser;
addOptional(p, 'Margin', 20, @(x) isnumeric(x) && isscalar(x));
addOptional(p, 'ColorTraj', colors);
if SSMDim == 2
    addOptional(p, 'ColorSSM', [-1 -1 -1]);
else
    addOptional(p, 'ColorSSM', [0 0 1]);
end
addOptional(p, 'NFT', []);
parse(p, varargin{:});
colororder(p.Results.ColorTraj)

% Two dimensional SSM case
if SSMDim == 2
    % Plot manifold
    plot2dimSSM(plotInds, cat(2,etaData{:,2}), SSMFunction, p.Results.Margin, 50, p.Results.ColorSSM, p.Results.NFT);
    legend('SSM')
    % Plot trajectories
    for iTraj = 1:size(yData,1)
        if isnumeric(plotInds) == 1
            if length(plotInds) == 3
                plot3(yData{iTraj,2}(plotInds(1),:), yData{iTraj,2}(plotInds(2),:), yData{iTraj,2}(plotInds(3),:), 'LineWidth', 2, 'DisplayName', ['Traj. ' num2str(iTraj)])
                xlabel(['$y_{' num2str(plotInds(1)) '}$'], 'Interpreter', 'latex')
                ylabel(['$y_{' num2str(plotInds(2)) '}$'], 'Interpreter', 'latex')
                zlabel(['$y_{' num2str(plotInds(3)) '}$'], 'Interpreter', 'latex')
            elseif length(plotInds) == 1
                if isempty(p.Results.NFT) == 1
                    plot3(etaData{iTraj,2}(1,:), etaData{iTraj,2}(2,:), yData{iTraj,2}(plotInds,:), 'LineWidth', 2, 'DisplayName', ['Traj. ' num2str(iTraj)])
                    xlabel('$\eta_1$', 'Interpreter', 'latex')
                    ylabel('$\eta_2$', 'Interpreter', 'latex')
                    zlabel(['$y_{' num2str(plotInds) '}$'], 'Interpreter', 'latex')
                else
                    plot3(real(etaData{iTraj,2}(1,:)), imag(etaData{iTraj,2}(1,:)), yData{iTraj,2}(plotInds,:), 'LineWidth', 2, 'DisplayName', ['Traj. ' num2str(iTraj)])
                    xlabel('$\rho\cos\theta$', 'Interpreter', 'latex')
                    ylabel('$\rho\sin\theta$', 'Interpreter', 'latex')
                    zlabel(['$y_{' num2str(plotInds) '}$'], 'Interpreter', 'latex')
                end
            end
        else
            yPlot = plotInds(yData{iTraj,2});
            if size(yPlot,1) == 3
                plot3(yPlot(1,:), yPlot(2,:), yPlot(3,:), 'LineWidth', 2, 'DisplayName', ['Traj. ' num2str(iTraj)])
                xlabel('$g_1$', 'Interpreter', 'latex')
                ylabel('$g_2$', 'Interpreter', 'latex')
                zlabel('$g_3$', 'Interpreter', 'latex')
            elseif size(yPlot,1) == 1
                if isempty(p.Results.NFT) == 1
                    plot3(etaData{iTraj,2}(1,:), etaData{iTraj,2}(2,:), yPlot, 'LineWidth', 2, 'DisplayName', ['Traj. ' num2str(iTraj)])
                    xlabel('$\eta_1$', 'Interpreter', 'latex')
                    ylabel('$\eta_2$', 'Interpreter', 'latex')
                    zlabel('$g$', 'Interpreter', 'latex')
                else
                    plot3(real(etaData{iTraj,2}(1,:)), imag(etaData{iTraj,2}(1,:)), yPlot, 'LineWidth', 2, 'DisplayName', ['Traj. ' num2str(iTraj)])
                    xlabel('$\rho\cos\theta$', 'Interpreter', 'latex')
                    ylabel('$\rho\sin\theta$', 'Interpreter', 'latex')
                    zlabel('$g$', 'Interpreter', 'latex')
                end
            end
        end
    end
% One dimensional SSM case    
elseif SSMDim == 1
    % Plot manifold
    etaDataMat = cat(2,etaData{:,2});
    etaMin = min(etaDataMat) * (100+p.Results.Margin)/100;
    etaMax = max(etaDataMat) * (100+p.Results.Margin)/100;
    parametrization = linspace(etaMin, etaMax);
    if isempty(p.Results.NFT) == 1
        dVals = SSMFunction(parametrization);
    else
        dVals = SSMFunction(p.Results.NFT(parametrization));
    end
    mfdColor4 = 0.1; mfdColor = p.Results.ColorSSM; mfdWidth = 7;
    % Plot trajectories
    for iTraj = 1:size(yData,1)
        if isnumeric(plotInds) == 1
            if length(plotInds) == 3
                if iTraj == 1
                    mfdplot = plot3(dVals(plotInds(1),:), dVals(plotInds(2),:), dVals(plotInds(3),:), 'Color', mfdColor, 'LineWidth', mfdWidth, 'DisplayName', 'SSM');
                    mfdplot.Color(4) = mfdColor4; view(3)
                    xlabel(['$y_{' num2str(plotInds(1)) '}$'], 'Interpreter', 'latex')
                    ylabel(['$y_{' num2str(plotInds(2)) '}$'], 'Interpreter', 'latex')
                    zlabel(['$y_{' num2str(plotInds(3)) '}$'], 'Interpreter', 'latex')
                end
                plot3(yData{iTraj,2}(plotInds(1),:), yData{iTraj,2}(plotInds(2),:), yData{iTraj,2}(plotInds(3),:), 'LineWidth', 2, 'DisplayName', ['Traj. ' num2str(iTraj)])
            elseif length(plotInds) == 2
                if iTraj == 1
                    mfdplot = plot3(parametrization, dVals(plotInds(1),:), dVals(plotInds(2),:), 'Color', mfdColor, 'LineWidth', mfdWidth, 'DisplayName', 'SSM');
                    mfdplot.Color(4) = mfdColor4; view(3)
                    ylabel(['$y_{' num2str(plotInds(1)) '}$'], 'Interpreter', 'latex')
                    zlabel(['$y_{' num2str(plotInds(2)) '}$'], 'Interpreter', 'latex')
                    if isempty(p.Results.NFT) == 1
                        xlabel('$\eta_1$', 'Interpreter', 'latex')
                    else
                        xlabel('$z_1$', 'Interpreter', 'latex')
                    end
                end
                plot3(etaData{iTraj,2}(1,:), yData{iTraj,2}(plotInds(1),:), yData{iTraj,2}(plotInds(2),:), 'LineWidth', 2, 'DisplayName', ['Traj. ' num2str(iTraj)])
            else
                if iTraj == 1
                    mfdplot = plot(parametrization, dVals(plotInds(1),:), 'Color', mfdColor, 'LineWidth', mfdWidth, 'DisplayName', 'SSM');
                    mfdplot.Color(4) = mfdColor4;
                    ylabel(['$y_{' num2str(plotInds(1)) '}$'], 'Interpreter', 'latex')
                    if isempty(p.Results.NFT) == 1
                        xlabel('$\eta_1$', 'Interpreter', 'latex')
                    else
                        xlabel('$z_1$', 'Interpreter', 'latex')
                    end
                end
                plot(etaData{iTraj,2}(1,:), yData{iTraj,2}(plotInds(1),:), 'LineWidth', 2, 'DisplayName', ['Traj. ' num2str(iTraj)])
            end
        else
            yPlot = plotInds(yData{iTraj,2});
            if size(yPlot,1) == 3
                if iTraj == 1
                    gMfd = plotInds(dVals);
                    mfdplot = plot3(gMfd(1,:), gMfd(2,:), gMfd(3,:), 'Color', mfdColor, 'LineWidth', mfdWidth, 'DisplayName', 'SSM');
                    mfdplot.Color(4) = mfdColor4; view(3)
                    xlabel('$g_1$', 'Interpreter', 'latex')
                    ylabel('$g_2$', 'Interpreter', 'latex')
                    zlabel('$g_3$', 'Interpreter', 'latex')
                end
                plot3(yPlot(1,:), yPlot(2,:), yPlot(3,:), 'LineWidth', 2, 'DisplayName', ['Traj. ' num2str(iTraj)])
            elseif size(yPlot,1) == 2
                if iTraj == 1
                    gMfd = plotInds(dVals);
                    mfdplot = plot3(parametrization, gMfd(1,:), gMfd(2,:), 'Color', mfdColor, 'LineWidth', mfdWidth, 'DisplayName', 'SSM');
                    mfdplot.Color(4) = mfdColor4; view(3)
                    ylabel('$g_1$', 'Interpreter', 'latex')
                    zlabel('$g_2$', 'Interpreter', 'latex')
                    if isempty(p.Results.NFT) == 1
                        xlabel('$\eta_1$', 'Interpreter', 'latex')
                    else
                        xlabel('$z_1$', 'Interpreter', 'latex')
                    end
                end
                plot3(etaData{iTraj,2}(1,:), yPlot(1,:), yPlot(2,:), 'LineWidth', 2, 'DisplayName', ['Traj. ' num2str(iTraj)])
            else
                if iTraj == 1
                    gMfd = plotInds(dVals);
                    mfdplot = plot(parametrization, gMfd(1,:), 'Color', mfdColor, 'LineWidth', mfdWidth, 'DisplayName', 'SSM');
                    mfdplot.Color(4) = mfdColor4;
                    ylabel('$g$', 'Interpreter', 'latex')
                    if isempty(p.Results.NFT) == 1
                        xlabel('$\eta_1$', 'Interpreter', 'latex')
                    else
                        xlabel('$z_1$', 'Interpreter', 'latex')
                    end
                end
                plot(etaData{iTraj,2}(1,:), yPlot(1,:), 'LineWidth', 2, 'DisplayName', ['Traj. ' num2str(iTraj)])
            end
        end
    end
    
else
    disp('SSM plotting only available for 1D and 2D manifolds. Use PlotTrajectories or custom plotting codes to visualize trajectories.')
end

end

function h = plot2dimSSM(cPlot, Q, IMParam, outDomainIncr, ...
                                                 gridEvals, colorSurf, NFT)

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

if isempty(NFT) == 1
    % Construct transformation to actual coordinates and map the standard mesh
    % to actual coordinates
    [B,~] = rcoordinatesStandardization(Q);
    U = B(U_std);
    x_SSM = IMParam(U);
else
   U = U_std*max(abs(Q(1,:)));
   x_SSM = IMParam(NFT(transformationComplexConj(U(1,:)+1i*U(2,:))));
end

% Coordinates to plot
if isnumeric(cPlot) == 1
    if length(cPlot) == 3
        c1 = x_SSM(cPlot(1),:);
        c2 = x_SSM(cPlot(2),:);
        c3 = x_SSM(cPlot(3),:);
    elseif length(cPlot) == 1
        c1 = U(1,:);
        c2 = U(2,:);
        c3 = x_SSM(cPlot(1),:);
    end
else
    g_SSM = cPlot(x_SSM);
    if size(g_SSM,1) == 3
        c1 = g_SSM(1,:);
        c2 = g_SSM(2,:);
        c3 = g_SSM(3,:);
    elseif size(g_SSM,1) == 1
        c1 = U(1,:);
        c2 = U(2,:);
        c3 = g_SSM(1,:);
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

This is a preview of the livescript `vortexshedding.mlx`.

# Vortex Shedding behind a Cylinder

POD (Principal Orthogonal Decomposition) trajectories (projections of velocity and pressure data to a linear subspace) of a flow past a cylinder (2D simulations from CFD) at Reynolds number equal to 70. Initial conditions are close to the origin, which is unstable, from which these trajectories converge to a stable limit cycle (representing the vortex shedding regime.

![image_0.png
](README_images/image_0.png
)

```matlab:Code
clearvars
close all
clc

% Load and display data
load data
idx = [1 3 5]; nTraj = size(PODTraj,1);
customFigure('subPlot',[length(idx) 1]); colors = winter(nTraj);
for ii = 1:nTraj
    for jj = 1:length(idx)
        subplot(length(idx),1,jj); hold on; grid on; box on;
        plot(PODTraj{ii,1},PODTraj{ii,2}(idx(jj),:),'Linewidth',1,'Color',colors(ii,:))
        xlabel('$t \, [$s$]$','Interpreter','latex')
        ylabel(['$a_{' num2str(idx(jj)) '}$'],'Interpreter','latex')
        xlim([0 200])
    end
end
```

![figure_0.png
](README_images/figure_0.png
)

# Datadriven SSM

The first two POD coordinates are used to parametrize the unstable manifold connecting the unstable steady flow to the stable limit cycle. Hence we impose a user-defined chart in this example, and we look for a polynomial parametrization for the unstable manifold over these POD coordinates.

```matlab:Code
% Data embedding: there is no need of delay embedding in this case, but the
% function checks the data to be in suitable format
SSM_dim = 2; 
[yData,opts_embd] = coordinatesEmbedding(PODTraj,SSM_dim);
```

```text:Output
The embedding coordinates consist of the measured states.
```

```matlab:Code
% Split in training and testing dataset
indTrain = 1:8;
indTest = 9;

% Polynomial degree of the parametrization
ParamDeg = 18; 
V = [eye(SSM_dim); zeros(size(yData{1,2},1)-SSM_dim,SSM_dim)];
IMInfo = IMGeometry(yData(indTrain,:),SSM_dim,ParamDeg,'chart',@(x) V'*x);

% Compute relative root-mean-squared error on test set, normalized w.r.t.
% the maximum norm point and expressed in percentage
errorTestSet = computeParametrizationErrors(IMInfo, yData(indTest,:))*100
```

```text:Output
errorTestSet = 
      0.21774

```

# Plot and validation

Now that we have computed the eigenspace of the manifold, we pass to the reduced coordinates $y$ by projecting all trajectories onto the 2-dim. leading space POD. 

```matlab:Code
etaData = projectTrajectories(IMInfo, yData);
indPlot = 1:3:length(indTrain);
plotSSMWithTrajectories(IMInfo, [1 2 5], yData(indTrain(indPlot),:),'Margin',0,'Colors',winter(length(indPlot)))
view(-100,20); zlabel('$a_3$','Interpreter','latex')
```

![figure_1.png
](README_images/figure_1.png
)

# Reduced Dynamics

To accurately reconstruct the dynamics, we need a high-order normal form of the form $\dot{\rho} =c\left(\rho \right)\rho$, $\dot{\theta} =\omega \left(\rho \right)$, which is sought from data.

```matlab:Code
ROMOrder = 11;
RDInfo = IMDynamicsFlow(etaData(indTrain,:), ...
    'R_PolyOrd', ROMOrder, 'style', 'normalform','MaxIter',2e3);
```

```text:Output
Estimation of the reduced dynamics...  Done. 
Estimation of the reduced dynamics in normal form...
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
     0           1      0.000140317                        0.0033
     1           3      4.13791e-05             10        0.00171  
     2           5      1.68909e-05       0.344258        0.00155  
     3           6      8.32405e-06              1       0.000949  
     4           7      1.35293e-06              1       0.000248  
     5           8      9.22708e-07              1       0.000101  
     6           9      8.05183e-07              1       6.76e-05  
     7          11      7.73253e-07       0.321918        4.3e-05  
     8          12      7.36067e-07              1       1.83e-05  
     9          13      7.23117e-07              1       3.01e-06  
    10          15      7.20034e-07             10       1.24e-05  
    11          16      7.17192e-07              1       8.69e-06  
    12          17      7.12627e-07              1       7.84e-06  
    13          18      7.09503e-07              1       8.36e-06 
    ...
   786         791      2.66441e-08              1       1.27e-09  
   787         792      2.66441e-08              1       1.27e-09  
   788         793      2.66441e-08              1       1.27e-09  
   789         794      2.66441e-08              1       1.27e-09  
   790         795      2.66441e-08              1       1.27e-09  
   791         796      2.66441e-08              1       1.27e-09  
   792         797      2.66441e-08              1       1.27e-09  
   793         798       2.6644e-08              1       1.26e-09  
   794         799      2.66438e-08              1       1.26e-09  
   795         800      2.66432e-08              1       1.25e-09  
   796         801      2.66417e-08              1       1.23e-09  
   797         802      2.66379e-08              1       1.22e-09  
   798         803      2.66281e-08              1        1.8e-09  
   799         804       2.6603e-08              1       2.83e-09  
                       ...
```

We transform the truncated initial condition of our test trajectory according to the obtained change of coordinates, and integrate our reduced order evolution rule to predict the development of the trajectory. 

```matlab:Code
[yRec, etaRec, zRec] = advect(IMInfo, RDInfo, yData);
```

# Evaluation of reduced dynamics

The normalized mean trajectory error NMTE is computed as the average distance of the predicted trajectory to the measured one in the observable space, normalized by the observable vector with maximum norm in the data. 

```matlab:Code
normedTrajDist = computeTrajectoryErrors(yRec, yData);
NMTE = mean(normedTrajDist(indTest))*100 % in percentage
```

```text:Output
NMTE = 
       5.4115

```

We plot the true test set trajectory in the reduced coordinates and compare it to the prediction. 

```matlab:Code
plotReducedCoordinates(etaData(indTest,:), etaRec(indTest,:))
legend({'Test set (truncated)', 'Prediction'})
```

![figure_2.png
](README_images/figure_2.png
)

We also plot the measured and predicted tip displacement. The reduced model seems to do well on previously unseen data, provided that it is close to the 2D manifold.

```matlab:Code
coordPlot = 1;
plotTrajectories(yData(indTest,:), yRec(indTest,:), 'm', 'PlotCoordinate',coordPlot,'DisplayName', {'Test set', 'Prediction'})
xlabel('$t \, [$s$]$','Interpreter','latex')
ylabel(['$a_{ ' num2str(coordPlot) '}$'],'Interpreter','latex')
```

![figure_3.png
](README_images/figure_3.png
)

# Frequency and amplitude dynamics trends in time

```matlab:Code
time = zRec{indTest,1};
amplitudeNormalForm = abs(zRec{indTest,2}(1,:));
instDamping = RDInfo.conjugateDynamics.damping(amplitudeNormalForm);
instFrequency = RDInfo.conjugateDynamics.frequency(amplitudeNormalForm);

fig = customFigure();
yyaxis left
plot(time,instDamping,'Linewidth',2)
ylabel('$c \, [$1/s$]$','Interpreter','latex')
yyaxis right
plot(time,instFrequency,'Linewidth',2)
xlabel('$t \, [$s$]$','Interpreter','latex')
ylabel('$\omega \, [$rad/s$]$','Interpreter','latex')
xlim([0 200])
```

![figure_4.png
](README_images/figure_4.png
)

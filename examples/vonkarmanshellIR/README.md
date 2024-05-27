This is a preview of the livescript `vonkarmanshellIR.mlx`.

# Shallow-curved shell structure with geometric nonlinearities
  

See [1] for the details of this model, and [2] for the description of this example.

[1] Jain, S., \& Tiso, P. (2018). Simulation-free hyper-reduction for geometrically nonlinear structural dynamics: a quadratic manifold lifting approach. *Journal of Computational and Nonlinear Dynamics*, *13*(7), 071003. [https://doi.org/10.1115/1.4040021](https://doi.org/10.1115/1.4040021)

[2] Cenedese, M., Marconi, J., Haller, G., \& Jain, S. (2023). Data-assisted non-intrusive model reduction for forced nonlinear finite elements models. Preprint: [arXiv: 2311.17865](https://arxiv.org/abs/2311.17865) 

The finite element code taken from the following package:

Jain, S., Marconi, J., Tiso P. (2020). YetAnotherFEcode (Version v1.1). Zenodo. [http://doi.org/10.5281/zenodo.4011282](http://doi.org/10.5281/zenodo.4011282)

![image_0.png
](README_images/image_0.png
)

# Example setup

The $N$-degree of freedom dynamical system is of the form

$$
{M\ddot{q} }+{C\dot{q} }+{Kq}+f(q,{\dot{q} })=0
$$

where $f=\mathcal{O}(|q|^2 ,|{\dot{q} }|^2 ,|q||{\dot{q} }|)$ represents the nonlinearities and $M$, $C$, and $K$ are the $n\times n$ mass, stiffness, and damping matrices, respectively.

We rewrite the system in first-order form as

$$
{\dot{x} }=Ax+G(x)=F(x)
$$

with

  

$x=\left\lbrack \begin{array}{c}
q\\
\dot{q} 
\end{array}\right\rbrack ,~~A=\left\lbrack \begin{array}{cc}
0 & I\\
-M^{-1} K & -M^{-1} C
\end{array}\right\rbrack ,~~G(x)=\left\lbrack \begin{array}{c}
0\\
-M^{-1} f(x)
\end{array}\right\rbrack$.

```matlab:Code
clearvars 
close all
clc

% Setup colors
colors = colororder; colSSMT = 5; colSSML = 7; colFOM = 1;
```

# Generate model

```matlab:Code
nDiscretization = 10; % Discretization parameter (#DOFs is proportional to the square of this number)
epsilon = 0.1; % converge at order 5
[M,C,K,fnl,f_0,outdof,Model] = buildModel(nDiscretization);
```

```text:Output
Building FE model
```

![figure_0.png
](README_images/figure_0.png
)

```text:Output
Assembling M,C,K matrices
Applying boundary conditions
Solving undamped eigenvalue problem
```

![figure_1.png
](README_images/figure_1.png
)

![figure_2.png
](README_images/figure_2.png
)

```text:Output
Using Rayleigh damping
a = 2x1    
    0.0000
    0.3981

Assembling external force vector
Getting nonlinearity coefficients
Loaded tensors from storage
Total time spent on model assembly = 00:00:12
```

```matlab:Code
n = length(M); % number of degrees of freedom
disp(['Number of degrees of freedom = ' num2str(n)])
```

```text:Output
Number of degrees of freedom = 1320
```

```matlab:Code
disp(['Phase space dimensionality = ' num2str(2*n)])
```

```text:Output
Phase space dimensionality = 2640
```

# Linear modal analysis

```matlab:Code
if ~isfile('linpart.mat')
    m = 10;
    [W,A,V,lambda] = linearpart(M,C,K,m);
    save('linpart.mat',"W","A","V","lambda")
else
    load('linpart.mat')
end
```

# Define master modes and linear part of the dynamics

We initialize the base properties of the SSM, i.e., its linear part, which we know from the linear dynamics of the model. In this case, we target the slow four-dimensional SSM of the system, which features an internal resonance. 

```matlab:Code
masterModes = [1 3 2 4]; % Modal displacements and modal velocities
Ve = V(:,masterModes); % Mode shape
We = W(masterModes,:); % Projection to mode shape
Ae = full(We*A*Ve) % Reduced, linearized dynamics
```

```text:Output
Ae = 4x4    
1.0e+04 *

         0         0    0.0001    0.0000
         0         0    0.0000    0.0001
   -2.2268   -0.0000   -0.0001   -0.0000
   -0.0000   -8.9270    0.0000   -0.0001

```

```matlab:Code
SSMDim = length(masterModes);

displacementVector = zeros(1,n); displacementVector(outdof) = 1;
loadVector = f_0;  %  could also be set as modal ones
```

# Compare linear and nonlinear response via modal displacement

We characterize the linear and nonlinear regimes via a static modal analysis, which serves to pick appropriate initial conditions for the trajectories we need to learn the SSM from data.

Displacement along the first mode

```matlab:Code
iMode = 1; scalingfactor1 = 1e-1; nsteps = 50; outdof1 = outdof(1);
[phi1, relativeDiffForceNorm] = modal_analysis(Model,scalingfactor1,nsteps,outdof1,true,iMode);
```

```text:Output
Solving undamped eigenvalue problem
Mode shape
   -0.0941
    ...

Eigenfrequency
  149.2239
```

![figure_3.png
](README_images/figure_3.png
)

![figure_4.png
](README_images/figure_4.png
)

![figure_5.png
](README_images/figure_5.png
)

![figure_6.png
](README_images/figure_6.png
)

```text:Output
Displacement at output DOF: -0.0026277
```

Pick up two initial trajectories that has high expected nonlinear content

```matlab:Code
indIC1 = [nsteps, nsteps-1];
IC1 = [phi1*(scalingfactor1*indIC1/nsteps);zeros(n,length(indIC1))];
```

Displacement along the second mode

```matlab:Code
iMode = 2; scalingfactor2 = 1e-1; nsteps = 50; outdof2 = outdof(2);
[phi2, relativeDiffForceNorm2] = modal_analysis(Model,scalingfactor2,nsteps,outdof2,true,iMode);
```

```text:Output
Solving undamped eigenvalue problem
Mode shape
    0.0745
   ...

Eigenfrequency
  298.7812
```

![figure_7.png
](README_images/figure_7.png
)

![figure_8.png
](README_images/figure_8.png
)

![figure_9.png
](README_images/figure_9.png
)

![figure_10.png
](README_images/figure_10.png
)

```text:Output
Displacement at output DOF: 0.0031676
```

Pick up two initial trajectories that has high expected nonlinear content

```matlab:Code
indIC2 = [nsteps, nsteps-1];
IC2 = [phi2*(scalingfactor2*indIC2/nsteps);zeros(n,length(indIC2))];
```

Pick combination of initial conditions along the both the modes for training and testing

```matlab:Code
a = rand; b = rand; 
ICs = [IC1, IC2, a*IC1 + b*IC2];
indTrain = [1 3 5 6];
indTest = [2 4];
```

Define the linear regime at 1 % relative force

```matlab:Code
linearDisplacementReference = scalingfactor1*(sum(relativeDiffForceNorm<1)+1)/nsteps;
nonlinearDisplacementReference = scalingfactor1*max(indIC1)/nsteps;
desiredAmplitudeDecay = nonlinearDisplacementReference/linearDisplacementReference;

```

# **Generate decaying trajectories via time integration**

We define observables and timescales. The computation of integration time is estimated from the linear decay that gets from the defined nonlinear amplitude to linear regime. We set the sampling time to capture approximately a fixed number points per period on the faster time scale. Then, we integrate using the initial conditions we obtained from the static analysis. Here, we use a pre-computed data set to avoid excessive computations.

```matlab:Code
observable = @(x) x; % Observe the full phase space
slowTimeScale = 2*pi/abs(lambda(1));
fastTimeScale = 2*pi/abs(lambda(2));
% The computation of integration time is estimated from the linear decay that gets 
% from the nonlinear amplitude to linear regime.
newSimulation = false;
if newSimulation
    numberPeriodsSlow = floor(log(desiredAmplitudeDecay)/...
        (2*pi*(-real(lambda(1))/abs(lambda(1)))))
    endTime = numberPeriodsSlow*slowTimeScale;
    % Set the sampling time to capture approximately 50 points per period on the 
    % faster time scale.
    numberPeriodsFast = floor(endTime/fastTimeScale);
    numberPointsPerPeriod = 50;
    nSamp = numberPeriodsFast*numberPointsPerPeriod+1;
    dt = endTime/(nSamp-1);
    xData = integrateTrajectoriesGalphaDirect(Model, endTime, ICs, nSamp, observable);
    loadShape = f_0;
    DataInfo = struct('nElements', Model.Mesh.nElements, 'loadShape', loadShape);
    save('dataVKDecay2DGalphaModal.mat', 'DataInfo', 'xData', 'dt', 'endTime', 'nSamp','-v7.3')
else
    load("../../data/vonkarmanshellIR/dataVKDecay2DGalphaModal.mat")
end

```

# Visualize data

```matlab:Code
% This step is here not mandatory as we know that the SSM exists in the
% phase space (no embedding needed)
xData = coordinatesEmbedding(xData, SSMDim, 'ForceEmbedding', 1);
```

```text:Output
The embedding coordinates consist of the measured states.
```

Data filtering: We need to make sure that the data that we use to identify the manifold lies close to it. We can do this by plotting a spectrogram of the observables of interest. In general, there may be many vibratory modes present at first, but the faster ones quickly die out.

```matlab:Code
indPlot1 = indTrain(1);
indPlot2 = indTest(1);

% Visualize spectrograms for the first outdof
fig = customFigure();
fontName = 'helvetica';
fontSize = 16;
tiledlayout(3,1);
nexttile
showSpectrogram(xData(1,:), outdof(1), 1);
ylim([0,70])%abs(lambda(1))/2/pi*5])
set(gca, 'fontname', fontName)
set(gca, 'fontsize', fontSize)
xticks([])
nexttile
showSpectrogram(xData(3,:), outdof(1), 1);
ylim([0,70])%abs(lambda(1))/2/pi*5])
set(gca, 'fontname', fontName)
set(gca, 'fontsize', fontSize)
ylabel('frequency [Hz]')
xticks([])
nexttile
showSpectrogram(xData(5,:), outdof(1), 1);
ylim([0,70])%abs(lambda(1))/2/pi*5])
xlabel('time [s]')
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'power spectral density [1/Hz]';
set(gca, 'fontname', fontName)
set(gca, 'fontsize', fontSize)
```

![figure_11.png
](README_images/figure_11.png
)

We plot the observables of interest over time for closer inspection. 

```matlab:Code
fig = customFigure();
for iTraj = indTrain
    if iTraj == indTrain(1)
        plot(xData{iTraj,1}, xData{iTraj,2}(outdof(1),:),'DisplayName','Training data');
    else
        plot(xData{iTraj,1}, xData{iTraj,2}(outdof(1),:),'HandleVisibility','off');
    end
end
for iTraj = indTest
    if iTraj == indTest(1)
        plot(xData{iTraj,1}, xData{iTraj,2}(outdof(1),:),':','DisplayName','Test data');
    else
        plot(xData{iTraj,1}, xData{iTraj,2}(outdof(1),:),':','HandleVisibility','off');
    end
end
legend('location','NE')
xlabel('$t \, [$s$]$','Interpreter','latex'); ylabel('$q_{\mathrm{out}} \, [$m$]$','Interpreter','latex'); 
title('Generated data')
```

![figure_12.png
](README_images/figure_12.png
)

# Truncate transient data from trajectories

We must however remove the first transient to fulfill the assumption that trajectories lie close to the SSM. We keep only the time interval |sliceInt|.

```matlab:Code
sliceInt = [5*slowTimeScale, endTime];
xDataTrunc = sliceTrajectories(xData, sliceInt);
```

# Datadriven manifold fitting

The measured trajectories are initialized to lie close to the manifold of interest that is tangent at the origin to the eigenspace spanned by the columns of $V_e$. 

As we also know the projection $W_e$ to this eigenspace, we define the modal coordinates as $y=W_e x$. These are the reduced coordinates for our graph style parametrization of the manifold, gauranteed to exists near the origin. We then use the data to learn the nonlinear feature of the manifold geometry, represented via polynomials. Indeed, we seek the $2N\times m_M$ polynomial coefficient matrix $H$ of the manifold on the form

> $x=V_e y+H{{\phi }}_{m,2:M} (y)$,

where the function ${{\phi }}_{m,2:M} (y)$ computes a vector of all $m_M$ monomials from orders 2 up to $M$ of an $m$-vector $y$. From SSM theory, the tangent space of the manifold is $V_e$. The coefficients $\mathbf{H}$are obtained via least squares regression.

```matlab:Code
SSMOrder = 2;

% Get projection or modal coordinates 
SSMDim = size(Ve,2);
yDataTrunc = xDataTrunc;
nTraj = size(xDataTrunc,1);
for iTraj = 1:nTraj
    yDataTrunc{iTraj,2} = We*xDataTrunc{iTraj,2};    
end

% Plot reduced coordinates
plotReducedCoordinates(yDataTrunc);
if SSMDim>2
   view(3) 
end
```

![figure_13.png
](README_images/figure_13.png
)

```matlab:Code

% Compute nonlinear part of the parametrization
IMInfo = IMGeometry(xDataTrunc(indTrain,:), SSMDim,SSMOrder,...
         'reducedCoordinates',yDataTrunc(indTrain,:),'Ve',Ve,'outdof',outdof); 
IMInfo.chart.map = @(x) We*x;                          

% Parametrization error on test trajectory
normedTrajDist = computeTrajectoryErrors(liftTrajectories(IMInfo,...
    yDataTrunc), xDataTrunc);
staticNMTE = mean(normedTrajDist(indTest))*100; % in percentage

disp(['Reconstruction error = ' num2str(staticNMTE) '%'])
```

```text:Output
Reconstruction error = 0.94102%
```

```matlab:Code

% Visualize the trajectories in different coordinates
fig = customFigure();
for iTraj = indTrain
    plot3(xData{iTraj,2}(outdof(1),:),xData{iTraj,2}(outdof(1)+size(M,1),:),xData{iTraj,2}(outdof(2),:))
end
xlabel('$q_{A}$ [m]','interpreter','latex')
ylabel('$\dot{q}_{A}$ [m/s]','interpreter','latex')
zlabel('$q_{B}$ [m]','interpreter','latex')
legend('Trajectory 1','Trajectory 2','Trajectory 3')
view(3)
```

![figure_14.png
](README_images/figure_14.png
)

```matlab:Code

fig = customFigure('subPlot',[3 2]);
for iTraj = [1 3 5]
    subplot(3,2,iTraj)
    plot(xData{iTraj,1}, xData{iTraj,2}(outdof(1),:));
    xlabel('$t \, [$s$]$','Interpreter','latex'); 
    ylabel('$q_{out,1} \, [$m$]$','Interpreter','latex'); 
    subplot(3,2,iTraj+1)
    plot(xData{iTraj,1}, xData{iTraj,2}(outdof(2),:));
    xlabel('$t \, [$s$]$','Interpreter','latex'); 
    ylabel('$q_{out,2} \, [$m$]$','Interpreter','latex'); 
end
```

![figure_15.png
](README_images/figure_15.png
)

# Reduced dynamics on the manifold

We compute a model for the reduced dynamics with the truncated training data projected onto the manifold. The function `IMDynamicsMech` finds the dynamics (considering that we know the linear part)

> $\dot{y} =A_e y+H_r {{\phi }}_{m,2:M} (y)$,

where ${\phi }(y)$ again computes a vector of all monomials of $u$, and $H_r$ is a matrix of polynomial coefficients of the form

$$
H_r =\left\lbrack \begin{array}{c}
0\\
W_r 
\end{array}\right\rbrack
$$

Then, we find from data the extended normal form, which expedites system characterization and computation of periodic forced responses.

```matlab:Code
ROMOrder = 3;
freqNorm = [1 2];

RDInfo = IMDynamicsMech(yDataTrunc(indTrain,:), ...
    'R_PolyOrd', 1,'N_PolyOrd', ROMOrder, 'style', 'normalform', ...
    'R_coeff',Ae,'rescale',1,'frequencies_norm',freqNorm,'MaxIter',5e3);
```

```text:Output
Estimation of the reduced dynamics...  Done. 
Estimation of the reduced dynamics in normal form...
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
     0           1          6.16336                          29.1
     1           3           4.2627     0.00343599           27.1  
     2           4          3.00363              1           13.9  
     3           5         0.876018              1           7.36  
     4           6          0.67696              1           3.61  
     5           7         0.630143              1           1.83  
     6           8         0.587024              1            2.2  
     7           9         0.504308              1           4.57  
     8          10         0.434066              1           4.08  
     9          11         0.389986              1           1.81  
    10          12         0.365842              1           1.43  
    11          13         0.343188              1           2.02  
    ... 
   789         791        0.0252652              1       0.000291  
   790         792        0.0252651              1       0.000253  
   791         793         0.025265              1       0.000216  
   792         794        0.0252648              1       0.000192  
   793         795        0.0252643              1       0.000204  
   794         796        0.0252637              1       0.000255  
   795         797        0.0252633              1       0.000319  
   796         798        0.0252631              1       0.000298  
   797         799        0.0252631              1       0.000275  
   798         800        0.0252631              1        0.00027  
   799         801        0.0252631              1       0.000262  
                       ...
```

![figure_16.png
](README_images/figure_16.png
)

```matlab:Code

% We transform the truncated initial condition of our test trajectory according to 
% the obtained change of coordinates, and integrate our reduced order evolution rule 
% to predict the development of the trajectory. 
[yRec, etaRec, zRec] = advect(IMInfo, RDInfo, xDataTrunc);

% Evaluation of reduced dynamics
% The error NMTE is computed as the average distance of the predicted trajectory 
% to the measured one in the full state space.
normedTrajDist = computeTrajectoryErrors(yRec, xDataTrunc);
NMTE = mean(normedTrajDist(indTest))*100;
disp(['Normalized mean trajectory error = ' num2str(NMTE) '%'])
```

```text:Output
Normalized mean trajectory error = 7.1901%
```

```matlab:Code

% We plot the true test set trajectory in the reduced coordinates and compare it to 
% the prediction. 
plotReducedCoordinates(yDataTrunc(indTest,:), etaRec(indTest,:))
legend({'Test set (truncated)', 'Prediction'})
if size(Ae,1)==2
    % Plot SSM with trajectories in the normal form reduced coordinates
    plotSSMandTrajectories(IMInfo, outdof, xDataTrunc(indTest,:), ...
        zRec(indTest,:), 'NFT', RDInfo.transformation.map)
    view(-100,20); legend('off')
else
    view(3)
end
```

![figure_17.png
](README_images/figure_17.png
)

We plot the model predictions in physical coordinates. The reduced model seems to do well on previously unseen data, provided that it is close to the manifold.

```matlab:Code
fig = customFigure('subPlot',[3 1]); cc = 0;
for iTraj = [indTrain(1:end-1)+1]
    cc = cc + 1;
    subplot(3,1,cc)
    plot(yDataTrunc{iTraj,1},yDataTrunc{iTraj,2}(1,:),'k','Linewidth',1,'Color',colors(colFOM,:))
    plot(etaRec{iTraj,1},etaRec{iTraj,2}(1,:),':','Linewidth',1.5,'Color',colors(colSSML,:))
    xlabel('$t \, [$s$]$','Interpreter','latex'); 
    ylabel('$u_{1}$ [m]','interpreter','latex'); xlim([0 8])
end
legend('Training trajectory','Prediction')
```

![figure_18.png
](README_images/figure_18.png
)

```matlab:Code
fig = customFigure('subPlot',[3 1]); cc = 0;
for iTraj = [indTrain(1:end-1)+1]  
    cc = cc + 1;
    subplot(3,1,cc)
    plot(yDataTrunc{iTraj,1},yDataTrunc{iTraj,2}(2,:),'k','Linewidth',1,'Color',colors(colFOM,:))
    plot(etaRec{iTraj,1},etaRec{iTraj,2}(2,:),':','Linewidth',1.5,'Color',colors(colSSML,:))
    xlabel('$t \, [$s$]$','Interpreter','latex'); 
    ylabel('$u_{2}$ [m]','interpreter','latex'); xlim([0 8])
end
legend('Training trajectory','Prediction')
```

![figure_19.png
](README_images/figure_19.png
)

```matlab:Code
fig = customFigure('subPlot',[3 1]); cc = 0;
for iTraj = [indTrain(1:end-1)+1]
    cc = cc + 1;
    subplot(3,1,cc)
    plot(xDataTrunc{iTraj,1},xDataTrunc{iTraj,2}(outdof(1),:),'k','Linewidth',1,'Color',colors(colFOM,:))
    plot(yRec{iTraj,1},yRec{iTraj,2}(outdof(1),:),':','Linewidth',1.5,'Color',colors(colSSML,:))
    xlabel('$t \, [$s$]$','Interpreter','latex'); 
    ylabel('$q_{A}$ [m]','interpreter','latex'); xlim([0 8])
end
legend('Training trajectory','Prediction')
```

![figure_20.png
](README_images/figure_20.png
)

```matlab:Code
fig = customFigure('subPlot',[3 1]); cc = 0;
for iTraj = [indTrain(1:end-1)+1]  
    cc = cc + 1;
    subplot(3,1,cc)
    plot(xDataTrunc{iTraj,1},xDataTrunc{iTraj,2}(outdof(2),:),'k','Linewidth',1,'Color',colors(colFOM,:))
    plot(yRec{iTraj,1},yRec{iTraj,2}(outdof(2),:),':','Linewidth',1.5,'Color',colors(colSSML,:))
    xlabel('$t \, [$s$]$','Interpreter','latex'); 
    ylabel('$q_{B}$ [m]','interpreter','latex'); xlim([0 8])
end
legend('Training trajectory','Prediction')
```

![figure_21.png
](README_images/figure_21.png
)

We also check the goodness of fit by evaluating the 1-step prediction error, defined as

> ${{{\mathrm{P}\mathrm{E}}}_1^{{\mathrm{t}\mathrm{e}\mathrm{s}\mathrm{t}}} (t)=\frac{100}{N_{{\mathrm{t}\mathrm{e}\mathrm{s}\mathrm{t}}} \|y \|}\sum_{k=1}^{N_{{\mathrm{t}\mathrm{e}\mathrm{s}\mathrm{t}}} } \|y_k (t)-\hat{y} (\Delta t;y_k (t-\Delta t))\|}$,

where $N_{{\mathrm{t}\mathrm{e}\mathrm{s}\mathrm{t}}}$ is the number of test trajectories, $y_k (t)$ is the $k$-th testing trajectory in reduced coordinates, $\hat{y} (\Delta t;y_k (t-\Delta t))$ is the reduced dynamics advection of a sampling time interval $\Delta t$ from the initial condition $y_k (t-\Delta t)$. We note that, compared to purely linear reduced dynamics, the nonlinear model reduces the 1-step prediction error of one order of magnitude.

```matlab:Code
step1error;
```

![figure_22.png
](README_images/figure_22.png
)

# Adding forcing to the ROM

Outer directions: either consider all outer modes or a subset

```matlab:Code
nModes = size(V,2);
outerModes = setdiff(1:nModes,masterModes);
Vo = V(:,outerModes); 
Wo = W(outerModes,:); Lo = full(Wo*A*Vo);
```

Forcing vector

```matlab:Code
forcingVectors = [zeros(n,1); M\loadVector];
```

Construct time periodic SSM model

```matlab:Code
[IMInfoF,RDInfoF] = forcedSSMROM(IMInfo,RDInfo,'nForcingFrequencies',1,...
         'forcingVectors',forcingVectors,'We',We,'Lo',Lo,'Vo',Vo, 'Wo',Wo);
```

# Generate Frequency Responses

We compute them also with SSMTool in order to compare the results (see the papers above for additional validations and comparisons).

```matlab:Code
epsilon = [0.02 0.07];
mFreqs = [1 2];
resonantModes = [1 2 3 4];
omegaSpan = [0.92 1.07]*imag(lambda(1));
% SSMLearn
[FRCSSMLearn] = continuationFRCep(IMInfoF, RDInfoF, epsilon, omegaSpan,@(x) x(outdof,:), mFreqs,resonantModes, 'SSMLearnFRC');
```

```text:Output
 Run='SSMLearnFRC.ep': Continue equilibria along primary branch.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          3.83e-01  2.11e+02    0.0    0.0    0.0
   1   1  1.00e+00  1.16e+00  4.28e-02  2.11e+02    0.0    0.0    0.0
   2   1  1.00e+00  4.82e-02  1.76e-03  2.11e+02    0.0    0.0    0.0
   3   1  1.00e+00  2.10e-03  8.06e-08  2.11e+02    0.0    0.1    0.0
   4   1  1.00e+00  2.07e-07  5.37e-14  2.11e+02    0.0    0.1    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   2.1111e+02      1  EP      1.4922e+02   1.8343e-02   2.2716e-02   2.7527e+00   2.8531e+00   2.0000e-02
    1  00:00:00   2.1108e+02      2  HB      1.4920e+02   1.8611e-02   2.2644e-02   2.7246e+00   2.8534e+00   2.0000e-02
   10  00:00:00   2.1072e+02      3          1.4895e+02   2.2752e-02   2.2534e-02   2.5629e+00   2.9487e+00   2.0000e-02
   20  00:00:00   2.0985e+02      4          1.4831e+02   3.5560e-02   2.6254e-02   2.7519e+00   3.6570e+00   2.0000e-02
   27  00:00:01   2.0925e+02      5  SN      1.4782e+02   4.2027e-02   2.5641e-02   3.5620e+00   5.3611e+00   2.0000e-02
   27  00:00:01   2.0925e+02      6  FP      1.4782e+02   4.2026e-02   2.5640e-02   3.5621e+00   5.3613e+00   2.0000e-02
   30  00:00:01   2.0946e+02      7          1.4792e+02   3.2116e-02   1.5783e-02   4.0766e+00   6.3797e+00   2.0000e-02
   31  00:00:01   2.0951e+02      8  FP      1.4793e+02   2.7737e-02   1.1897e-02   4.2192e+00   6.6628e+00   2.0000e-02
   31  00:00:01   2.0951e+02      9  SN      1.4793e+02   2.7707e-02   1.1870e-02   4.2201e+00   6.6646e+00   2.0000e-02
   40  00:00:01   2.0913e+02     10          1.4763e+02   1.5980e-02   3.2760e-03   4.4943e+00   7.2481e+00   2.0000e-02
   50  00:00:01   2.0782e+02     11          1.4669e+02   9.3748e-03   7.3938e-04   4.5942e+00   7.5062e+00   2.0000e-02
   60  00:00:01   2.0606e+02     12          1.4543e+02   6.2255e-03   2.2325e-04   4.6354e+00   7.6234e+00   2.0000e-02
   70  00:00:01   2.0387e+02     13          1.4387e+02   4.4116e-03   8.0581e-05   4.6586e+00   7.6909e+00   2.0000e-02
   80  00:00:01   2.0122e+02     14          1.4200e+02   3.2674e-03   3.3012e-05   4.6731e+00   7.7336e+00   2.0000e-02
   90  00:00:01   1.9763e+02     15          1.3945e+02   2.4154e-03   1.3418e-05   4.6839e+00   7.7655e+00   2.0000e-02
   97  00:00:02   1.9458e+02     16  EP      1.3729e+02   1.9780e-03   7.3905e-06   4.6894e+00   7.7820e+00   2.0000e-02

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:02   2.1111e+02     17  EP      1.4922e+02   1.8343e-02   2.2716e-02   2.7527e+00   2.8531e+00   2.0000e-02
    6  00:00:02   2.1151e+02     18  HB      1.4950e+02   1.8250e-02   2.3812e-02   3.2144e+00   2.9204e+00   2.0000e-02
   10  00:00:02   2.1181e+02     19          1.4970e+02   2.1725e-02   2.4537e-02   3.4275e+00   2.8879e+00   2.0000e-02
   20  00:00:02   2.1218e+02     20          1.4997e+02   2.7325e-02   2.5882e-02   3.4501e+00   2.6563e+00   2.0000e-02
   30  00:00:02   2.1296e+02     21          1.5056e+02   3.6380e-02   2.4707e-02   2.6465e+00   8.2567e-01   2.0000e-02
   31  00:00:02   2.1297e+02     22  SN      1.5057e+02   3.5178e-02   2.3051e-02   2.5376e+00   6.0728e-01   2.0000e-02
   31  00:00:02   2.1297e+02     23  FP      1.5057e+02   3.5176e-02   2.3048e-02   2.5374e+00   6.0688e-01   2.0000e-02
   34  00:00:02   2.1293e+02     24  SN      1.5055e+02   2.7960e-02   1.4853e-02   2.1555e+00  -1.5187e-01   2.0000e-02
   34  00:00:02   2.1293e+02     25  FP      1.5055e+02   2.7913e-02   1.4805e-02   2.1537e+00  -1.5559e-01   2.0000e-02
   40  00:00:02   2.1325e+02     26          1.5078e+02   1.6811e-02   4.5493e-03   1.8185e+00  -8.6492e-01   2.0000e-02
   50  00:00:02   2.1439e+02     27          1.5158e+02   1.0126e-02   1.0607e-03   1.7039e+00  -1.1703e+00   2.0000e-02
   60  00:00:02   2.1609e+02     28          1.5279e+02   6.6264e-03   2.9448e-04   1.6570e+00  -1.3121e+00   2.0000e-02
   70  00:00:03   2.1822e+02     29          1.5429e+02   4.6568e-03   1.0103e-04   1.6317e+00  -1.3894e+00   2.0000e-02
   80  00:00:03   2.2078e+02     30          1.5610e+02   3.4313e-03   4.0083e-05   1.6162e+00  -1.4369e+00   2.0000e-02
   90  00:00:03   2.2418e+02     31          1.5850e+02   2.5437e-03   1.6228e-05   1.6049e+00  -1.4711e+00   2.0000e-02
   94  00:00:03   2.2583e+02     32  EP      1.5967e+02   2.2603e-03   1.1364e-05   1.6014e+00  -1.4820e+00   2.0000e-02
```

![figure_23.png
](README_images/figure_23.png
)

![figure_24.png
](README_images/figure_24.png
)

```text:Output
 Run='SSMLearnFRC.ep.ep': Continue equilibria along primary branch.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          2.90e+00  2.11e+02    0.0    0.0    0.0
   1   1  5.56e-01  1.45e+00  1.55e+00  2.11e+02    0.0    0.0    0.0
   2   1  1.00e+00  3.47e-01  3.69e-01  2.11e+02    0.0    0.0    0.0
   3   1  1.00e+00  4.28e-01  8.54e-03  2.11e+02    0.0    0.0    0.0
   4   1  1.00e+00  2.10e-02  2.43e-04  2.11e+02    0.0    0.0    0.0
   5   1  1.00e+00  2.03e-04  2.39e-09  2.11e+02    0.0    0.0    0.0
   6   1  1.00e+00  5.22e-10  4.23e-16  2.11e+02    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   2.1111e+02      1  EP      1.4922e+02   2.8755e-02   5.7950e-02   2.7578e+00   2.7856e+00   7.0000e-02
   10  00:00:00   2.1041e+02      2          1.4875e+02   3.9200e-02   5.1795e-02   2.2443e+00   2.4604e+00   7.0000e-02
   12  00:00:00   2.1033e+02      3  HB      1.4869e+02   4.0585e-02   5.1558e-02   2.2244e+00   2.4552e+00   7.0000e-02
   20  00:00:00   2.0849e+02      4          1.4738e+02   7.1446e-02   5.9413e-02   2.2215e+00   2.7182e+00   7.0000e-02
   30  00:00:00   2.0447e+02      5          1.4449e+02   1.4061e-01   9.6554e-02   2.9535e+00   4.2710e+00   7.0000e-02
   37  00:00:00   2.0416e+02      6  SN      1.4424e+02   1.4483e-01   9.7485e-02   3.2565e+00   4.8802e+00   7.0000e-02
   37  00:00:00   2.0416e+02      7  FP      1.4424e+02   1.4483e-01   9.7485e-02   3.2565e+00   4.8802e+00   7.0000e-02
   40  00:00:00   2.0421e+02      8          1.4427e+02   1.4316e-01   9.5709e-02   3.3651e+00   5.0969e+00   7.0000e-02
   50  00:00:00   2.0587e+02      9          1.4538e+02   1.0654e-01   6.6614e-02   3.9683e+00   6.2877e+00   7.0000e-02
   60  00:00:00   2.0773e+02     10          1.4664e+02   4.8206e-02   1.9331e-02   4.4842e+00   7.2874e+00   7.0000e-02
   61  00:00:00   2.0773e+02     11  FP      1.4664e+02   4.7767e-02   1.8979e-02   4.4874e+00   7.2937e+00   7.0000e-02
   61  00:00:00   2.0773e+02     12  SN      1.4664e+02   4.7767e-02   1.8979e-02   4.4874e+00   7.2937e+00   7.0000e-02
   70  00:00:01   2.0768e+02     13          1.4659e+02   4.1374e-02   1.3996e-02   4.5310e+00   7.3828e+00   7.0000e-02
   80  00:00:01   2.0706e+02     14          1.4615e+02   2.9549e-02   6.1493e-03   4.5985e+00   7.5328e+00   7.0000e-02
   90  00:00:01   2.0532e+02     15          1.4490e+02   1.9521e-02   1.9398e-03   4.6425e+00   7.6463e+00   7.0000e-02
  100  00:00:01   2.0303e+02     16          1.4328e+02   1.3987e-02   7.3207e-04   4.6636e+00   7.7061e+00   7.0000e-02
  110  00:00:01   2.0018e+02     17          1.4126e+02   1.0390e-02   3.0343e-04   4.6768e+00   7.7447e+00   7.0000e-02
  120  00:00:01   1.9617e+02     18          1.3842e+02   7.6526e-03   1.2204e-04   4.6867e+00   7.7741e+00   7.0000e-02
  124  00:00:01   1.9458e+02     19  EP      1.3729e+02   6.9256e-03   9.0609e-05   4.6894e+00   7.7819e+00   7.0000e-02

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:01   2.1111e+02     20  EP      1.4922e+02   2.8755e-02   5.7950e-02   2.7578e+00   2.7856e+00   7.0000e-02
   10  00:00:01   2.1211e+02     21          1.4989e+02   3.8986e-02   5.6405e-02   3.9112e+00   3.6151e+00   7.0000e-02
   12  00:00:01   2.1218e+02     22  HB      1.4994e+02   4.0262e-02   5.6272e-02   3.9341e+00   3.6253e+00   7.0000e-02
   20  00:00:01   2.1308e+02     23          1.5058e+02   5.7456e-02   5.9098e-02   4.0163e+00   3.5550e+00   7.0000e-02
   30  00:00:02   2.1621e+02     24          1.5283e+02   1.2211e-01   9.3709e-02   3.3346e+00   2.0377e+00   7.0000e-02
   36  00:00:02   2.1651e+02     25  SN      1.5306e+02   1.2605e-01   9.3767e-02   2.9570e+00   1.2774e+00   7.0000e-02
   36  00:00:02   2.1651e+02     26  FP      1.5306e+02   1.2605e-01   9.3767e-02   2.9570e+00   1.2774e+00   7.0000e-02
   40  00:00:02   2.1636e+02     27          1.5296e+02   1.1905e-01   8.6392e-02   2.6972e+00   7.6037e-01   7.0000e-02
   50  00:00:02   2.1469e+02     28          1.5179e+02   5.0912e-02   2.4258e-02   1.8440e+00  -9.0324e-01   7.0000e-02
   52  00:00:02   2.1468e+02     29  SN      1.5179e+02   4.8704e-02   2.2252e-02   1.8252e+00  -9.4046e-01   7.0000e-02
   52  00:00:02   2.1468e+02     30  FP      1.5179e+02   4.8704e-02   2.2252e-02   1.8252e+00  -9.4046e-01   7.0000e-02
   60  00:00:02   2.1475e+02     31          1.5184e+02   4.1361e-02   1.5776e-02   1.7671e+00  -1.0588e+00   7.0000e-02
   70  00:00:02   2.1537e+02     32          1.5228e+02   2.9615e-02   6.8929e-03   1.6917e+00  -1.2274e+00   7.0000e-02
   80  00:00:02   2.1714e+02     33          1.5352e+02   1.9557e-02   2.1096e-03   1.6452e+00  -1.3513e+00   7.0000e-02
   90  00:00:02   2.1942e+02     34          1.5514e+02   1.4021e-02   7.8062e-04   1.6238e+00  -1.4143e+00   7.0000e-02
  100  00:00:03   2.2228e+02     35          1.5716e+02   1.0420e-02   3.1922e-04   1.6105e+00  -1.4543e+00   7.0000e-02
  110  00:00:03   2.2583e+02     36  EP      1.5967e+02   7.9143e-03   1.3931e-04   1.6014e+00  -1.4819e+00   7.0000e-02
```

![figure_25.png
](README_images/figure_25.png
)

![figure_26.png
](README_images/figure_26.png
)

```matlab:Code
% SSMTool
[FRCSSMTool1] = SSMToolFRCFE(M,C,K,fnl,f_0,outdof,epsilon(1),resonantModes,4,omegaSpan,mFreqs,'SSMToolFRC');
```

```text:Output
Due to high-dimensionality, we compute only the first 5 eigenvalues with the smallest magnitude. These would also be used to compute the spectral quotients
Assuming a proportional damping hypthesis with symmetric matrices
modal damping ratio for 1 mode is 2.000000e-03
modal damping ratio for 2 mode is 2.000000e-03
modal damping ratio for 3 mode is 2.102524e-03
modal damping ratio for 4 mode is 2.141338e-03
modal damping ratio for 5 mode is 2.369557e-03
the left eigenvectors may be incorrect in case of asymmetry of matrices

 The first 10 nonzero eigenvalues are given as 
   1.0e+02 *

  -0.0030 + 1.4922i
  -0.0030 - 1.4922i
  -0.0060 + 2.9878i
  -0.0060 - 2.9878i
  -0.0071 + 3.3973i
  -0.0071 - 3.3973i
  -0.0076 + 3.5356i
  -0.0076 - 3.5356i
  -0.0101 + 4.2617i
  -0.0101 - 4.2617i

(near) outer resonance detected for the following combination of master eigenvalues
     2     0     0     0
     0     0     2     1
     1     1     1     0
     0     2     0     0
     0     0     1     2
     1     1     0     1
     2     0     0     0
     0     0     2     1
     1     1     1     0
     0     2     0     0
     0     0     1     2
     1     1     0     1
     1     0     1     0
     0     1     2     0
     3     0     0     0
     0     1     0     1
     0     3     0     0
     1     0     0     2

These are in resonance with the follwing eigenvalues of the slave subspace
   1.0e+02 *

  -0.0071 + 3.3973i
  -0.0071 + 3.3973i
  -0.0071 + 3.3973i
  -0.0071 - 3.3973i
  -0.0071 - 3.3973i
  -0.0071 - 3.3973i
  -0.0076 + 3.5356i
  -0.0076 + 3.5356i
  -0.0076 + 3.5356i
  -0.0076 - 3.5356i
  -0.0076 - 3.5356i
  -0.0076 - 3.5356i
  -0.0101 + 4.2617i
  -0.0101 + 4.2617i
  -0.0101 + 4.2617i
  -0.0101 - 4.2617i
  -0.0101 - 4.2617i
  -0.0101 - 4.2617i

sigma_out = 3
(near) inner resonance detected for the following combination of master eigenvalues
     0     1     1     0
     1     0     1     1
     2     1     0     0
     1     0     0     1
     0     1     1     1
     1     2     0     0
     2     0     0     0
     0     0     2     1
     1     1     1     0
     0     2     0     0
     0     0     1     2
     1     1     0     1

These are in resonance with the follwing eigenvalues of the master subspace
   1.0e+02 *

  -0.0030 + 1.4922i
  -0.0030 + 1.4922i
  -0.0030 + 1.4922i
  -0.0030 - 1.4922i
  -0.0030 - 1.4922i
  -0.0030 - 1.4922i
  -0.0060 + 2.9878i
  -0.0060 + 2.9878i
  -0.0060 + 2.9878i
  -0.0060 - 2.9878i
  -0.0060 - 2.9878i
  -0.0060 - 2.9878i

sigma_in = 3
Due to (near) outer resonance, the exisitence of the manifold is questionable and the underlying computation may suffer.
Attempting manifold computation
Manifold computation time at order 2 = 00:00:02
Estimated memory usage at order  2 = 1.61E+01 MB
Manifold computation time at order 3 = 00:00:03
Estimated memory usage at order  3 = 2.10E+01 MB
Manifold computation time at order 4 = 00:00:08
Estimated memory usage at order  4 = 3.26E+01 MB

 Run='SSMToolFRC.ep': Continue equilibria along primary branch.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          5.49e-01  2.11e+02    0.0    0.0    0.0
   1   1  1.00e+00  4.60e-01  7.51e-02  2.11e+02    0.0    0.0    0.0
   2   1  1.00e+00  8.03e-02  1.21e-02  2.11e+02    0.0    0.0    0.0
   3   1  1.00e+00  1.63e-02  1.08e-04  2.11e+02    0.0    0.0    0.0
   4   1  1.00e+00  3.78e-04  1.70e-08  2.11e+02    0.0    0.0    0.0
   5   1  1.00e+00  1.94e-08  2.38e-16  2.11e+02    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   2.1113e+02      1  EP      1.4922e+02   3.3933e-03   2.3719e-03   1.1624e+00   4.4064e+00   2.0000e-02
    1  00:00:00   2.1110e+02      2  HB      1.4920e+02   3.4421e-03   2.3640e-03   1.1359e+00   4.4068e+00   2.0000e-02
   10  00:00:00   2.1074e+02      3          1.4895e+02   4.2261e-03   2.3469e-03   9.7590e-01   4.5042e+00   2.0000e-02
   20  00:00:00   2.0984e+02      4          1.4828e+02   6.7184e-03   2.7463e-03   1.1819e+00   5.2470e+00   2.0000e-02
   27  00:00:00   2.0928e+02      5  SN      1.4781e+02   7.8773e-03   2.6854e-03   1.9696e+00   6.8988e+00   2.0000e-02
   27  00:00:00   2.0928e+02      6  FP      1.4781e+02   7.8772e-03   2.6854e-03   1.9696e+00   6.8989e+00   2.0000e-02
   30  00:00:00   2.0951e+02      7          1.4791e+02   5.9507e-03   1.6383e-03   2.5024e+00   7.9505e+00   2.0000e-02
   31  00:00:00   2.0956e+02      8  FP      1.4793e+02   5.1002e-03   1.2215e-03   2.6499e+00   8.2419e+00   2.0000e-02
   31  00:00:00   2.0956e+02      9  SN      1.4793e+02   5.0957e-03   1.2194e-03   2.6506e+00   8.2434e+00   2.0000e-02
   40  00:00:00   2.0920e+02     10          1.4764e+02   2.9728e-03   3.4923e-04   2.9186e+00   8.8107e+00   2.0000e-02
   50  00:00:00   2.0791e+02     11          1.4670e+02   1.7392e-03   7.8631e-05   3.0206e+00   9.0722e+00   2.0000e-02
   60  00:00:00   2.0614e+02     12          1.4544e+02   1.1523e-03   2.3607e-05   3.0623e+00   9.1907e+00   2.0000e-02
   70  00:00:01   2.0396e+02     13          1.4389e+02   8.1610e-04   8.5094e-06   3.0856e+00   9.2585e+00   2.0000e-02
   80  00:00:01   2.0132e+02     14          1.4202e+02   6.0429e-04   3.4840e-06   3.1002e+00   9.3014e+00   2.0000e-02
   90  00:00:01   1.9774e+02     15          1.3948e+02   4.4672e-04   1.4161e-06   3.1110e+00   9.3334e+00   2.0000e-02
   97  00:00:01   1.9465e+02     16  EP      1.3729e+02   3.6478e-04   7.7348e-07   3.1166e+00   9.3501e+00   2.0000e-02

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:01   2.1113e+02     17  EP      1.4922e+02   3.3933e-03   2.3719e-03   1.1624e+00   4.4064e+00   2.0000e-02
    5  00:00:01   2.1153e+02     18  HB      1.4950e+02   3.3527e-03   2.4974e-03   1.6228e+00   4.4714e+00   2.0000e-02
   10  00:00:01   2.1180e+02     19          1.4969e+02   3.9243e-03   2.5742e-03   1.8299e+00   4.4456e+00   2.0000e-02
   20  00:00:01   2.1213e+02     20          1.4993e+02   4.8543e-03   2.7052e-03   1.8702e+00   4.2468e+00   2.0000e-02
   30  00:00:01   2.1292e+02     21          1.5053e+02   6.7610e-03   2.7518e-03   1.1888e+00   2.6337e+00   2.0000e-02
   32  00:00:01   2.1294e+02     22  SN      1.5055e+02   6.3508e-03   2.3818e-03   9.3791e-01   2.1273e+00   2.0000e-02
   32  00:00:01   2.1294e+02     23  FP      1.5055e+02   6.3505e-03   2.3815e-03   9.3776e-01   2.1270e+00   2.0000e-02
   34  00:00:02   2.1291e+02     24  FP      1.5054e+02   5.1970e-03   1.6084e-03   5.9856e-01   1.4514e+00   2.0000e-02
   34  00:00:02   2.1291e+02     25  SN      1.5054e+02   5.1911e-03   1.6048e-03   5.9722e-01   1.4487e+00   2.0000e-02
   40  00:00:02   2.1320e+02     26          1.5075e+02   3.1706e-03   5.0883e-04   2.5427e-01   7.2470e-01   2.0000e-02
   50  00:00:02   2.1427e+02     27          1.5151e+02   1.9232e-03   1.2144e-04   1.3548e-01   4.1072e-01   2.0000e-02
   60  00:00:02   2.1595e+02     28          1.5270e+02   1.2519e-03   3.3180e-05   8.6284e-02   2.6254e-01   2.0000e-02
   70  00:00:02   2.1805e+02     29          1.5419e+02   8.7671e-04   1.1257e-05   6.0176e-02   1.8256e-01   2.0000e-02
   80  00:00:02   2.2059e+02     30          1.5598e+02   6.4446e-04   4.4331e-06   4.4187e-02   1.3370e-01   2.0000e-02
   90  00:00:02   2.2391e+02     31          1.5833e+02   4.7812e-04   1.7989e-06   3.2770e-02   9.8950e-02   2.0000e-02
   95  00:00:02   2.2581e+02     32  EP      1.5967e+02   4.1685e-04   1.1893e-06   2.8567e-02   8.6194e-02   2.0000e-02
```

![figure_27.png
](README_images/figure_27.png
)

![figure_28.png
](README_images/figure_28.png
)

```text:Output
the forcing frequency 1.3729e+02 is nearly resonant with the eigenvalue -2.9845e-01 + i1.4922e+02
...
the forcing frequency 1.5967e+02 is nearly resonant with the eigenvalue -2.9845e-01 + i1.4922e+02
```

![figure_29.png
](README_images/figure_29.png
)

![figure_30.png
](README_images/figure_30.png
)

```matlab:Code
[FRCSSMTool2] = SSMToolFRCFE(M,C,K,fnl,f_0,outdof,epsilon(2),resonantModes,4,omegaSpan,mFreqs,'SSMToolFRC');
```

```text:Output
Due to high-dimensionality, we compute only the first 5 eigenvalues with the smallest magnitude. These would also be used to compute the spectral quotients
Assuming a proportional damping hypthesis with symmetric matrices
modal damping ratio for 1 mode is 2.000000e-03
modal damping ratio for 2 mode is 2.000000e-03
modal damping ratio for 3 mode is 2.102524e-03
modal damping ratio for 4 mode is 2.141338e-03
modal damping ratio for 5 mode is 2.369557e-03
the left eigenvectors may be incorrect in case of asymmetry of matrices

 The first 10 nonzero eigenvalues are given as 
   1.0e+02 *

  -0.0030 + 1.4922i
  -0.0030 - 1.4922i
  -0.0060 + 2.9878i
  -0.0060 - 2.9878i
  -0.0071 + 3.3973i
  -0.0071 - 3.3973i
  -0.0076 + 3.5356i
  -0.0076 - 3.5356i
  -0.0101 + 4.2617i
  -0.0101 - 4.2617i

(near) outer resonance detected for the following combination of master eigenvalues
     2     0     0     0
     0     0     2     1
     1     1     1     0
     0     2     0     0
     0     0     1     2
     1     1     0     1
     2     0     0     0
     0     0     2     1
     1     1     1     0
     0     2     0     0
     0     0     1     2
     1     1     0     1
     1     0     1     0
     0     1     2     0
     3     0     0     0
     0     1     0     1
     0     3     0     0
     1     0     0     2

These are in resonance with the follwing eigenvalues of the slave subspace
   1.0e+02 *

  -0.0071 + 3.3973i
  -0.0071 + 3.3973i
  -0.0071 + 3.3973i
  -0.0071 - 3.3973i
  -0.0071 - 3.3973i
  -0.0071 - 3.3973i
  -0.0076 + 3.5356i
  -0.0076 + 3.5356i
  -0.0076 + 3.5356i
  -0.0076 - 3.5356i
  -0.0076 - 3.5356i
  -0.0076 - 3.5356i
  -0.0101 + 4.2617i
  -0.0101 + 4.2617i
  -0.0101 + 4.2617i
  -0.0101 - 4.2617i
  -0.0101 - 4.2617i
  -0.0101 - 4.2617i

sigma_out = 3
(near) inner resonance detected for the following combination of master eigenvalues
     0     1     1     0
     1     0     1     1
     2     1     0     0
     1     0     0     1
     0     1     1     1
     1     2     0     0
     2     0     0     0
     0     0     2     1
     1     1     1     0
     0     2     0     0
     0     0     1     2
     1     1     0     1

These are in resonance with the follwing eigenvalues of the master subspace
   1.0e+02 *

  -0.0030 + 1.4922i
  -0.0030 + 1.4922i
  -0.0030 + 1.4922i
  -0.0030 - 1.4922i
  -0.0030 - 1.4922i
  -0.0030 - 1.4922i
  -0.0060 + 2.9878i
  -0.0060 + 2.9878i
  -0.0060 + 2.9878i
  -0.0060 - 2.9878i
  -0.0060 - 2.9878i
  -0.0060 - 2.9878i

sigma_in = 3
Due to (near) outer resonance, the exisitence of the manifold is questionable and the underlying computation may suffer.
Attempting manifold computation
Manifold computation time at order 2 = 00:00:01
Estimated memory usage at order  2 = 1.61E+01 MB
Manifold computation time at order 3 = 00:00:02
Estimated memory usage at order  3 = 2.10E+01 MB
Manifold computation time at order 4 = 00:00:08
Estimated memory usage at order  4 = 3.26E+01 MB

 Run='SSMToolFRC.ep': Continue equilibria along primary branch.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          3.82e+00  2.11e+02    0.0    0.0    0.0
   1   1  1.00e+00  7.34e-01  1.29e+00  2.11e+02    0.0    0.0    0.0
   2   1  1.00e+00  1.02e+00  7.34e-01  2.11e+02    0.0    0.0    0.0
   3   1  5.40e-01  2.17e+00  3.94e-01  2.11e+02    0.0    0.0    0.0
   4   1  1.00e+00  6.52e-01  4.40e-02  2.11e+02    0.0    0.0    0.0
   5   1  1.00e+00  1.43e-01  2.44e-03  2.11e+02    0.0    0.0    0.0
   6   1  1.00e+00  1.21e-03  1.30e-06  2.11e+02    0.0    0.0    0.0
   7   1  1.00e+00  3.00e-06  1.15e-12  2.11e+02    0.0    0.0    0.0
   8   1  1.00e+00  4.70e-13  8.54e-17  2.11e+02    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   2.1113e+02      1  EP      1.4922e+02   5.3502e-03   6.0161e-03   1.1237e+00   4.2968e+00   7.0000e-02
   10  00:00:00   2.1045e+02      2          1.4875e+02   7.2896e-03   5.3427e-03   6.4656e-01   3.9989e+00   7.0000e-02
   12  00:00:00   2.1036e+02      3  HB      1.4869e+02   7.5694e-03   5.3122e-03   6.2580e-01   3.9938e+00   7.0000e-02
   20  00:00:00   2.0836e+02      4          1.4727e+02   1.3802e-02   6.1589e-03   6.3466e-01   4.2799e+00   7.0000e-02
   30  00:00:00   2.0410e+02      5          1.4420e+02   2.7433e-02   9.8985e-03   1.3263e+00   5.7521e+00   7.0000e-02
   38  00:00:00   2.0367e+02      6  SN      1.4386e+02   2.8611e-02   1.0092e-02   1.6541e+00   6.4119e+00   7.0000e-02
   38  00:00:00   2.0367e+02      7  FP      1.4386e+02   2.8611e-02   1.0092e-02   1.6541e+00   6.4119e+00   7.0000e-02
   40  00:00:00   2.0369e+02      8          1.4387e+02   2.8455e-02   1.0007e-02   1.7187e+00   6.5408e+00   7.0000e-02
   50  00:00:00   2.0479e+02      9          1.4460e+02   2.4097e-02   8.3022e-03   2.1732e+00   7.4385e+00   7.0000e-02
   60  00:00:00   2.0776e+02     10          1.4662e+02   9.0585e-03   2.0526e-03   2.9044e+00   8.8465e+00   7.0000e-02
   62  00:00:00   2.0777e+02     11  FP      1.4662e+02   8.7863e-03   1.9331e-03   2.9149e+00   8.8672e+00   7.0000e-02
   62  00:00:00   2.0777e+02     12  SN      1.4662e+02   8.7863e-03   1.9331e-03   2.9149e+00   8.8672e+00   7.0000e-02
   70  00:00:00   2.0773e+02     13          1.4658e+02   7.7410e-03   1.4856e-03   2.9534e+00   8.9447e+00   7.0000e-02
   80  00:00:01   2.0727e+02     14          1.4625e+02   5.7767e-03   7.4189e-04   3.0160e+00   9.0800e+00   7.0000e-02
   90  00:00:01   2.0554e+02     15          1.4502e+02   3.7106e-03   2.2076e-04   3.0671e+00   9.2079e+00   7.0000e-02
  100  00:00:01   2.0330e+02     16          1.4342e+02   2.6441e-03   8.2301e-05   3.0894e+00   9.2705e+00   7.0000e-02
  110  00:00:01   2.0052e+02     17          1.4145e+02   1.9642e-03   3.4157e-05   3.1030e+00   9.3101e+00   7.0000e-02
  120  00:00:01   1.9663e+02     18          1.3869e+02   1.4479e-03   1.3782e-05   3.1132e+00   9.3401e+00   7.0000e-02
  125  00:00:01   1.9465e+02     19  EP      1.3729e+02   1.2773e-03   9.4836e-06   3.1166e+00   9.3501e+00   7.0000e-02

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:01   2.1113e+02     20  EP      1.4922e+02   5.3502e-03   6.0161e-03   1.1237e+00   4.2968e+00   7.0000e-02
   10  00:00:01   2.1215e+02     21          1.4991e+02   7.1561e-03   6.0057e-03   2.3200e+00   5.1556e+00   7.0000e-02
   11  00:00:01   2.1220e+02     22  HB      1.4994e+02   7.3280e-03   5.9990e-03   2.3374e+00   5.1634e+00   7.0000e-02
   20  00:00:01   2.1309e+02     23          1.5057e+02   1.0418e-02   6.3430e-03   2.4254e+00   5.0963e+00   7.0000e-02
   30  00:00:01   2.1601e+02     24          1.5269e+02   2.1375e-02   1.0049e-02   1.7629e+00   3.6130e+00   7.0000e-02
   36  00:00:01   2.1630e+02     25  SN      1.5291e+02   2.2072e-02   1.0037e-02   1.3740e+00   2.8295e+00   7.0000e-02
   36  00:00:01   2.1630e+02     26  FP      1.5291e+02   2.2072e-02   1.0037e-02   1.3740e+00   2.8295e+00   7.0000e-02
   40  00:00:02   2.1612e+02     27          1.5280e+02   2.0652e-02   9.0545e-03   1.0827e+00   2.2501e+00   7.0000e-02
   50  00:00:02   2.1463e+02     28  SN      1.5177e+02   9.0011e-03   2.4033e-03   2.5747e-01   6.3957e-01   7.0000e-02
   50  00:00:02   2.1463e+02     29  FP      1.5177e+02   9.0010e-03   2.4033e-03   2.5747e-01   6.3956e-01   7.0000e-02
   50  00:00:02   2.1464e+02     30          1.5177e+02   8.7598e-03   2.2739e-03   2.4617e-01   6.1693e-01   7.0000e-02
   60  00:00:02   2.1476e+02     31          1.5186e+02   7.2124e-03   1.4808e-03   1.7968e-01   4.7951e-01   7.0000e-02
   70  00:00:02   2.1575e+02     32          1.5256e+02   4.8206e-03   5.1516e-04   1.0160e-01   2.9716e-01   7.0000e-02
   80  00:00:02   2.1768e+02     33          1.5392e+02   3.2793e-03   1.6686e-04   6.5428e-02   1.9670e-01   7.0000e-02
   90  00:00:02   2.2009e+02     34          1.5563e+02   2.3854e-03   6.4151e-05   4.6988e-02   1.4184e-01   7.0000e-02
  100  00:00:02   2.2320e+02     35          1.5783e+02   1.7729e-03   2.6213e-05   3.4787e-02   1.0499e-01   7.0000e-02
  107  00:00:02   2.2581e+02     36  EP      1.5967e+02   1.4594e-03   1.4579e-05   2.8607e-02   8.6276e-02   7.0000e-02
```

![figure_31.png
](README_images/figure_31.png
)

![figure_32.png
](README_images/figure_32.png
)

```text:Output
the forcing frequency 1.3729e+02 is nearly resonant with the eigenvalue -2.9845e-01 + i1.4922e+02
...
the forcing frequency 1.5967e+02 is nearly resonant with the eigenvalue -2.9845e-01 + i1.4922e+02
```

![figure_33.png
](README_images/figure_33.png
)

![figure_34.png
](README_images/figure_34.png
)

```matlab:Code
FRCSSMTool = struct('F1',FRCSSMTool1.F1,'F2',FRCSSMTool2.F1);
```

# **Plot results**

```matlab:Code
fig = customFigure('subPlot',[2 1]);
subplot(211)
plotFRC(FRCSSMLearn, colors(colSSML,:), '- SSMLearn','freqscale',2*pi)
plotFRC(FRCSSMTool, colors(colSSMT,:), '- Tool','freqscale',2*pi)
xlabel('$\Omega$ [Hz]','interpreter','latex')
ylabel('amp($q_{A}$) [m]','interpreter','latex')
legend('off')
xlim([22.5 25])
subplot(212)
plotFRC(FRCSSMLearn, colors(colSSML,:), '- SSMLearn','outamp',2,'freqscale',2*pi)
plotFRC(FRCSSMTool, colors(colSSMT,:), '- SSMTool','outamp',2,'freqscale',2*pi)
xlabel('$\Omega$ [Hz]','interpreter','latex')
ylabel('amp($q_{B}$) [m]','interpreter','latex')
xlim([22.5 25])
ylim([0 11e-4])
legend('NumColumns',2)

% Linear Response
LinResp = @(w) (-M*w^2+1i*w*C+K)\(M*forcingVectors(n+1:end));
omegaVec = linspace(omegaSpan(1),omegaSpan(2),1001);
ampVec = zeros(2,length(omegaVec));
for iW = 1:length(omegaVec)
    X0 = abs(LinResp(omegaVec(iW)));
    ampVec(:,iW) = X0(outdof);
end
for iEps = 1:length(epsilon)
    if iEps == 1; visib = 'on'; else; visib = 'off'; end 
    subplot(211)
    plot(omegaVec/2/pi,ampVec(1,:)*epsilon(iEps),'k:','DisplayName', 'FRC - Linear', 'Linewidth', 1.5,'HandleVisibility',visib)
    subplot(212)
    plot(omegaVec/2/pi,ampVec(2,:)*epsilon(iEps),'k:','DisplayName', 'FRC - Linear', 'Linewidth', 1.5,'HandleVisibility',visib)
end
```

![figure_35.png
](README_images/figure_35.png
)

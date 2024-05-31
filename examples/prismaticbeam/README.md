This is a preview of the livescript `prismaticbeam.mlx`.

# Fitting a 4D SSM for an internally resonant (1:3) prismatic beam
  

See [1] for the details of this model, and [2] for the description of this example.

[1] Li, M., Jain, S., \& Haller, G. (2021). Nonlinear analysis of forced mechanical systems with internal resonance using spectral submanifolds-Part I: Periodic response and forced response curve. *Nonlinear Dynamics* 110, 1005-1043. [DOI: 10.1007/s11071-022-07714-x](https://doi.org/10.1007/s11071-022-07714-x)

[2] Cenedese, M., Marconi, J., Haller, G., \& Jain, S. (2023). Data-assisted non-intrusive model reduction for forced nonlinear finite elements models. Preprint: [arXiv: 2311.17865](https://arxiv.org/abs/2311.17865) 

```matlab:Code
clearvars
close all
format shortg
clc

% Setup colors
colors = colororder; colSSMT = 5; colSSML = 7; colFOM = 1;
```

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
epsilon = 1e-4;
c  = 100;
f = 5/epsilon;
n = 10;               % number of modes
Fext = zeros(n,1);    % excitation at modal coordinate
Fext(1) = f;
[M,C,K,fnl,fext] = buildModel(c,Fext,epsilon,n);
```

```text:Output
Getting nonlinearity coefficients
Loaded coefficients from storage
```

```matlab:Code
outdof = [1];
```

Preliminaries

```matlab:Code
if ~isfile('init.mat')
    [F, lambda] = functionFromTensors(M, C, K, fnl); d_r1 = -real(lambda(1))/abs(lambda(1))*100
    m = 10;
    [W,A,V,lambda] = linearpart(M,C,K,m);
    save('init.mat',"W","A","V","F","lambda")
else
    load('init.mat')
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
            0            0            1            0
            0            0            0            1
      -14.864            0        -0.02            0
            0      -156.07            0        -0.02

```

```matlab:Code
SSMDim = length(masterModes);

% Load and displacement vector: midpoint displacement
displacementVector = zeros(1,n); displacementVector(outdof) = 1;
loadVector = fext;  %  could also be set as modal ones
```

# Compare linear and nonlinear response via modal displacement

We characterize the linear and nonlinear regimes via a static modal analysis, which serves to pick appropriate initial conditions for the trajectories we need to learn the SSM from data.

```matlab:Code
Model = struct('M', M ,'K', K, 'F',F); 
```

Displacement along the first mode

```matlab:Code
iMode = 1; scalingfactor1 = 1e2; nsteps = 50; outdof1 = 1;
[phi1, relativeDiffForceNorm] = modal_analysis(Model,scalingfactor1,nsteps,outdof1,false,iMode);
```

```text:Output
Solving undamped eigenvalue problem
Mode shape
     1
     0
     0
     0
     0
     0
     0
     0
     0
     0

Eigenfrequency
       3.8553
```

![figure_0.png
](README_images/figure_0.png
)

![figure_1.png
](README_images/figure_1.png
)

![figure_2.png
](README_images/figure_2.png
)

![figure_3.png
](README_images/figure_3.png
)

```text:Output
Displacement at output DOF: 100
```

Pick up two initial trajectories that has high expected nonlinear content

```matlab:Code
indIC1 = [nsteps, nsteps-1];
IC1 = [phi1*(scalingfactor1*indIC1/nsteps);zeros(n,length(indIC1))];
```

Displacement along the second mode

```matlab:Code
iMode = 2; scalingfactor2 = 1e2; nsteps = 50; outdof2 = 2;
[phi2, relativeDiffForceNorm2] = modal_analysis(Model,scalingfactor2,nsteps,outdof2,false,iMode);
```

```text:Output
Solving undamped eigenvalue problem
Mode shape
     0
     1
     0
     0
     0
     0
     0
     0
     0
     0

Eigenfrequency
       12.493
```

![figure_4.png
](README_images/figure_4.png
)

![figure_5.png
](README_images/figure_5.png
)

![figure_6.png
](README_images/figure_6.png
)

![figure_7.png
](README_images/figure_7.png
)

```text:Output
Displacement at output DOF: 100
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
indTrain = [1 3 5];
indTest = [2 4 6];
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
newSimulation = false;
observable = @(x) x; % Observe the full phase space
slowTimeScale = 2*pi/abs(lambda(1));
fastTimeScale = 2*pi/abs(lambda(round(SSMDim/2)));
if newSimulation
    % Set integration time to get approximately the desired ampltiude decay
    % to get to the linear regime (guess based on linearized damping)
    numberPeriodsSlow = floor(log(desiredAmplitudeDecay)/...
        (2*pi*(-real(lambda(1))/abs(lambda(1)))))
    endTime = numberPeriodsSlow*slowTimeScale;
    % Set the sampling time to capture approximately 50 points per period 
    % on the faster time scale
    numberPeriodsFast = floor(endTime/fastTimeScale);
    numberPointsPerPeriod = 50;
    nSamp = numberPeriodsFast*numberPointsPerPeriod+1;
    dt = endTime/(nSamp-1);
    % Integrate
    xData = integrateTrajectories(F, endTime, ICs, nSamp, observable);
    DataInfo = struct('n', n, 'loadvector', loadVector);
    save('dataPrismaticDecayModal.mat', 'DataInfo', 'xData', 'dt', 'endTime', 'nSamp')
else
    load dataPrismaticDecayModal.mat
    if n ~= DataInfo.n
       error('The loaded data comes from a model with a different number of elements.') 
    end
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
indPlot1 = indTrain(2);
indPlot2 = indTest(2);

showSpectrogram(xData(indPlot1,:), outdof);
ylim([0,abs(lambda(1))/2/pi*5])
```

![figure_8.png
](README_images/figure_8.png
)

We plot the observables of interest over time for closer inspection. 

```matlab:Code
customFigure();
plot(xData{indPlot1,1}, xData{indPlot1,2}(outdof,:), xData{indPlot2,1}, xData{indPlot2,2}(outdof,:), ':');
xlabel('$t \, [$s$]$','Interpreter','latex'); ylabel('$u \, [$m$]$','Interpreter','latex'); 
legend({'Trajectory 1', 'Trajectory 2'})
title('Generated data')
```

![figure_9.png
](README_images/figure_9.png
)

# Truncate transient data from trajectories

We must however remove the first transient to fulfill the assumption that trajectories lie close to the SSM. We keep only the time interval |sliceInt|.

```matlab:Code
sliceInt = [20*slowTimeScale, endTime];
xDataTrunc = sliceTrajectories(xData, sliceInt);
```

# Datadriven manifold fitting

The measured trajectories are initialized to lie close to the manifold of interest that is tangent at the origin to the eigenspace spanned by the columns of $V_e$. 

As we also know the projection $W_e$ to this eigenspace, we define the modal coordinates as $y=W_e x$. These are the reduced coordinates for our graph style parametrization of the manifold, gauranteed to exists near the origin. We then use the data to learn the nonlinear feature of the manifold geometry, represented via polynomials. Indeed, we seek the $2N\times m_M$ polynomial coefficient matrix $H$ of the manifold on the form

> $x=V_e y+H{{\phi }}_{m,2:M} (y)$,

where the function ${{\phi }}_{m,2:M} (y)$ computes a vector of all $m_M$ monomials from orders 2 up to $M$ of an $m$-vector $y$. From SSM theory, the tangent space of the manifold is $V_e$. The coefficients $H$are obtained via least squares regression.

```matlab:Code
SSMOrder = 3;

% Get projection or modal coordinates 
SSMDim = size(Ve,2);
yDataTrunc = xDataTrunc;
nTraj = size(xDataTrunc,1);
for iTraj = 1:nTraj
    yDataTrunc{iTraj,2} = We*xDataTrunc{iTraj,2};    
end

% Plot reduced coordinates
plotReducedCoordinates(yDataTrunc);
legend({'Test set trajectory', 'Training set trajectory'})
if SSMDim>2
   view(3) 
end
```

![figure_10.png
](README_images/figure_10.png
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
Reconstruction error = 0.43707%
```

```matlab:Code

if ~isempty(outdof) && SSMDim<=2
idxPlot = [outdof]; % 3D Plot: eta_1, eta_2 and idxPlot coordinate
plotSSMandTrajectories(IMInfo, idxPlot, xDataTrunc(indTest,:), yDataTrunc(indTest,:))
view(-100,20); legend('off')
end
```

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
freqNorm = [1 3];

RDInfo = IMDynamicsMech(yDataTrunc(indTrain,:), ...
    'R_PolyOrd', 1,'N_PolyOrd', ROMOrder, 'style', 'normalform', ...
    'R_coeff',Ae,'rescale',1,'frequencies_norm',freqNorm,'MaxIter',5e3);
```

```text:Output
Estimation of the reduced dynamics...  Done. 
Estimation of the reduced dynamics in normal form...
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
     0           1       0.00317276                        0.0382
     1           2       0.00196118              1          0.026  
     2           3       0.00146781              1         0.0155  
     3           5       0.00143775        0.22715        0.00273  
     4           6       0.00142977              1        0.00151  
     5           7       0.00142439              1        0.00151  
     6           8       0.00141316              1        0.00362  
     7           9       0.00140279              1        0.00285  
     8          10       0.00139264              1        0.00198  
     9          11       0.00137469              1         0.0022  
    10          12       0.00136157              1        0.00245  
    11          13       0.00132991              1        0.00262  
    12          14       0.00130493              1        0.00707  
    13          16       0.00128548       0.425172        0.00423  
    14          17        0.0012607              1        0.00345  
    15          18       0.00121545              1        0.00311  
    16          19       0.00117923              1        0.00389  
    17          20        0.0010628              1        0.00811  
    ... 
   791         794      4.72275e-07              1       4.22e-08  
   792         795      4.72274e-07              1       4.25e-08  
   793         796      4.72274e-07              1       4.25e-08  
   794         797      4.72274e-07              1       4.25e-08  
   795         798      4.72274e-07              1       4.25e-08  
   796         799      4.72274e-07              1       4.25e-08  
   797         800      4.72274e-07              1       4.24e-08  
   798         801      4.72274e-07              1       4.24e-08  
   799         802      4.72272e-07              1       4.23e-08  
                       ...
```

![figure_11.png
](README_images/figure_11.png
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
Normalized mean trajectory error = 9.7356%
```

```matlab:Code

% We plot the true test set trajectory in the reduced coordinates and compare it to the prediction. 
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

![figure_12.png
](README_images/figure_12.png
)

We plot the model predictions in physical coordinates. The reduced model seems to do well on previously unseen data, provided that it is close to the manifold.

```matlab:Code
customFigure('subPlot',[3 2]); ipos = 0;
for iRow = 1:3
    for iCol = 1:2
        ipos = ipos+1;
        subplot(3,2,ipos)
        plot(xData{indTest(iRow),1},xData{indTest(iRow),2}(iCol,:),'k','Linewidth',1,'Color',colors(colFOM,:))
        plot(yRec{indTest(iRow),1},yRec{indTest(iRow),2}(iCol,:),':','Linewidth',2,'Color',colors(colSSML,:))
        xlabel('time')
        ylabel(['$u_' num2str(iCol) '$'],'interpreter','latex')
        xlim([xData{indTest(iRow),1}(1) xData{indTest(iRow),1}(end)])
    end
end
legend('Test trajectory','Prediction')
```

![figure_13.png
](README_images/figure_13.png
)

# Adding forcing to the ROM

```matlab:Code
% Outer directions: either consider all outer modes or a subset
numberOuterModes = (n-SSMDim)/2;
Vo = V(:,SSMDim+[1:2*numberOuterModes]); 
Wo = W(SSMDim+[1:2*numberOuterModes],:); Lo = full(Wo*A*Vo);
% Forcing vector
forcingVectors = [zeros(n,1); M\loadVector];
% Construct time periodic SSM model
[IMInfoF,RDInfoF] = forcedSSMROM(IMInfo,RDInfo,'nForcingFrequencies',1,...
         'forcingVectors',forcingVectors,'We',We,'Lo',Lo,'Vo',Vo, 'Wo',Wo);
```

# Generate Frequency Responses via SSMLearn \& SSMTool

We compute them also with SSMTool in order to compare the results (see the papers above for additional validations and comparisons).

```matlab:Code
mFreqs = [1 3];
resonantModes = [1 2 3 4];
omegaSpan = [0.975 1.06]*imag(lambda(1));
epsilon = [.4 .7];
outdofs = [outdof1 outdof2];
[FRCSSMLearn] = continuationFRCep(IMInfoF, RDInfoF, epsilon, omegaSpan,@(x) x(outdofs,:), mFreqs,resonantModes, 'SSMLearnFRC');
```

```text:Output
 Run='SSMLearnFRC.ep': Continue equilibria along primary branch.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          9.26e-01  1.13e+01    0.0    0.0    0.0
   1   1  6.03e-02  8.07e+00  8.86e-01  1.14e+01    0.0    0.0    0.0
   2   2  3.82e-02  5.03e+00  6.54e-01  1.14e+01    0.0    0.0    0.0
   3   3  2.50e-01  1.94e-02  1.01e-01  1.14e+01    0.0    0.0    0.0
   4   1  1.00e+00  5.31e-02  3.53e-02  1.14e+01    0.0    0.0    0.0
   5   1  1.00e+00  8.53e-03  8.98e-04  1.14e+01    0.0    0.0    0.0
   6   1  1.00e+00  5.96e-05  8.00e-07  1.14e+01    0.0    0.0    0.0
   7   1  1.00e+00  9.30e-09  7.08e-13  1.14e+01    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   1.1385e+01      1  EP      3.8649e+00   5.3285e-02   8.6382e-05   2.6224e-01   7.0511e+00   4.0000e-01
   10  00:00:00   1.4645e+01      2          3.8388e+00   2.5367e-02   8.6507e-06   1.0953e+00   9.5510e+00   4.0000e-01
   20  00:00:00   1.6039e+01      3          3.7801e+00   7.2415e-03   1.7086e-07   1.4417e+00   1.0592e+01   4.0000e-01
   22  00:00:00   1.6148e+01      4  EP      3.7589e+00   5.6816e-03   7.8221e-08   1.4702e+00   1.0678e+01   4.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   1.1385e+01      5  EP      3.8649e+00   5.3285e-02   8.6382e-05   2.6224e-01   7.0511e+00   4.0000e-01
   10  00:00:00   8.2116e+00      6          3.8713e+00   4.3268e-02   4.7432e-05  -6.6586e-01   4.2665e+00   4.0000e-01
   20  00:00:00   6.6155e+00      7          3.9045e+00   1.1103e-02   9.1040e-07  -1.3654e+00   2.1663e+00   4.0000e-01
   30  00:00:00   6.6035e+00      8          4.0770e+00   2.4853e-03   3.0385e-08  -1.5231e+00   1.6679e+00   4.0000e-01
   31  00:00:00   6.6111e+00      9  EP      4.0866e+00   2.3815e-03   3.0063e-08  -1.5250e+00   1.6575e+00   4.0000e-01
```

![figure_14.png
](README_images/figure_14.png
)

![figure_15.png
](README_images/figure_15.png
)

```text:Output
 Run='SSMLearnFRC.ep.ep': Continue equilibria along primary branch.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          8.96e-01  8.04e+00    0.0    0.0    0.0
   1   4  1.25e-01  2.62e-01  9.58e-01  8.03e+00    0.0    0.0    0.0
   2   1  4.71e-01  1.13e+00  8.72e-01  8.00e+00    0.0    0.0    0.0
   3   1  6.82e-02  6.20e+00  8.07e-01  7.99e+00    0.0    0.0    0.0
   4   3  3.96e-02  2.10e+00  7.13e-01  7.99e+00    0.0    0.0    0.0
   5   4  1.25e-01  1.05e-01  2.92e-01  7.99e+00    0.0    0.0    0.0
   6   2  5.00e-01  1.34e-01  6.93e-02  7.99e+00    0.0    0.0    0.0
   7   1  1.00e+00  9.32e-02  5.42e-03  7.99e+00    0.0    0.0    0.0
   8   1  1.00e+00  1.48e-02  4.13e-04  7.99e+00    0.0    0.0    0.0
   9   1  1.00e+00  2.32e-04  2.86e-07  7.99e+00    0.0    0.0    0.0
  10   1  1.00e+00  7.62e-08  5.59e-14  7.99e+00    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   7.9863e+00      1  EP      3.8887e+00   9.1965e-02  -4.7028e-04   3.0996e-01   4.0520e+00   7.0000e-01
   10  00:00:00   1.0777e+01      2          3.8430e+00   4.2756e-02  -4.1772e-05   1.1143e+00   6.4663e+00   7.0000e-01
   20  00:00:00   1.2025e+01      3          3.7781e+00   1.2283e-02  -8.2928e-07   1.4457e+00   7.4624e+00   7.0000e-01
   22  00:00:00   1.2111e+01      4  EP      3.7589e+00   9.9137e-03  -4.1548e-07   1.4705e+00   7.5370e+00   7.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   7.9863e+00      5  EP      3.8887e+00   9.1965e-02  -4.7028e-04   3.0996e-01   4.0520e+00   7.0000e-01
    6  00:00:00   6.7812e+00      6  SN      3.8962e+00   9.5701e-02  -5.4253e-04  -1.2517e-01   2.7463e+00   7.0000e-01
    6  00:00:00   6.7811e+00      7  FP      3.8962e+00   9.5702e-02  -5.4253e-04  -1.2520e-01   2.7462e+00   7.0000e-01
   10  00:00:00   5.8932e+00      8          3.8891e+00   7.8486e-02  -2.9565e-04  -6.1814e-01   1.2676e+00   7.0000e-01
   14  00:00:00   5.7194e+00      9  FP      3.8826e+00   5.0834e-02  -7.9744e-05  -1.0133e+00   8.2234e-02   7.0000e-01
   14  00:00:00   5.7194e+00     10  SN      3.8826e+00   5.0826e-02  -7.9703e-05  -1.0134e+00   8.1928e-02   7.0000e-01
   20  00:00:00   6.0210e+00     11          3.9017e+00   2.1200e-02  -6.2601e-06  -1.3467e+00  -9.1887e-01   7.0000e-01
   30  00:00:00   6.5103e+00     12          4.0629e+00   4.6452e-03  -1.7084e-07  -1.5200e+00  -1.4592e+00   7.0000e-01
   31  00:00:00   6.5534e+00     13  EP      4.0866e+00   4.1685e-03  -1.6119e-07  -1.5250e+00  -1.4841e+00   7.0000e-01
```

![figure_16.png
](README_images/figure_16.png
)

![figure_17.png
](README_images/figure_17.png
)

```matlab:Code
[FRCSSMTool] = SSMToolFRCFE(M,C,K,fnl,fext,outdofs,epsilon,resonantModes,3,omegaSpan,mFreqs, 'SSMToolFRC');
```

```text:Output
Due to high-dimensionality, we compute only the first 5 eigenvalues with the smallest magnitude. These would also be used to compute the spectral quotients
Assuming a proportional damping hypthesis with symmetric matrices
modal damping ratio for 1 mode is 2.593810e-03
modal damping ratio for 2 mode is 8.004681e-04
modal damping ratio for 3 mode is 3.837148e-04
modal damping ratio for 4 mode is 2.243713e-04
modal damping ratio for 5 mode is 1.470421e-04

 The first 10 nonzero eigenvalues are given as 
        -0.01 +     3.8553i
        -0.01 -     3.8553i
        -0.01 +     12.493i
        -0.01 -     12.493i
        -0.01 +     26.061i
        -0.01 -     26.061i
        -0.01 +     44.569i
        -0.01 -     44.569i
        -0.01 +     68.008i
        -0.01 -     68.008i

sigma_out = 1
sigma_in = 1
Manifold computation time at order 2 = 00:00:00
Estimated memory usage at order  2 = 6.78E-02 MB
Manifold computation time at order 3 = 00:00:00
Estimated memory usage at order  3 = 1.01E-01 MB

 Run='SSMToolFRCeps.ep': Continue equilibria with varied epsilon.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          5.15e-01  1.89e+01    0.0    0.0    0.0
   1   2  5.00e-01  5.68e-01  3.29e-02  1.89e+01    0.0    0.0    0.0
   2   1  1.00e+00  2.56e-01  1.51e-02  1.90e+01    0.0    0.0    0.0
   3   1  1.00e+00  2.15e-03  2.54e-04  1.90e+01    0.0    0.0    0.0
   4   1  1.00e+00  3.89e-06  6.92e-08  1.90e+01    0.0    0.0    0.0
   5   1  1.00e+00  1.15e-09  5.38e-15  1.90e+01    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE           eps         rho1         rho2          th1          th2           om
    0  00:00:00   1.8963e+01      1  EP      3.8985e-01   9.9110e+00   3.2178e-03   5.3819e+00   6.7103e+00   3.8553e+00
    4  00:00:00   1.8379e+01      2  EP      3.6000e-01   9.4643e+00   2.8038e-03   5.3376e+00   6.5773e+00   3.8553e+00

 STEP      TIME        ||U||  LABEL  TYPE           eps         rho1         rho2          th1          th2           om
    0  00:00:00   1.8963e+01      3  EP      3.8985e-01   9.9110e+00   3.2178e-03   5.3819e+00   6.7103e+00   3.8553e+00
    2  00:00:00   1.9152e+01      4  UZ      4.0000e-01   1.0055e+01   3.3593e-03   5.3960e+00   6.7525e+00   3.8553e+00
   10  00:00:00   2.2675e+01      5          6.4546e-01   1.2728e+01   6.7839e-03   5.6294e+00   7.4528e+00   3.8553e+00
   12  00:00:00   2.3278e+01      6  UZ      7.0000e-01   1.3186e+01   7.5373e-03   5.6632e+00   7.5543e+00   3.8553e+00
   13  00:00:00   2.3992e+01      7  EP      7.7000e-01   1.3729e+01   8.4993e-03   5.7009e+00   7.6673e+00   3.8553e+00

 Run='SSMToolFRCeps1.ep': Continue equilibria with varied omega at eps equal to 4.000000e-01.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          2.83e-17  1.95e+01    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   1.9532e+01      1  EP      3.8553e+00   1.0055e+01   3.3593e-03   5.3960e+00   6.7525e+00   4.0000e-01
   10  00:00:00   1.7792e+01      2          3.8448e+00   7.2697e+00   1.2322e-03   5.6881e+00   7.6293e+00   4.0000e-01
   20  00:00:00   1.6599e+01      3          3.8244e+00   3.8582e+00   1.7361e-04   5.9811e+00   8.5088e+00   4.0000e-01
   30  00:00:00   1.6556e+01      4          3.7616e+00   1.3734e+00   6.6149e-06   6.1771e+00   9.0982e+00   4.0000e-01
   31  00:00:00   1.6560e+01      5  EP      3.7589e+00   1.3364e+00   6.0550e-06   6.1800e+00   9.1069e+00   4.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   1.9532e+01      6  EP      3.8553e+00   1.0055e+01   3.3593e-03   5.3960e+00   6.7525e+00   4.0000e-01
   10  00:00:00   2.1335e+01      7          3.8656e+00   1.2593e+01   6.7974e-03   4.9537e+00   5.4253e+00   4.0000e-01
   20  00:00:00   2.1252e+01      8          3.8692e+00   1.2959e+01   7.4913e-03   4.6730e+00   4.5831e+00   4.0000e-01
   30  00:00:00   2.0112e+01      9          3.8708e+00   1.2382e+01   6.5764e-03   4.4105e+00   3.7954e+00   4.0000e-01
   40  00:00:00   1.5772e+01     10          3.8720e+00   9.3749e+00   2.8795e-03   3.9495e+00   2.4125e+00   4.0000e-01
   50  00:00:00   1.1508e+01     11          3.8774e+00   5.9951e+00   7.7027e-04   3.6221e+00   1.4300e+00   4.0000e-01
   60  00:00:01   8.1807e+00     12          3.9052e+00   2.5734e+00   6.7623e-05   3.3413e+00   5.8640e-01   4.0000e-01
   69  00:00:01   7.3817e+00     13  EP      4.0866e+00   5.6020e-01   2.3275e-06   3.1848e+00   8.6694e-02   4.0000e-01

 Run='SSMToolFRCeps2.ep': Continue equilibria with varied omega at eps equal to 7.000000e-01.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          1.12e-15  2.36e+01    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   2.3585e+01      1  EP      3.8553e+00   1.3186e+01   7.5373e-03   5.6632e+00   7.5543e+00   7.0000e-01
   10  00:00:00   2.0915e+01      2          3.8442e+00   1.0289e+01   3.4742e-03   5.8127e+00   8.0029e+00   7.0000e-01
   20  00:00:00   1.8390e+01      3          3.8272e+00   6.7932e+00   9.5340e-04   5.9792e+00   8.5030e+00   7.0000e-01
   30  00:00:00   1.6955e+01      4          3.7880e+00   3.2937e+00   9.7595e-05   6.1375e+00   8.9790e+00   7.0000e-01
   34  00:00:00   1.6790e+01      5  EP      3.7589e+00   2.3317e+00   3.2153e-05   6.1803e+00   9.1078e+00   7.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   2.3585e+01      6  EP      3.8553e+00   1.3186e+01   7.5373e-03   5.6632e+00   7.5543e+00   7.0000e-01
   10  00:00:00   2.6589e+01      7          3.8662e+00   1.6073e+01   1.4055e-02   5.4962e+00   7.0530e+00   7.0000e-01
   20  00:00:00   3.0451e+01      8          3.8801e+00   1.9519e+01   2.6136e-02   5.2478e+00   6.3072e+00   7.0000e-01
   30  00:00:00   3.3818e+01      9          3.8949e+00   2.2524e+01   4.1886e-02   4.8352e+00   5.0690e+00   7.0000e-01
   40  00:00:00   3.3864e+01     10          3.8970e+00   2.2687e+01   4.3097e-02   4.6849e+00   4.6179e+00   7.0000e-01
   48  00:00:00   3.3540e+01     11  SN      3.8974e+00   2.2522e+01   4.2244e-02   4.5886e+00   4.3291e+00   7.0000e-01
   48  00:00:00   3.3540e+01     12  FP      3.8974e+00   2.2522e+01   4.2244e-02   4.5886e+00   4.3291e+00   7.0000e-01
   50  00:00:00   3.3211e+01     13          3.8973e+00   2.2321e+01   4.1139e-02   4.5307e+00   4.1552e+00   7.0000e-01
   60  00:00:01   2.9099e+01     14          3.8919e+00   1.9518e+01   2.7234e-02   4.1769e+00   3.0939e+00   7.0000e-01
   70  00:00:01   2.4273e+01     15          3.8861e+00   1.6072e+01   1.5037e-02   3.9284e+00   2.3488e+00   7.0000e-01
   80  00:00:01   1.9538e+01     16          3.8831e+00   1.2592e+01   7.2127e-03   3.7298e+00   1.7528e+00   7.0000e-01
   83  00:00:01   1.8505e+01     17  SN      3.8830e+00   1.1816e+01   5.9660e-03   3.6892e+00   1.6311e+00   7.0000e-01
   83  00:00:01   1.8505e+01     18  FP      3.8830e+00   1.1816e+01   5.9659e-03   3.6892e+00   1.6311e+00   7.0000e-01
   90  00:00:01   1.4986e+01     19          3.8848e+00   9.1008e+00   2.7559e-03   3.5542e+00   1.2259e+00   7.0000e-01
  100  00:00:01   1.0851e+01     20          3.8971e+00   5.6032e+00   6.7537e-04   3.3911e+00   7.3594e-01   7.0000e-01
  110  00:00:01   7.8623e+00     21          3.9631e+00   2.1040e+00   4.7582e-05   3.2344e+00   2.6194e-01   7.0000e-01
  114  00:00:01   7.4910e+00     22  EP      4.0866e+00   9.8056e-01   1.2480e-05   3.1848e+00   8.6731e-02   7.0000e-01
Calculate FRC in physical domain at epsilon 4.000000e-01
the forcing frequency 3.7589e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
...
```

![figure_18.png
](README_images/figure_18.png
)

![figure_19.png
](README_images/figure_19.png
)

![figure_20.png
](README_images/figure_20.png
)

![figure_21.png
](README_images/figure_21.png
)

![figure_22.png
](README_images/figure_22.png
)

![figure_23.png
](README_images/figure_23.png
)

![figure_24.png
](README_images/figure_24.png
)

![figure_25.png
](README_images/figure_25.png
)

# **Plot results**

Phases.

```matlab:Code
% Phase plot outdof 1
plotfreqscale = omegaSpan(1);
customFigure; 
plotFRC(FRCSSMLearn, colors(colSSML,:), '- SSMLearn','y', 'Phase','freqscale',plotfreqscale)   
plotFRC(FRCSSMTool, colors(colSSMT,:), '- SSMTool','y', 'Phase','freqscale',plotfreqscale)
xlabel('forcing frequency [Hz]')
ylabel('phase mode 1')
xlim([1 1.06])
ylim([-pi 0 ])
set(gca,'YTick',[-pi -3/4*pi -pi/2 -pi/4 0])
set(gca,'YTickLabel',{'-\pi', '-3\pi/4', '-\pi/2', '-\pi/4', 0})
```

![figure_26.png
](README_images/figure_26.png
)

Amplitudes.

```matlab:Code
% Amplitude plot outdof 1 & 2
customFigure('subPlot',[ 2 1]); colors = colororder;
subplot(211)
plotFRC(FRCSSMLearn, colors(colSSML,:), '- SSMLearn','freqscale',plotfreqscale)
plotFRC(FRCSSMTool, colors(colSSMT,:), '- SSMTool','freqscale',plotfreqscale)
xlabel('$\Omega/\omega_{1}$ [-]','interpreter','latex')
ylabel('amp($u_{1}$) [-]','interpreter','latex')
legend('off')
xlim([1 1.06])
ylim([0 50])
subplot(212)
plotFRC(FRCSSMLearn, colors(colSSML,:), '- SSMLearn','freqscale',plotfreqscale,'outamp',2)
plotFRC(FRCSSMTool, colors(colSSMT,:), '- SSMTool','freqscale',plotfreqscale,'outamp',2)
xlabel('$\Omega/\omega_{1}$ [-]','interpreter','latex')
ylabel('amp($u_{2}$) [-]','interpreter','latex')
xlim([1 1.06])
```

![figure_27.png
](README_images/figure_27.png
)

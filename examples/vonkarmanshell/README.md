This is a preview of the livescript `vonkarmanshell.mlx`.

# Shallow-curved shell structure with geometric nonlinearities
  

See [1] for the details of this model, and [2] the description of this example.

[1] Jain, S., \& Tiso, P. (2018). Simulation-free hyper-reduction for geometrically nonlinear structural dynamics: a quadratic manifold lifting approach. *Journal of Computational and Nonlinear Dynamics*, *13*(7), 071003. [https://doi.org/10.1115/1.4040021](https://doi.org/10.1115/1.4040021)

[2] Cenedese, M., Marconi, J., Haller, G., \& Jain, S. (2023). Data-assisted non-intrusive model reduction for forced nonlinear finite elements models. Preprint: [arXiv: 2311.17865](https://arxiv.org/abs/2311.17865) 

The finite element code taken from the following package:

Jain, S., Marconi, J., Tiso P. (2020). YetAnotherFEcode (Version v1.1). Zenodo. [http://doi.org/10.5281/zenodo.4011282](http://doi.org/10.5281/zenodo.4011282)

See the `README` of the main repository to retrieve simulations data.

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

![figure_3.png
](README_images/figure_3.png
)

![figure_4.png
](README_images/figure_4.png
)

```text:Output
Using Rayleigh damping
a = 2x1    
   8.6312e-06
      0.40215

Assembling external force vector
Getting nonlinearity coefficients
Loaded tensors from storage
Total time spent on model assembly = 00:00:17
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

```matlab:Code
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
% set(DS.Options,'Emax',5,'Nmax',10,'notation','tensor')
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

We initialize the base properties of the SSM, i.e., its linear part, which we know from the linear dynamics of the model. In this case, we target the slow two-dimensional SSM of the system. 

```matlab:Code
masterModes = [1 2]; % Modal displacements and modal velocities
Ve = V(:,masterModes); % Mode shape
We = W(masterModes,:); % Projection to mode shape
Ae = full(We*A*Ve) % Reduced, linearized dynamics
```

```text:Output
Ae = 2x2    
            0            1
       -21743     -0.58982

```

```matlab:Code
SSMDim = length(masterModes);

displacementVector = zeros(1,n); displacementVector(outdof) = 1;
loadVector = f_0;  %  could also be set as modal ones
```

# Load static analysis and plot

We characterize the linear and nonlinear regimes via a static analysis, which serves to pick appropriate initial conditions for the trajectories we need to learn the SSM from data.

```matlab:Code
load("StaticDisplacements.mat")
nsteps = size(uNonlinear,2);
loadCoefficients = linspace(0,100,nsteps);

% Define observed deformation 
uNonlinearOut = uNonlinear(outdof,:);
uLinearOut = uLinear(outdof,:);
relativeDispDiff = 100*abs(uNonlinearOut-uLinearOut)./abs(uNonlinearOut);

customFigure;
plot(loadCoefficients,uNonlinearOut,'r--','DisplayName','Linear')
plot(loadCoefficients,uLinearOut,'-k','DisplayName','Nonlinear')
legend('location','best')
ylabel('$$u_{\mathrm{out}}$$ [m]','Interpreter','latex'); 
xlabel('Normalized load'); 
title('Static loading analysis');
axis tight
```

![figure_5.png
](README_images/figure_5.png
)

```matlab:Code

Z = [uNonlinear; zeros(n,nsteps)];
```

Pick up two initial trajectories that has high expected nonlinear content

```matlab:Code
indIC = [nsteps, nsteps-1];
IC = Z(:,indIC);
indTrain = 1;
indTest = 2;
```

Define the linear regime at 0.05 % relative displacement

```matlab:Code
linearDisplacementReference = uNonlinearOut(sum(relativeDispDiff<(0.05))+1);
nonlinearDisplacementReference = max(abs(uNonlinearOut(indIC)));
desiredAmplitudeDecay = nonlinearDisplacementReference/linearDisplacementReference;
```

# **Generate decaying trajectories via time integration**

We define observables and timescales. The computation of integration time is estimated from the linear decay that gets from the defined nonlinear amplitude to linear regime. We set the sampling time to capture approximately a fixed number points per period on the faster time scale. Then, we integrate using the initial conditions we obtained from the static analysis. Here, we use a pre-computed data set to avoid excessive computations.

```matlab:Code
observable = @(x) x; % Observe the full phase space
slowTimeScale = 2*pi/abs(lambda(1));
fastTimeScale = 2*pi/abs(lambda(round(SSMDim/2)));
% The computation of integration time is estimated from the linear decay 
% that gets from the nonlinear amplitude to linear regime.
numberPeriodsSlow = floor(log(desiredAmplitudeDecay)/...
    (2*pi*(-real(lambda(1))/abs(lambda(1)))))
```

```text:Output
numberPeriodsSlow = 
   201

```

```matlab:Code
endTime = numberPeriodsSlow*slowTimeScale;
% Set the sampling time to capture approximately 100 points per period on 
% the faster time scale
numberPeriodsFast = floor(endTime/fastTimeScale);
numberPointsPerPeriod = 100;
nSamp = numberPeriodsFast*numberPointsPerPeriod+1;
dt = endTime/(nSamp-1);
% xData = integrateTrajectoriesGalphaDirect(Model, DS, endTime, IC, nSamp, observable);
% % xData = integrateTrajectoriesNewmarkDirect(Model, DS, endTime, IC, nSamp, observable);
% DataInfo = struct('nElements', Model.Mesh.nElements, 'loadShape', loadShape);
% save('dataVKDecay2DGalpha.mat', 'DataInfo', 'xData', 'dt', 'endTime', 'nSamp')
load("../../data/vonkarmanshell/dataVKDecay2DGalphaStatic.mat")
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
showSpectrogram(xData(indTrain,:), outdof);
ylim([0,abs(lambda(1))/2/pi*5])
```

![figure_6.png
](README_images/figure_6.png
)

We plot theobservables of interest over time for closer inspection. 

```matlab:Code
customFigure();
plot(xData{1,1}, xData{1,2}(outdof,:), xData{2,1}, xData{2,2}(outdof,:), ':');
xlabel('$t \, [$s$]$','Interpreter','latex'); ylabel('$u \, [$m$]$','Interpreter','latex'); 
legend({'Trajectory 1', 'Trajectory 2'})
title('Generated data')
```

![figure_7.png
](README_images/figure_7.png
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

where the function ${{\phi }}_{m,2:M} (y)$ computes a vector of all $m_M$ monomials from orders 2 up to $M$ of an $m$-vector $y$. From SSM theory, the tangent space of the manifold is $V_e$. The coefficients $H$are obtained via least squares regression.

```matlab:Code
SSMOrder = 7;

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
```

![figure_8.png
](README_images/figure_8.png
)

```matlab:Code
if SSMDim>2
   view(3) 
end

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
Reconstruction error = 1.0349%
```

# Plot SSM

We plot the SSM geometry in different coordinates.

```matlab:Code
idxPlot = [outdof]; % 3D Plot: eta_1, eta_2 and idxPlot coordinate
plotSSMandTrajectories(IMInfo, idxPlot, xDataTrunc(indTest,:), yDataTrunc(indTest,:))
view(-100,20); legend('off')
```

![figure_9.png
](README_images/figure_9.png
)

```matlab:Code
gPlot = @(x) W(5,:)*x; % 3D Plot: eta_1, eta_2 and gPlot(x) 
plotSSMandTrajectories(IMInfo, gPlot, xDataTrunc(indTest,:), yDataTrunc(indTest,:))
view(-100,20); legend('off')
```

![figure_10.png
](README_images/figure_10.png
)

```matlab:Code
gPlot = @(x) W([3],:)*x; % 3D Plot with the three values of gPlot(x) 
plotSSMandTrajectories(IMInfo, gPlot, xDataTrunc(indTest,:), yDataTrunc(indTest,:))
view(-100,20); legend('off')
xlabel('$u_1$','interpreter','latex')
ylabel('$\dot{u}_1$','interpreter','latex')
zlabel('$u_2$','interpreter','latex')
```

![figure_11.png
](README_images/figure_11.png
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
ROMOrder = 7;

RDInfo = IMDynamicsMech(yDataTrunc(indTrain,:), ...
    'R_PolyOrd', 1,'N_PolyOrd', ROMOrder, 'style', 'normalform', ...
    'R_coeff',Ae,'rescale',1,'MaxIter',5e3);
```

```text:Output
Estimation of the reduced dynamics...  Done. 
Estimation of the reduced dynamics in normal form...
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
     0           1          2.77856                            24
     1           3          1.87811     0.00416699           21.5  
     2           4          1.32201              1           12.5  
     3           6          1.15664       0.446756           4.69  
     4           7           1.0076              1           1.87  
     5           8          0.97563              1          0.358  
     6           9         0.974841              1          0.356  
     7          10         0.969736              1          0.518  
     8          11         0.966971              1           0.54  
     9          12         0.961994              1          0.614  
    10          13         0.956859              1           0.61  
    11          14         0.951337              1           0.66  
    12          15         0.947452              1          0.653  
    ... 
   490         493         0.275881              1       0.000648  
   491         494         0.275881              1        0.00054  
   492         495         0.275881              1       0.000234  
   493         496         0.275881              1       3.62e-05  
   494         497         0.275881              1       9.98e-07  
   495         498         0.275881              1          1e-06  
   496         499         0.275881              1          1e-06  
   497         500         0.275881              1       9.99e-07  

Local minimum possible.

fminunc stopped because the size of the current step is less than
the value of the step size tolerance.

<stopping criteria details>
Plotting figure with the polar normal form equations ... Done. 
```

![figure_12.png
](README_images/figure_12.png
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
Normalized mean trajectory error = 6.3316%
```

```matlab:Code

% We plot the true test set trajectory in the reduced coordinates and compare it to the prediction. 
plotReducedCoordinates(yDataTrunc(indTest,:), etaRec(indTest,:))
legend({'Test set (truncated)', 'Prediction'})
```

![figure_13.png
](README_images/figure_13.png
)

```matlab:Code
if size(Ae,1)==2
    % Plot SSM with trajectories in the normal form reduced coordinates
    plotSSMandTrajectories(IMInfo, outdof, xDataTrunc(indTest,:), ...
        zRec(indTest,:), 'NFT', RDInfo.transformation.map)
    view(-100,20); legend('off')
else
    view(3)
end
```

![figure_14.png
](README_images/figure_14.png
)

We plot the model predictions in physical coordinates. The reduced model seems to do well on previously unseen data, provided that it is close to the manifold.

```matlab:Code
plotTrajectories(xData(indTest,:), yRec(indTest,:), 'm','PlotCoordinate', outdof(1), 'DisplayName', {'Test set', 'Prediction'})
xlabel('$t \, [$s$]$','Interpreter','latex')
ylabel('$q_A \, [$m$]$','Interpreter','latex')
```

![figure_15.png
](README_images/figure_15.png
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
epsilon = [0.1 0.2];
mFreqs = [1];
resonantModes = [1 2];
omegaSpan = [0.92 1.07]*imag(lambda(1));
% SSMLearn
[FRCSSMLearn] = continuationFRCep(IMInfoF, RDInfoF, epsilon, omegaSpan,@(x) x(outdof,:), mFreqs,resonantModes, 'SSMLearnFRC');
```

```text:Output
 Run='SSMLearnFRC.ep': Continue equilibria along primary branch.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          3.70e-04  2.09e+02    0.0    0.0    0.0
   1   1  1.00e+00  8.40e-04  5.48e-08  2.09e+02    0.1    0.1    0.0
   2   1  1.00e+00  1.12e-07  2.88e-15  2.09e+02    0.1    0.2    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1          th1          eps
    0  00:00:00   2.0867e+02      1  EP      1.4745e+02   1.1416e-01   5.3150e+00   1.0000e-01
   10  00:00:00   2.0702e+02      2          1.4625e+02   1.8703e-01   6.3141e+00   1.0000e-01
   14  00:00:00   2.0700e+02      3  FP      1.4623e+02   1.8456e-01   6.4671e+00   1.0000e-01
   14  00:00:00   2.0700e+02      4  SN      1.4623e+02   1.8456e-01   6.4671e+00   1.0000e-01
   20  00:00:00   2.0718e+02      5          1.4634e+02   1.6262e-01   6.8494e+00   1.0000e-01
   27  00:00:01   2.0752e+02      6  SN      1.4656e+02   1.0226e-01   7.3309e+00   1.0000e-01
   27  00:00:01   2.0752e+02      7  FP      1.4656e+02   1.0226e-01   7.3309e+00   1.0000e-01
   30  00:00:01   2.0751e+02      8          1.4655e+02   9.1414e-02   7.3955e+00   1.0000e-01
   40  00:00:01   2.0733e+02      9          1.4641e+02   6.6737e-02   7.5309e+00   1.0000e-01
   50  00:00:01   2.0458e+02     10          1.4445e+02   2.0906e-02   7.7574e+00   1.0000e-01
   60  00:00:01   1.9959e+02     11          1.4092e+02   9.5960e-03   7.8109e+00   1.0000e-01
   70  00:00:01   1.9460e+02     12          1.3738e+02   6.2298e-03   7.8267e+00   1.0000e-01
   75  00:00:01   1.9217e+02     13  EP      1.3566e+02   5.3207e-03   7.8310e+00   1.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1          th1          eps
    0  00:00:01   2.0867e+02     14  EP      1.4745e+02   1.1416e-01   5.3150e+00   1.0000e-01
   10  00:00:01   2.1263e+02     15          1.5028e+02   2.2000e-02   4.8181e+00   1.0000e-01
   20  00:00:01   2.1763e+02     16          1.5381e+02   9.8609e-03   4.7607e+00   1.0000e-01
   30  00:00:01   2.2262e+02     17          1.5735e+02   6.3429e-03   4.7442e+00   1.0000e-01
   32  00:00:02   2.2323e+02     18  EP      1.5778e+02   6.0789e-03   4.7430e+00   1.0000e-01
```

![figure_16.png
](README_images/figure_16.png
)

![figure_17.png
](README_images/figure_17.png
)

```text:Output
 Run='SSMLearnFRC.ep.ep': Continue equilibria along primary branch.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          2.97e-03  2.09e+02    0.0    0.0    0.0
   1   1  1.00e+00  1.89e-03  1.16e-06  2.09e+02    0.0    0.0    0.0
   2   1  1.00e+00  4.36e-07  6.70e-14  2.09e+02    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1          th1          eps
    0  00:00:00   2.0866e+02      1  EP      1.4745e+02   1.4968e-01   5.1079e+00   2.0000e-01
   10  00:00:00   2.0472e+02      2          1.4464e+02   3.0103e-01   5.8985e+00   2.0000e-01
   19  00:00:00   2.0404e+02      3  SN      1.4414e+02   3.1763e-01   6.3740e+00   2.0000e-01
   19  00:00:00   2.0404e+02      4  FP      1.4414e+02   3.1763e-01   6.3740e+00   2.0000e-01
   20  00:00:00   2.0404e+02      5          1.4414e+02   3.1704e-01   6.3988e+00   2.0000e-01
   30  00:00:00   2.0427e+02      6          1.4429e+02   3.0375e-01   6.6406e+00   2.0000e-01
   40  00:00:00   2.0665e+02      7          1.4593e+02   1.4312e-01   7.4832e+00   2.0000e-01
   47  00:00:00   2.0669e+02      8  FP      1.4596e+02   1.2482e-01   7.5387e+00   2.0000e-01
   47  00:00:00   2.0669e+02      9  SN      1.4596e+02   1.2482e-01   7.5387e+00   2.0000e-01
   50  00:00:00   2.0668e+02     10          1.4595e+02   1.1640e-01   7.5631e+00   2.0000e-01
   60  00:00:00   2.0653e+02     11          1.4584e+02   9.3203e-02   7.6271e+00   2.0000e-01
   70  00:00:00   2.0315e+02     12          1.4344e+02   3.1435e-02   7.7818e+00   2.0000e-01
   80  00:00:01   1.9816e+02     13          1.3990e+02   1.6632e-02   7.8169e+00   2.0000e-01
   90  00:00:01   1.9317e+02     14          1.3637e+02   1.1324e-02   7.8294e+00   2.0000e-01
   93  00:00:01   1.9217e+02     15  EP      1.3566e+02   1.0644e-02   7.8310e+00   2.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1          th1          eps
    0  00:00:01   2.0866e+02     16  EP      1.4745e+02   1.4968e-01   5.1079e+00   2.0000e-01
   10  00:00:01   2.1276e+02     17          1.5037e+02   4.1978e-02   4.8139e+00   2.0000e-01
   20  00:00:01   2.1776e+02     18          1.5390e+02   1.9408e-02   4.7601e+00   2.0000e-01
   30  00:00:01   2.2275e+02     19          1.5744e+02   1.2562e-02   4.7439e+00   2.0000e-01
   31  00:00:01   2.2323e+02     20  EP      1.5778e+02   1.2153e-02   4.7430e+00   2.0000e-01
```

![figure_18.png
](README_images/figure_18.png
)

![figure_19.png
](README_images/figure_19.png
)

```matlab:Code
% SSMTool
FRCSSMTool = SSMToolFRCFE(M,C,K,fnl,f_0,outdof,epsilon,resonantModes,ROMOrder,omegaSpan,mFreqs,'SSMToolFRC');
```

```text:Output
Due to high-dimensionality, we compute only the first 5 eigenvalues with the smallest magnitude. These would also be used to compute the spectral quotients
Assuming a proportional damping hypthesis with symmetric matrices
modal damping ratio for 1 mode is 2.000000e-03
modal damping ratio for 2 mode is 2.000000e-03
modal damping ratio for 3 mode is 2.269789e-03
modal damping ratio for 4 mode is 2.721500e-03
modal damping ratio for 5 mode is 2.909530e-03
the left eigenvectors may be incorrect in case of asymmetry of matrices

 The first 10 nonzero eigenvalues are given as 
     -0.29491 +     147.45i
     -0.29491 -     147.45i
     -0.63196 +     315.98i
     -0.63196 -     315.98i
     -0.93785 +     413.19i
     -0.93785 -     413.19i
      -1.4836 +     545.15i
      -1.4836 -     545.15i
      -1.7341 +     596.01i
      -1.7341 -     596.01i

(near) outer resonance detected for the following combination of master eigenvalues
     2     0
     3     1
     0     2
     1     3
     3     0
     4     1
     0     3
     1     4
     4     0
     0     4
     4     0
     0     4

These are in resonance with the follwing eigenvalues of the slave subspace
     -0.63196 +     315.98i
     -0.63196 +     315.98i
     -0.63196 -     315.98i
     -0.63196 -     315.98i
     -0.93785 +     413.19i
     -0.93785 +     413.19i
     -0.93785 -     413.19i
     -0.93785 -     413.19i
      -1.4836 +     545.15i
      -1.4836 -     545.15i
      -1.7341 +     596.01i
      -1.7341 -     596.01i

sigma_out = 5
(near) inner resonance detected for the following combination of master eigenvalues
     2     1
     3     2
     1     2
     2     3

These are in resonance with the follwing eigenvalues of the master subspace
     -0.29491 +     147.45i
     -0.29491 +     147.45i
     -0.29491 -     147.45i
     -0.29491 -     147.45i

sigma_in = 5
Due to (near) outer resonance, the exisitence of the manifold is questionable and the underlying computation may suffer.
Attempting manifold computation
Manifold computation time at order 2 = 00:00:00
Estimated memory usage at order  2 = 1.30E+01 MB
Manifold computation time at order 3 = 00:00:01
Estimated memory usage at order  3 = 1.41E+01 MB
Manifold computation time at order 4 = 00:00:01
Estimated memory usage at order  4 = 1.52E+01 MB
Manifold computation time at order 5 = 00:00:01
Estimated memory usage at order  5 = 1.70E+01 MB
Manifold computation time at order 6 = 00:00:03
Estimated memory usage at order  6 = 1.91E+01 MB
Manifold computation time at order 7 = 00:00:04
Estimated memory usage at order  7 = 2.23E+01 MB

 Run='SSMToolFRCeps.ep': Continue equilibria with varied epsilon.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          2.55e-03  1.48e+02    0.0    0.0    0.0
   1   1  1.00e+00  9.10e-04  3.23e-07  1.48e+02    0.0    0.0    0.0
   2   1  1.00e+00  3.36e-07  3.10e-14  1.48e+02    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE           eps         rho1          th1           om
    0  00:00:00   1.4755e+02      1  EP      1.0000e-01   3.9582e-02   3.7344e+00   1.4745e+02
    1  00:00:00   1.4755e+02      2  UZ      1.0000e-01   3.9582e-02   3.7344e+00   1.4745e+02
    1  00:00:00   1.4755e+02      3  EP      9.0000e-02   3.7838e-02   3.7745e+00   1.4745e+02

 STEP      TIME        ||U||  LABEL  TYPE           eps         rho1          th1           om
    0  00:00:00   1.4755e+02      4  EP      1.0000e-01   3.9582e-02   3.7344e+00   1.4745e+02
    5  00:00:00   1.4754e+02      5  UZ      2.0000e-01   5.1962e-02   3.5277e+00   1.4745e+02
    6  00:00:00   1.4754e+02      6  EP      2.2000e-01   5.3830e-02   3.5059e+00   1.4745e+02

 Run='SSMToolFRCeps1.ep': Continue equilibria with varied omega at eps equal to 1.000000e-01.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          3.10e-14  2.09e+02    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1          th1          eps
    0  00:00:00   2.0860e+02      1  EP      1.4745e+02   3.9582e-02   3.7344e+00   1.0000e-01
   10  00:00:00   2.0688e+02      2          1.4621e+02   6.6380e-02   4.7322e+00   1.0000e-01
   14  00:00:00   2.0685e+02      3  FP      1.4619e+02   6.5602e-02   4.8807e+00   1.0000e-01
   14  00:00:00   2.0685e+02      4  SN      1.4619e+02   6.5602e-02   4.8807e+00   1.0000e-01
   20  00:00:00   2.0701e+02      5          1.4629e+02   5.8890e-02   5.2305e+00   1.0000e-01
   27  00:00:00   2.0742e+02      6  FP      1.4656e+02   3.5347e-02   5.7652e+00   1.0000e-01
   27  00:00:00   2.0742e+02      7  SN      1.4656e+02   3.5347e-02   5.7652e+00   1.0000e-01
   30  00:00:00   2.0742e+02      8          1.4655e+02   3.2416e-02   5.8143e+00   1.0000e-01
   40  00:00:00   2.0727e+02      9          1.4644e+02   2.4221e-02   5.9427e+00   1.0000e-01
   50  00:00:00   2.0497e+02     10          1.4480e+02   8.1970e-03   6.1714e+00   1.0000e-01
   60  00:00:00   1.9998e+02     11          1.4127e+02   3.5063e-03   6.2355e+00   1.0000e-01
   70  00:00:00   1.9499e+02     12          1.3773e+02   2.2316e-03   6.2528e+00   1.0000e-01
   76  00:00:00   1.9205e+02     13  EP      1.3566e+02   1.8391e-03   6.2582e+00   1.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1          th1          eps
    0  00:00:01   2.0860e+02     14  EP      1.4745e+02   3.9582e-02   3.7344e+00   1.0000e-01
   10  00:00:01   2.1260e+02     15          1.5030e+02   7.5474e-03   3.2445e+00   1.0000e-01
   20  00:00:01   2.1760e+02     16          1.5383e+02   3.3968e-03   3.1878e+00   1.0000e-01
   30  00:00:01   2.2260e+02     17          1.5737e+02   2.1876e-03   3.1713e+00   1.0000e-01
   32  00:00:01   2.2318e+02     18  EP      1.5778e+02   2.1012e-03   3.1702e+00   1.0000e-01

 Run='SSMToolFRCeps2.ep': Continue equilibria with varied omega at eps equal to 2.000000e-01.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          6.14e-14  2.09e+02    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1          th1          eps
    0  00:00:00   2.0859e+02      1  EP      1.4745e+02   5.1962e-02   3.5277e+00   2.0000e-01
   10  00:00:00   2.0460e+02      2          1.4461e+02   1.0533e-01   4.2499e+00   2.0000e-01
   20  00:00:00   2.0354e+02      3          1.4385e+02   1.1465e-01   4.7658e+00   2.0000e-01
   21  00:00:00   2.0354e+02      4  SN      1.4385e+02   1.1459e-01   4.7784e+00   2.0000e-01
   21  00:00:00   2.0354e+02      5  FP      1.4385e+02   1.1459e-01   4.7784e+00   2.0000e-01
   30  00:00:00   2.0364e+02      6          1.4391e+02   1.1264e-01   4.9310e+00   2.0000e-01
   40  00:00:00   2.0558e+02      7          1.4526e+02   7.9939e-02   5.6050e+00   2.0000e-01
   50  00:00:00   2.0659e+02      8          1.4596e+02   4.3746e-02   5.9666e+00   2.0000e-01
   51  00:00:00   2.0659e+02      9  FP      1.4596e+02   4.3317e-02   5.9701e+00   2.0000e-01
   51  00:00:00   2.0659e+02     10  SN      1.4596e+02   4.3317e-02   5.9701e+00   2.0000e-01
   60  00:00:00   2.0653e+02     11          1.4592e+02   3.6361e-02   6.0252e+00   2.0000e-01
   70  00:00:00   2.0486e+02     12          1.4473e+02   1.6268e-02   6.1717e+00   2.0000e-01
   80  00:00:00   1.9987e+02     13          1.4119e+02   6.9364e-03   6.2360e+00   2.0000e-01
   90  00:00:01   1.9488e+02     14          1.3766e+02   4.4293e-03   6.2531e+00   2.0000e-01
   96  00:00:01   1.9205e+02     15  EP      1.3566e+02   3.6791e-03   6.2582e+00   2.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1          th1          eps
    0  00:00:01   2.0859e+02     16  EP      1.4745e+02   5.1962e-02   3.5277e+00   2.0000e-01
   10  00:00:01   2.1271e+02     17          1.5037e+02   1.4494e-02   3.2408e+00   2.0000e-01
   20  00:00:01   2.1770e+02     18          1.5391e+02   6.7051e-03   3.1872e+00   2.0000e-01
   30  00:00:01   2.2270e+02     19          1.5744e+02   4.3408e-03   3.1711e+00   2.0000e-01
   31  00:00:01   2.2318e+02     20  EP      1.5778e+02   4.2007e-03   3.1702e+00   2.0000e-01
Calculate FRC in physical domain at epsilon 1.000000e-01
the forcing frequency 1.3566e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
...
the forcing frequency 1.5778e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
Calculate FRC in physical domain at epsilon 2.000000e-01
the forcing frequency 1.3566e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
...
the forcing frequency 1.5778e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
```

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

# **Plot results**

```matlab:Code
customFigure('subPlot',[2 1]); 
subplot(211)
plotFRC(FRCSSMLearn, colors(colSSML,:), '- SSMLearn','freqscale',2*pi)
plotFRC(FRCSSMTool, colors(colSSMT,:), '- Tool','freqscale',2*pi)
xlabel('$$\Omega$$ [Hz]','interpreter','latex')
ylabel('$$\mathrm{amp}(q_{A})$$ [m]','interpreter','latex')
legend('off')
%xlim([ceil(omegaSpan(1)*100)/100 floor(omegaSpan(2)*100)/100]/(2*pi))
xlim([22 25])
ylim([0 0.007])
subplot(212)
plotFRC(FRCSSMLearn, colors(colSSML,:), '- SSMLearn','y','Phase','freqscale',2*pi)
plotFRC(FRCSSMTool, colors(colSSMT,:), '- SSMTool','y','Phase','freqscale',2*pi)
xlabel('$$\Omega$$ [Hz]','interpreter','latex')
ylabel('$$\mathrm{phase}(q_{A})$$ [rad]','interpreter','latex')
ylim([-pi 0])
set(gca,'YTick',[-pi -3/4*pi -pi/2 -pi/4 0])
set(gca,'YTickLabel',{'-\pi', '-3\pi/4', '-\pi/2', '-\pi/4', '0'})
%xlim([ceil(omegaSpan(1)*100)/100 floor(omegaSpan(2)*100)/100]/(2*pi))
xlim([22 25])
```

![figure_24.png
](README_images/figure_24.png
)

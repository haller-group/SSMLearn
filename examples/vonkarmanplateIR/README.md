This is a preview of the livescript `vonkarmanplateIR.mlx`.

# A von Karman plate with 1:1 internal resonance

In this example, we consider a simply-supported von Karman square plate subject to harmonic excitation. Due to the geometric symmetry, the natural frequenices of the second and the third modes are the same. In other words, the system has 1:1 internal resonance between the two modes. We extract the forced response curve using SSM reduction. The particularity of this example is that the SSM is not a slow one, but an intermediate SSM and has dimension four.

 See [1] for the details of this model, and [2] for the description of this example.

[1] Li, M., Jain, S., \& Haller, G. (2021). Nonlinear analysis of forced mechanical systems with internal resonance using spectral submanifolds-Part I: Periodic response and forced response curve. *Nonlinear Dynamics* 110, 1005-1043. [DOI: 10.1007/s11071-022-07714-x](https://doi.org/10.1007/s11071-022-07714-x)

[2] Cenedese, M., Marconi, J., Haller, G., \& Jain, S. (2023). Data-assisted non-intrusive model reduction for forced nonlinear finite elements models. Preprint: [arXiv: 2311.17865](https://arxiv.org/abs/2311.17865) 

The finite element code taken from the following package:

Jain, S., Marconi, J., Tiso P. (2020). YetAnotherFEcode (Version v1.1). Zenodo. [http://doi.org/10.5281/zenodo.4011282](http://doi.org/10.5281/zenodo.4011282) 

See the `README` of the main repository to retrieve simulations data.

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

  
> $x=\left\lbrack \begin{array}{c}
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
l = 1; % length of domain [m]
b = 1;  % breadth of domain [m]
t = 1e-2; % thickness of plate [m]
w = 0.0; % curvature parameter (height of the midpoint relative to ends) [m]
% material properties
E     = 70e9;  % 70e9 % 200e9 % Young's modulus [Pa]
rho   = 2700; % 2700 % 7850 % density [kg/m^3]
nu    = 0.33;    % Poisson's ratio 
kappa = 1e5; % material damping modulus 1e8

% Mesh
nElements = 10;
nl = nElements;
nb = nElements;
bc = 'SSSS';

```

# Generate model

```matlab:Code
[M,C,K,fnl,~,outdof, Model] = buildModel(l,b,t,w,E,rho,nu,kappa,bc,nl,nb);
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
The circular natural frequencies of the first five modes:
omega = 5x1    
       306.69
       763.59
       767.75
       1218.8
       1531.1

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
Getting nonlinearity coefficients
Loaded tensors from storage
Assembling external force vector
```

```matlab:Code
n = length(M); % number of degrees of freedom
disp(['Number of degrees of freedom = ' num2str(n)])
```

```text:Output
Number of degrees of freedom = 606
```

```matlab:Code
disp(['Phase space dimensionality = ' num2str(2*n)])
```

```text:Output
Phase space dimensionality = 1212
```

```matlab:Code
f_0 = zeros(n,1);
f_0(outdof(1)) = 100;
epsilon = 0.5;

```

Preliminaries: import the linear part without recomputing it.

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

We initialize the base properties of the SSM, i.e., its linear part, which we know from the linear dynamics of the model. In this case, we target an intermediate four-dimensional SSM of the system, which features an internal $1:1$ resonance.

```matlab:Code
masterModes = [3 5 4 6]; % Modal displacements and modal velocities
SSMDim = length(masterModes);
Ve = V(:,masterModes); % Mode shape
We = W(masterModes,:); % Projection to mode shape
Ae = full(We*A*Ve) % Reduced, linearized dynamics
```

```text:Output
Ae = 4x4    
            0            0            1   2.0339e-13
            0            0   2.0345e-13            1
  -5.8308e+05  -1.2276e-07      -3.3323  -6.8996e-13
  -1.2034e-07  -5.8944e+05  -6.8166e-13      -3.3578

```

Load and displacement vector

```matlab:Code
displacementVector = zeros(1,n); displacementVector(outdof) = 1;
loadVector = f_0;  %  could also be set as modal ones
```

# Compare linear and nonlinear response via modal displacement

We characterize the linear and nonlinear regimes via a static modal analysis, which serves to pick appropriate initial conditions for the trajectories we need to learn the SSM from data.

Displacement along the first mode

```matlab:Code
iMode = 2; scalingfactor1 = 1e0; nsteps = 50; outdof1 = outdof(1);
[phi1, relativeDiffForceNorm] = modal_analysis(Model,scalingfactor1,nsteps,outdof1,true,iMode);
```

```text:Output
Solving undamped eigenvalue problem
Eigenfrequency
       763.59
```

![figure_5.png
](README_images/figure_5.png
)

![figure_6.png
](README_images/figure_6.png
)

![figure_7.png
](README_images/figure_7.png
)

![figure_8.png
](README_images/figure_8.png
)

```text:Output
Displacement at output DOF: 0.0037652
```

Pick up two initial trajectories that has high expected nonlinear content

```matlab:Code
indIC1 = [nsteps, nsteps-1];
IC1 = [phi1*(scalingfactor1*indIC1/nsteps);zeros(n,length(indIC1))];
```

Displacement along the second mode

```matlab:Code
iMode = 3; scalingfactor2 = 1e-1; nsteps = 50; outdof2 = outdof(1);
[phi2, relativeDiffForceNorm2] = modal_analysis(Model,scalingfactor2,nsteps,outdof2,true,iMode);
```

```text:Output
Solving undamped eigenvalue problem
Eigenfrequency
       767.75
```

![figure_9.png
](README_images/figure_9.png
)

![figure_10.png
](README_images/figure_10.png
)

![figure_11.png
](README_images/figure_11.png
)

![figure_12.png
](README_images/figure_12.png
)

```text:Output
Displacement at output DOF: 0.0023927
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
observable = @(x) x; % Observe the full phase space
slowTimeScale = 2*pi/abs(lambda(2));
fastTimeScale = 2*pi/abs(lambda(3));
% The computation of integration time is estimated from the linear decay that 
% gets from the nonlinear amplitude to linear regime.
newSimulation = false;
if newSimulation
    numberPeriodsSlow = floor(log(desiredAmplitudeDecay)/...
        (2*pi*(-real(lambda(1))/abs(lambda(1)))))
    endTime = numberPeriodsSlow*slowTimeScale;
    % Set the sampling time to capture approximately 50 points per period on the 
    % faster time scale
    numberPeriodsFast = floor(endTime/fastTimeScale);
    numberPointsPerPeriod = 50;
    nSamp = numberPeriodsFast*numberPointsPerPeriod+1;
    dt = endTime/(nSamp-1);
    xData = integrateTrajectoriesGalphaDirect(Model, endTime, ICs, nSamp, observable);
    loadShape = f_0;
    DataInfo = struct('nElements', Model.Mesh.nElements, 'loadShape', loadShape);
    save('dataVKDecay2DGalphaModal.mat', 'DataInfo', 'xData', 'dt', 'endTime', 'nSamp')
else
    load("../../data/vonkarmanplateIR/dataVKDecay2DGalphaModalData2.mat")
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
fig = customFigure();
fontName = 'helvetica';
fontSize = 16;
tiledlayout(3,1);
nexttile
showSpectrogram(xData(1,:), outdof(1), 1);
ylim([25,200])
set(gca, 'fontname', fontName)
set(gca, 'fontsize', fontSize)
xticks([])
nexttile
showSpectrogram(xData(3,:), outdof(1), 1);
ylim([25,200])
set(gca, 'fontname', fontName)
set(gca, 'fontsize', fontSize)
ylabel('frequency [Hz]')
xticks([])
nexttile
showSpectrogram(xData(5,:), outdof(1), 1);
ylim([25,200])
xlabel('time [s]')
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'power spectral density [1/Hz]';
set(gca, 'fontname', fontName)
set(gca, 'fontsize', fontSize)
```

![figure_13.png
](README_images/figure_13.png
)

We plot the observables of interest over time for closer inspection. 

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

![figure_14.png
](README_images/figure_14.png
)

# Truncate transient data from trajectories

We must however remove the first transient to fulfill the assumption that trajectories lie close to the SSM. We keep only the time interval |sliceInt|.

```matlab:Code
sliceInt = [2*slowTimeScale, endTime];
xDataTrunc = sliceTrajectories(xData, sliceInt);
```

# Datadriven manifold fitting

The measured trajectories are initialized to lie close to the manifold of interest that is tangent at the origin to the eigenspace spanned by the columns of $V_e$. 

As we also know the projection $W_e$ to this eigenspace, we define the modal coordinates as $y=W_e x$. These are the reduced coordinates for our graph style parametrization of the manifold, gauranteed to exists near the origin. We then use the data to learn the nonlinear feature of the manifold geometry, represented via polynomials. Indeed, we seek the $2N\times m_M$ polynomial coefficient matrix $H$ of the manifold on the form

> $x=V_e y+H{{\phi }}_{m,2:M} (y)$,

where the function ${{\phi }}_{m,2:M} (y)$ computes a vector of all $m_M$ monomials from orders 2 up to $M$ of an $m$-vector $y$. From SSM theory, the tangent space of the manifold is $V_e$. The coefficients $\mathbf{H}$are obtained via least squares regression.

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

![figure_15.png
](README_images/figure_15.png
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
Reconstruction error = 0.29665%
```

```matlab:Code

% Plot physical coordinates
customFigure;
for iTraj = indTrain
    plot3(xData{iTraj,2}(outdof(1),:),xData{iTraj,2}(outdof(1)+size(M,1),:),xData{iTraj,2}(outdof(2),:))
end
xlabel('$q_{A}$ [m]','interpreter','latex')
ylabel('$\dot{q}_{A}$ [m/s]','interpreter','latex')
zlabel('$q_{B}$ [m]','interpreter','latex')
legend('Trajectory 1','Trajectory 2','Trajectory 3')
view(3)
```

![figure_16.png
](README_images/figure_16.png
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
ROMOrder = 5;
freqNorm = [1 1];

RDInfo = IMDynamicsMech(yDataTrunc(indTrain,:), ...
    'R_PolyOrd', 1,'N_PolyOrd', ROMOrder, 'style', 'normalform', ...
    'R_coeff',Ae,'rescale',1,'frequencies_norm',freqNorm,'MaxIter',5e3);
```

```text:Output
Estimation of the reduced dynamics...  Done. 
Estimation of the reduced dynamics in normal form...
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
     0           1          64.5768                           155
     1           3          40.5944     0.00122686            127  
     2           5          35.6058       0.482621            136  
     3           7          32.3071       0.349418           17.9  
     4           8          31.6722              1           16.6  
     5           9          29.4008              1           43.8  
     6          11          29.1837       0.219698           17.3  
     7          12          28.8577              1           11.8  
     8          13          28.5823              1           2.27  
     9          14          28.5556              1           2.26  
    10          15          28.4543              1           3.51  
    11          16          28.3995              1           3.22  
   ...  
   791         797          4.94294              1          0.155  
   792         798          4.93065              1          0.209  
   793         799          4.92435              1          0.185  
   794         800          4.92302              1          0.177  
   795         801           4.9227              1          0.169  
   796         802          4.92214              1          0.156  
   797         803          4.92069              1          0.134  
   798         804          4.91714              1         0.0949  
   799         805          4.90903              1         0.0903  
                       ...
```

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
Normalized mean trajectory error = 5.7168%
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
plotTrajectories(xData(indTest,:), yRec(indTest,:), 'm','PlotCoordinate',...
    outdof(1), 'DisplayName', {'Test set', 'Prediction'})
ylabel('$u \, [$m$]$','Interpreter','latex')
```

![figure_18.png
](README_images/figure_18.png
)

```matlab:Code

fig = customFigure('subPlot',[3 2]);
for iTraj = indTest
    subplot(3,2,iTraj-1)
    plot(yDataTrunc{iTraj,1},yDataTrunc{iTraj,2}(1,:),'k','Linewidth',1,'Color',colors(colFOM,:))
    plot(etaRec{iTraj,1},etaRec{iTraj,2}(1,:),':','Linewidth',1.5,'Color',colors(colSSML,:))
    xlabel('$t \, [$s$]$','Interpreter','latex'); 
    ylabel('$u_2$','Interpreter','latex'); xlim([0 2])
    if iTraj == indTest(2); ylim([-1 1]*max(abs(yDataTrunc{iTraj,2}(2,:)))); end
    subplot(3,2,iTraj)
    plot(yDataTrunc{iTraj,1},yDataTrunc{iTraj,2}(2,:),'k','Linewidth',1,'Color',colors(colFOM,:))
    plot(etaRec{iTraj,1},etaRec{iTraj,2}(2,:),':','Linewidth',1.5,'Color',colors(colSSML,:))
    xlabel('$t \, [$s$]$','Interpreter','latex'); 
    ylabel('$u_3$','Interpreter','latex'); xlim([0 2])
    if iTraj == indTest(1); ylim([-1 1]*max(abs(yDataTrunc{iTraj,2}(1,:)))); end
end
```

![figure_19.png
](README_images/figure_19.png
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
epsilon = [.2 .4];
mFreqs = [1 1];
resonantModes = [1 2 3 4];
omegaSpan =  [0.95 1.1]*imag(lambda(2));
% SSMLearn
[FRCSSMLearn] = continuationFRCep(IMInfoF, RDInfoF, epsilon, omegaSpan,@(x) x(outdof,:), mFreqs,resonantModes, 'SSMLearnvonKarmanPlateIR');
```

```text:Output
 Run='SSMLearnvonKarmanPlateIR.ep': Continue equilibria along primary branch.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          7.92e-03  1.08e+03    0.0    0.0    0.0
   1   1  1.00e+00  5.62e-03  3.62e-05  1.08e+03    0.0    0.0    0.0
   2   1  1.00e+00  2.90e-05  2.14e-10  1.08e+03    0.0    0.0    0.0
   3   1  1.00e+00  2.36e-10  1.20e-15  1.08e+03    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   1.0799e+03      1  EP      7.6359e+02   1.4414e-02   3.3645e-02   1.8524e-01   4.3780e+00   2.0000e-01
   10  00:00:00   1.0760e+03      2          7.6081e+02   7.6587e-03   2.2740e-02   1.0582e+00   4.4859e+00   2.0000e-01
   20  00:00:00   1.0710e+03      3          7.5729e+02   4.0024e-03   1.5684e-02   1.3180e+00   4.5572e+00   2.0000e-01
   30  00:00:00   1.0660e+03      4          7.5375e+02   2.6428e-03   1.1864e-02   1.4061e+00   4.5957e+00   2.0000e-01
   40  00:00:00   1.0610e+03      5          7.5022e+02   1.9632e-03   9.5177e-03   1.4494e+00   4.6193e+00   2.0000e-01
   50  00:00:00   1.0560e+03      6          7.4668e+02   1.5593e-03   7.9396e-03   1.4749e+00   4.6352e+00   2.0000e-01
   60  00:00:00   1.0510e+03      7          7.4315e+02   1.2924e-03   6.8079e-03   1.4917e+00   4.6465e+00   2.0000e-01
   70  00:00:01   1.0460e+03      8          7.3961e+02   1.1032e-03   5.9575e-03   1.5037e+00   4.6550e+00   2.0000e-01
   80  00:00:01   1.0410e+03      9          7.3608e+02   9.6224e-04   5.2954e-03   1.5125e+00   4.6616e+00   2.0000e-01
   90  00:00:01   1.0360e+03     10          7.3254e+02   8.5312e-04   4.7655e-03   1.5194e+00   4.6669e+00   2.0000e-01
  100  00:00:01   1.0310e+03     11          7.2900e+02   7.6620e-04   4.3318e-03   1.5249e+00   4.6713e+00   2.0000e-01
  110  00:00:01   1.0260e+03     12          7.2547e+02   6.9532e-04   3.9703e-03   1.5293e+00   4.6749e+00   2.0000e-01
  111  00:00:01   1.0259e+03     13  EP      7.2541e+02   6.9430e-04   3.9651e-03   1.5294e+00   4.6749e+00   2.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:01   1.0799e+03     14  EP      7.6359e+02   1.4414e-02   3.3645e-02   1.8524e-01   4.3780e+00   2.0000e-01
   10  00:00:01   1.0837e+03     15          7.6626e+02   1.0511e-02   5.1933e-02  -9.8992e-01   4.1644e+00   2.0000e-01
   20  00:00:01   1.0886e+03     16          7.6974e+02   5.2264e-03   8.2415e-02  -1.3523e+00   3.7411e+00   2.0000e-01
   30  00:00:01   1.0924e+03     17          7.7246e+02   3.6455e-03   9.9255e-02  -1.1458e+00   3.0464e+00   2.0000e-01
   38  00:00:02   1.0925e+03     18  FP      7.7250e+02   3.7252e-03   9.7795e-02  -1.1258e+00   2.9454e+00   2.0000e-01
   38  00:00:02   1.0925e+03     19  SN      7.7250e+02   3.7252e-03   9.7795e-02  -1.1258e+00   2.9454e+00   2.0000e-01
   40  00:00:02   1.0925e+03     20          7.7250e+02   3.7513e-03   9.7222e-02  -1.1227e+00   2.9178e+00   2.0000e-01
   50  00:00:02   1.0923e+03     21          7.7237e+02   3.9150e-03   9.2331e-02  -1.1249e+00   2.7541e+00   2.0000e-01
   60  00:00:02   1.0913e+03     22          7.7165e+02   3.8887e-03   6.2201e-02  -1.2665e+00   2.2444e+00   2.0000e-01
   62  00:00:02   1.0913e+03     23  SN      7.7165e+02   3.8715e-03   6.1349e-02  -1.2703e+00   2.2336e+00   2.0000e-01
   62  00:00:02   1.0913e+03     24  FP      7.7165e+02   3.8715e-03   6.1349e-02  -1.2703e+00   2.2336e+00   2.0000e-01
   70  00:00:02   1.0913e+03     25          7.7167e+02   3.7719e-03   5.6954e-02  -1.2893e+00   2.1789e+00   2.0000e-01
   80  00:00:02   1.0917e+03     26          7.7196e+02   3.4473e-03   4.6120e-02  -1.3325e+00   2.0522e+00   2.0000e-01
   90  00:00:02   1.0962e+03     27          7.7513e+02   2.3209e-03   2.2918e-02  -1.4215e+00   1.8041e+00   2.0000e-01
  100  00:00:02   1.1012e+03     28          7.7866e+02   1.7620e-03   1.5366e-02  -1.4576e+00   1.7271e+00   2.0000e-01
  110  00:00:03   1.1062e+03     29          7.8220e+02   1.4250e-03   1.1600e-02  -1.4790e+00   1.6892e+00   2.0000e-01
  120  00:00:03   1.1112e+03     30          7.8573e+02   1.1972e-03   9.3226e-03  -1.4934e+00   1.6663e+00   2.0000e-01
  130  00:00:03   1.1162e+03     31          7.8927e+02   1.0324e-03   7.7940e-03  -1.5038e+00   1.6509e+00   2.0000e-01
  140  00:00:03   1.1212e+03     32          7.9280e+02   9.0757e-04   6.6964e-03  -1.5116e+00   1.6399e+00   2.0000e-01
  150  00:00:03   1.1262e+03     33          7.9634e+02   8.0969e-04   5.8698e-03  -1.5178e+00   1.6317e+00   2.0000e-01
  160  00:00:03   1.1312e+03     34          7.9988e+02   7.3087e-04   5.2249e-03  -1.5227e+00   1.6252e+00   2.0000e-01
  170  00:00:03   1.1362e+03     35          8.0341e+02   6.6604e-04   4.7076e-03  -1.5268e+00   1.6200e+00   2.0000e-01
  180  00:00:03   1.1412e+03     36          8.0695e+02   6.1177e-04   4.2835e-03  -1.5302e+00   1.6158e+00   2.0000e-01
  190  00:00:03   1.1462e+03     37          8.1048e+02   5.6568e-04   3.9295e-03  -1.5331e+00   1.6123e+00   2.0000e-01
  200  00:00:04   1.1512e+03     38          8.1402e+02   5.2605e-04   3.6295e-03  -1.5356e+00   1.6093e+00   2.0000e-01
  210  00:00:04   1.1562e+03     39          8.1755e+02   4.9161e-04   3.3721e-03  -1.5377e+00   1.6067e+00   2.0000e-01
  220  00:00:04   1.1612e+03     40          8.2109e+02   4.6140e-04   3.1487e-03  -1.5396e+00   1.6045e+00   2.0000e-01
  230  00:00:04   1.1662e+03     41          8.2462e+02   4.3468e-04   2.9531e-03  -1.5413e+00   1.6025e+00   2.0000e-01
  240  00:00:04   1.1712e+03     42          8.2816e+02   4.1089e-04   2.7804e-03  -1.5428e+00   1.6008e+00   2.0000e-01
  250  00:00:04   1.1762e+03     43          8.3170e+02   3.8957e-04   2.6268e-03  -1.5442e+00   1.5992e+00   2.0000e-01
  260  00:00:04   1.1812e+03     44          8.3523e+02   3.7035e-04   2.4892e-03  -1.5454e+00   1.5979e+00   2.0000e-01
  270  00:00:04   1.1862e+03     45          8.3877e+02   3.5294e-04   2.3654e-03  -1.5465e+00   1.5966e+00   2.0000e-01
  274  00:00:04   1.1879e+03     46  EP      8.3995e+02   3.4747e-04   2.3266e-03  -1.5468e+00   1.5962e+00   2.0000e-01
```

![figure_20.png
](README_images/figure_20.png
)

![figure_21.png
](README_images/figure_21.png
)

```text:Output
 Run='SSMLearnvonKarmanPlateIR.ep.ep': Continue equilibria along primary branch.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          5.17e-02  1.08e+03    0.0    0.0    0.0
   1   1  1.00e+00  4.81e-03  1.44e-04  1.08e+03    0.0    0.0    0.0
   2   1  1.00e+00  5.84e-05  1.18e-09  1.08e+03    0.0    0.0    0.0
   3   1  1.00e+00  7.08e-10  1.69e-15  1.08e+03    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   1.0799e+03      1  EP      7.6359e+02   2.2158e-02   5.6418e-02   4.7644e-01   4.4423e+00   4.0000e-01
   10  00:00:00   1.0759e+03      2          7.6074e+02   1.3092e-02   4.1754e-02   1.1176e+00   4.5057e+00   4.0000e-01
   20  00:00:00   1.0709e+03      3          7.5721e+02   7.5834e-03   3.0256e-02   1.3294e+00   4.5629e+00   4.0000e-01
   30  00:00:00   1.0659e+03      4          7.5368e+02   5.1563e-03   2.3290e-02   1.4098e+00   4.5979e+00   4.0000e-01
   40  00:00:00   1.0609e+03      5          7.5014e+02   3.8716e-03   1.8822e-02   1.4510e+00   4.6204e+00   4.0000e-01
   50  00:00:00   1.0559e+03      6          7.4661e+02   3.0900e-03   1.5758e-02   1.4758e+00   4.6358e+00   4.0000e-01
   60  00:00:00   1.0509e+03      7          7.4307e+02   2.5678e-03   1.3539e-02   1.4923e+00   4.6469e+00   4.0000e-01
   70  00:00:00   1.0459e+03      8          7.3953e+02   2.1954e-03   1.1862e-02   1.5040e+00   4.6553e+00   4.0000e-01
   80  00:00:01   1.0409e+03      9          7.3600e+02   1.9167e-03   1.0552e-02   1.5128e+00   4.6618e+00   4.0000e-01
   90  00:00:01   1.0359e+03     10          7.3246e+02   1.7005e-03   9.5019e-03   1.5196e+00   4.6671e+00   4.0000e-01
  100  00:00:01   1.0309e+03     11          7.2893e+02   1.5280e-03   8.6408e-03   1.5250e+00   4.6714e+00   4.0000e-01
  110  00:00:01   1.0259e+03     12  EP      7.2541e+02   1.3878e-03   7.9260e-03   1.5294e+00   4.6750e+00   4.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:01   1.0799e+03     13  EP      7.6359e+02   2.2158e-02   5.6418e-02   4.7644e-01   4.4423e+00   4.0000e-01
   10  00:00:01   1.0836e+03     14          7.6623e+02   3.3209e-02   7.3624e-02  -7.9409e-01   4.3125e+00   4.0000e-01
   20  00:00:01   1.0885e+03     15          7.6969e+02   1.9248e-02   1.0045e-01  -1.4624e+00   4.1689e+00   4.0000e-01
   30  00:00:01   1.0935e+03     16          7.7321e+02   1.2230e-02   1.2677e-01  -1.6593e+00   4.0134e+00   4.0000e-01
   40  00:00:01   1.0985e+03     17          7.7674e+02   8.4958e-03   1.5010e-01  -1.7615e+00   3.8472e+00   4.0000e-01
   50  00:00:01   1.1035e+03     18          7.8027e+02   5.9079e-03   1.7080e-01  -1.7789e+00   3.6597e+00   4.0000e-01
   60  00:00:02   1.1085e+03     19          7.8379e+02   3.9690e-03   1.8897e-01  -1.5871e+00   3.4050e+00   4.0000e-01
   70  00:00:02   1.1108e+03     20          7.8542e+02   3.7724e-03   1.9523e-01  -1.1608e+00   3.1086e+00   4.0000e-01
   74  00:00:02   1.1108e+03     21  SN      7.8542e+02   3.8231e-03   1.9509e-01  -1.1398e+00   3.0900e+00   4.0000e-01
   74  00:00:02   1.1108e+03     22  FP      7.8542e+02   3.8231e-03   1.9509e-01  -1.1398e+00   3.0900e+00   4.0000e-01
   80  00:00:02   1.1107e+03     23          7.8541e+02   3.9099e-03   1.9475e-01  -1.1105e+00   3.0619e+00   4.0000e-01
   90  00:00:02   1.1105e+03     24          7.8524e+02   4.2011e-03   1.9307e-01  -1.0470e+00   2.9843e+00   4.0000e-01
  100  00:00:02   1.1067e+03     25          7.8255e+02   5.8058e-03   1.7344e-01  -1.0032e+00   2.6513e+00   4.0000e-01
  110  00:00:02   1.1017e+03     26          7.7903e+02   6.7861e-03   1.4495e-01  -1.1109e+00   2.3933e+00   4.0000e-01
  120  00:00:02   1.0967e+03     27          7.7551e+02   6.8438e-03   1.0411e-01  -1.2635e+00   2.1225e+00   4.0000e-01
  129  00:00:03   1.0953e+03     28  FP      7.7450e+02   6.0074e-03   7.3141e-02  -1.3541e+00   1.9478e+00   4.0000e-01
  129  00:00:03   1.0953e+03     29  SN      7.7450e+02   6.0074e-03   7.3141e-02  -1.3541e+00   1.9478e+00   4.0000e-01
  130  00:00:03   1.0953e+03     30          7.7450e+02   5.9407e-03   7.1465e-02  -1.3583e+00   1.9388e+00   4.0000e-01
  140  00:00:03   1.0956e+03     31          7.7474e+02   5.4466e-03   6.0484e-02  -1.3852e+00   1.8804e+00   4.0000e-01
  150  00:00:03   1.1003e+03     32          7.7806e+02   3.7724e-03   3.3882e-02  -1.4491e+00   1.7432e+00   4.0000e-01
  160  00:00:03   1.1053e+03     33          7.8159e+02   2.9790e-03   2.4591e-02  -1.4748e+00   1.6962e+00   4.0000e-01
  170  00:00:03   1.1103e+03     34          7.8513e+02   2.4760e-03   1.9443e-02  -1.4908e+00   1.6703e+00   4.0000e-01
  180  00:00:03   1.1153e+03     35          7.8866e+02   2.1219e-03   1.6109e-02  -1.5020e+00   1.6536e+00   4.0000e-01
  190  00:00:03   1.1203e+03     36          7.9220e+02   1.8576e-03   1.3762e-02  -1.5103e+00   1.6418e+00   4.0000e-01
  200  00:00:04   1.1253e+03     37          7.9573e+02   1.6523e-03   1.2015e-02  -1.5167e+00   1.6331e+00   4.0000e-01
  210  00:00:04   1.1303e+03     38          7.9927e+02   1.4881e-03   1.0664e-02  -1.5219e+00   1.6263e+00   4.0000e-01
  220  00:00:04   1.1353e+03     39          8.0280e+02   1.3537e-03   9.5868e-03  -1.5261e+00   1.6209e+00   4.0000e-01
  230  00:00:04   1.1403e+03     40          8.0634e+02   1.2416e-03   8.7076e-03  -1.5296e+00   1.6165e+00   4.0000e-01
  240  00:00:04   1.1453e+03     41          8.0988e+02   1.1467e-03   7.9763e-03  -1.5326e+00   1.6128e+00   4.0000e-01
  250  00:00:04   1.1503e+03     42          8.1341e+02   1.0653e-03   7.3585e-03  -1.5352e+00   1.6098e+00   4.0000e-01
  260  00:00:04   1.1553e+03     43          8.1695e+02   9.9467e-04   6.8296e-03  -1.5374e+00   1.6071e+00   4.0000e-01
  270  00:00:04   1.1603e+03     44          8.2048e+02   9.3285e-04   6.3716e-03  -1.5393e+00   1.6048e+00   4.0000e-01
  280  00:00:04   1.1653e+03     45          8.2402e+02   8.7826e-04   5.9713e-03  -1.5410e+00   1.6028e+00   4.0000e-01
  290  00:00:05   1.1703e+03     46          8.2755e+02   8.2971e-04   5.6183e-03  -1.5426e+00   1.6011e+00   4.0000e-01
  300  00:00:05   1.1753e+03     47          8.3109e+02   7.8625e-04   5.3047e-03  -1.5439e+00   1.5995e+00   4.0000e-01
  310  00:00:05   1.1803e+03     48          8.3462e+02   7.4712e-04   5.0243e-03  -1.5452e+00   1.5981e+00   4.0000e-01
  320  00:00:05   1.1853e+03     49          8.3816e+02   7.1170e-04   4.7720e-03  -1.5463e+00   1.5968e+00   4.0000e-01
  326  00:00:05   1.1879e+03     50  EP      8.3995e+02   6.9500e-04   4.6536e-03  -1.5468e+00   1.5962e+00   4.0000e-01
```

![figure_22.png
](README_images/figure_22.png
)

![figure_23.png
](README_images/figure_23.png
)

```matlab:Code
% SSMTool
[FRCSSMTool] = SSMToolFRCFE(M,C,K,fnl,M*forcingVectors(n+1:end),outdof,epsilon,[3:6],5,omegaSpan,mFreqs, 'SSMToolvonKarmanPlateIR');
```

```text:Output
Due to high-dimensionality, we compute only the first 5 eigenvalues with the smallest magnitude. These would also be used to compute the spectral quotients
Assuming a proportional damping hypthesis with symmetric matrices
modal damping ratio for 1 mode is 2.243693e-03
modal damping ratio for 2 mode is 2.181986e-03
modal damping ratio for 3 mode is 2.186754e-03
modal damping ratio for 4 mode is 2.847788e-03
modal damping ratio for 5 mode is 3.388843e-03
the left eigenvectors may be incorrect in case of asymmetry of matrices

 The first 10 nonzero eigenvalues are given as 
     -0.68812 +     306.69i
     -0.68812 -     306.69i
      -1.6662 +     763.59i
      -1.6662 -     763.59i
      -1.6789 +     767.75i
      -1.6789 -     767.75i
      -3.4708 +     1218.8i
      -3.4708 -     1218.8i
      -5.1888 +     1531.1i
      -5.1888 -     1531.1i

(near) outer resonance detected for the following combination of master eigenvalues
     0     0     1     1
     0     1     1     0
     1     0     0     1
     1     1     0     0
     0     0     1     1
     0     1     1     0
     1     0     0     1
     1     1     0     0
     0     0     2     0
     1     0     1     0
     2     0     0     0
     0     0     0     2
     0     1     0     1
     0     2     0     0
     0     0     2     0
     1     0     1     0
     2     0     0     0
     0     0     0     2
     0     1     0     1
     0     2     0     0

These are in resonance with the follwing eigenvalues of the slave subspace
     -0.68812 +     306.69i
     -0.68812 +     306.69i
     -0.68812 +     306.69i
     -0.68812 +     306.69i
     -0.68812 -     306.69i
     -0.68812 -     306.69i
     -0.68812 -     306.69i
     -0.68812 -     306.69i
      -3.4708 +     1218.8i
      -3.4708 +     1218.8i
      -3.4708 +     1218.8i
      -3.4708 -     1218.8i
      -3.4708 -     1218.8i
      -3.4708 -     1218.8i
      -5.1888 +     1531.1i
      -5.1888 +     1531.1i
      -5.1888 +     1531.1i
      -5.1888 -     1531.1i
      -5.1888 -     1531.1i
      -5.1888 -     1531.1i

sigma_out = 3
(near) inner resonance detected for the following combination of master eigenvalues
     0     0     2     1
     0     1     2     0
     1     0     1     1
     1     1     1     0
     2     0     0     1
     2     1     0     0
     0     0     1     2
     0     1     1     1
     0     2     1     0
     1     0     0     2
     1     1     0     1
     1     2     0     0
     0     0     2     1
     0     1     2     0
     1     0     1     1
     1     1     1     0
     2     0     0     1
     2     1     0     0
     0     0     1     2
     0     1     1     1
     0     2     1     0
     1     0     0     2
     1     1     0     1
     1     2     0     0

These are in resonance with the follwing eigenvalues of the master subspace
      -1.6662 +     763.59i
      -1.6662 +     763.59i
      -1.6662 +     763.59i
      -1.6662 +     763.59i
      -1.6662 +     763.59i
      -1.6662 +     763.59i
      -1.6662 -     763.59i
      -1.6662 -     763.59i
      -1.6662 -     763.59i
      -1.6662 -     763.59i
      -1.6662 -     763.59i
      -1.6662 -     763.59i
      -1.6789 +     767.75i
      -1.6789 +     767.75i
      -1.6789 +     767.75i
      -1.6789 +     767.75i
      -1.6789 +     767.75i
      -1.6789 +     767.75i
      -1.6789 -     767.75i
      -1.6789 -     767.75i
      -1.6789 -     767.75i
      -1.6789 -     767.75i
      -1.6789 -     767.75i
      -1.6789 -     767.75i

sigma_in = 3
Due to (near) outer resonance, the exisitence of the manifold is questionable and the underlying computation may suffer.
Attempting manifold computation
Manifold computation time at order 2 = 00:00:00
Estimated memory usage at order  2 = 2.89E+00 MB
Manifold computation time at order 3 = 00:00:00
Estimated memory usage at order  3 = 4.65E+00 MB
Manifold computation time at order 4 = 00:00:01
Estimated memory usage at order  4 = 9.22E+00 MB
Manifold computation time at order 5 = 00:00:03
Estimated memory usage at order  5 = 1.65E+01 MB

 Run='SSMToolvonKarmanPlateIReps.ep': Continue equilibria with varied epsilon.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          1.35e-01  7.64e+02    0.0    0.0    0.0
   1   1  1.00e+00  2.66e-02  1.77e-03  7.64e+02    0.0    0.0    0.0
   2   1  1.00e+00  3.45e-04  3.00e-07  7.64e+02    0.0    0.0    0.0
   3   1  1.00e+00  3.77e-08  1.07e-14  7.64e+02    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE           eps         rho1         rho2          th1          th2           om
    0  00:00:00   7.6367e+02      1  EP      2.0000e-01   2.9964e-03   6.9387e-03   4.9180e+00   5.9505e+00   7.6359e+02
    1  00:00:00   7.6367e+02      2  UZ      2.0000e-01   2.9964e-03   6.9387e-03   4.9180e+00   5.9505e+00   7.6359e+02
    1  00:00:00   7.6367e+02      3  EP      1.8000e-01   2.7577e-03   6.3612e-03   4.8856e+00   5.9430e+00   7.6359e+02

 STEP      TIME        ||U||  LABEL  TYPE           eps         rho1         rho2          th1          th2           om
    0  00:00:00   7.6367e+02      4  EP      2.0000e-01   2.9964e-03   6.9387e-03   4.9180e+00   5.9505e+00   7.6359e+02
    4  00:00:00   7.6368e+02      5  UZ      4.0000e-01   4.4916e-03   1.1520e-02   5.2282e+00   6.0162e+00   7.6359e+02
    5  00:00:00   7.6368e+02      6  EP      4.4000e-01   4.6544e-03   1.2246e-02   5.2799e+00   6.0262e+00   7.6359e+02

 Run='SSMToolvonKarmanPlateIReps1.ep': Continue equilibria with varied omega at eps equal to 2.000000e-01.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          1.07e-14  1.08e+03    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   1.0799e+03      1  EP      7.6359e+02   2.9964e-03   6.9387e-03   4.9180e+00   5.9505e+00   2.0000e-01
   10  00:00:00   1.0760e+03      2          7.6080e+02   1.5957e-03   4.7220e-03   5.7727e+00   6.0558e+00   2.0000e-01
   20  00:00:00   1.0710e+03      3          7.5728e+02   8.3992e-04   3.2666e-03   6.0289e+00   6.1262e+00   2.0000e-01
   30  00:00:00   1.0660e+03      4          7.5375e+02   5.5565e-04   2.4732e-03   6.1165e+00   6.1645e+00   2.0000e-01
   40  00:00:00   1.0610e+03      5          7.5021e+02   4.1304e-04   1.9849e-03   6.1597e+00   6.1880e+00   2.0000e-01
   50  00:00:00   1.0560e+03      6          7.4668e+02   3.2816e-04   1.6561e-03   6.1852e+00   6.2038e+00   2.0000e-01
   60  00:00:00   1.0510e+03      7          7.4314e+02   2.7205e-04   1.4202e-03   6.2020e+00   6.2151e+00   2.0000e-01
   70  00:00:00   1.0460e+03      8          7.3960e+02   2.3225e-04   1.2429e-03   6.2139e+00   6.2236e+00   2.0000e-01
   80  00:00:00   1.0410e+03      9          7.3607e+02   2.0258e-04   1.1048e-03   6.2228e+00   6.2303e+00   2.0000e-01
   90  00:00:00   1.0360e+03     10          7.3253e+02   1.7962e-04   9.9429e-04   6.2296e+00   6.2356e+00   2.0000e-01
  100  00:00:01   1.0310e+03     11          7.2900e+02   1.6132e-04   9.0382e-04   6.2351e+00   6.2399e+00   2.0000e-01
  110  00:00:01   1.0260e+03     12          7.2546e+02   1.4640e-04   8.2843e-04   6.2395e+00   6.2435e+00   2.0000e-01
  111  00:00:01   1.0260e+03     13  EP      7.2541e+02   1.4621e-04   8.2746e-04   6.2396e+00   6.2436e+00   2.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:01   1.0799e+03     14  EP      7.6359e+02   2.9964e-03   6.9387e-03   4.9180e+00   5.9505e+00   2.0000e-01
   10  00:00:01   1.0837e+03     15          7.6625e+02   2.3057e-03   1.0583e-02   3.7272e+00   5.7460e+00   2.0000e-01
   20  00:00:01   1.0886e+03     16          7.6974e+02   1.1758e-03   1.6617e-02   3.3192e+00   5.3576e+00   2.0000e-01
   30  00:00:01   1.0931e+03     17          7.7289e+02   7.2422e-04   2.0837e-02   3.5177e+00   4.6859e+00   2.0000e-01
   40  00:00:01   1.0932e+03     18  FP      7.7300e+02   7.4644e-04   2.0520e-02   3.5749e+00   4.5337e+00   2.0000e-01
   40  00:00:01   1.0932e+03     19  SN      7.7300e+02   7.4644e-04   2.0520e-02   3.5749e+00   4.5337e+00   2.0000e-01
   40  00:00:01   1.0932e+03     20          7.7300e+02   7.4815e-04   2.0491e-02   3.5768e+00   4.5259e+00   2.0000e-01
   50  00:00:01   1.0931e+03     21          7.7293e+02   7.7707e-04   1.9909e-02   3.5939e+00   4.4094e+00   2.0000e-01
   60  00:00:01   1.0917e+03     22          7.7192e+02   8.3040e-04   1.3891e-02   3.4720e+00   3.8694e+00   2.0000e-01
   68  00:00:02   1.0916e+03     23  SN      7.7185e+02   8.0096e-04   1.2249e-02   3.4346e+00   3.7684e+00   2.0000e-01
   68  00:00:02   1.0916e+03     24  FP      7.7185e+02   8.0096e-04   1.2249e-02   3.4346e+00   3.7684e+00   2.0000e-01
   70  00:00:02   1.0916e+03     25          7.7185e+02   7.9530e-04   1.2001e-02   3.4293e+00   3.7538e+00   2.0000e-01
   80  00:00:02   1.0917e+03     26          7.7191e+02   7.6468e-04   1.0838e-02   3.4050e+00   3.6872e+00   2.0000e-01
   90  00:00:02   1.0937e+03     27          7.7334e+02   5.9280e-04   6.5481e-03   3.3240e+00   3.4605e+00   2.0000e-01
  100  00:00:02   1.0987e+03     28          7.7688e+02   4.2268e-04   3.8503e-03   3.2684e+00   3.3270e+00   2.0000e-01
  110  00:00:02   1.1037e+03     29          7.8041e+02   3.3236e-04   2.7646e-03   3.2410e+00   3.2744e+00   2.0000e-01
  120  00:00:02   1.1087e+03     30          7.8395e+02   2.7437e-04   2.1605e-03   3.2235e+00   3.2452e+00   2.0000e-01
  130  00:00:02   1.1137e+03     31          7.8748e+02   2.3372e-04   1.7738e-03   3.2113e+00   3.2266e+00   2.0000e-01
  140  00:00:02   1.1187e+03     32          7.9102e+02   2.0359e-04   1.5047e-03   3.2023e+00   3.2137e+00   2.0000e-01
  150  00:00:02   1.1237e+03     33          7.9455e+02   1.8036e-04   1.3066e-03   3.1954e+00   3.2042e+00   2.0000e-01
  160  00:00:02   1.1287e+03     34          7.9809e+02   1.6189e-04   1.1546e-03   3.1899e+00   3.1969e+00   2.0000e-01
  170  00:00:03   1.1337e+03     35          8.0162e+02   1.4685e-04   1.0342e-03   3.1854e+00   3.1911e+00   2.0000e-01
  180  00:00:03   1.1387e+03     36          8.0516e+02   1.3437e-04   9.3660e-04   3.1817e+00   3.1865e+00   2.0000e-01
  190  00:00:03   1.1437e+03     37          8.0870e+02   1.2385e-04   8.5581e-04   3.1785e+00   3.1826e+00   2.0000e-01
  200  00:00:03   1.1487e+03     38          8.1223e+02   1.1485e-04   7.8786e-04   3.1758e+00   3.1793e+00   2.0000e-01
  210  00:00:03   1.1537e+03     39          8.1577e+02   1.0708e-04   7.2990e-04   3.1735e+00   3.1765e+00   2.0000e-01
  220  00:00:03   1.1587e+03     40          8.1930e+02   1.0028e-04   6.7988e-04   3.1715e+00   3.1742e+00   2.0000e-01
  230  00:00:03   1.1637e+03     41          8.2284e+02   9.4303e-05   6.3628e-04   3.1697e+00   3.1721e+00   2.0000e-01
  240  00:00:03   1.1687e+03     42          8.2637e+02   8.8995e-05   5.9793e-04   3.1681e+00   3.1702e+00   2.0000e-01
  250  00:00:03   1.1737e+03     43          8.2991e+02   8.4252e-05   5.6394e-04   3.1667e+00   3.1686e+00   2.0000e-01
  260  00:00:03   1.1787e+03     44          8.3344e+02   7.9990e-05   5.3360e-04   3.1654e+00   3.1671e+00   2.0000e-01
  270  00:00:04   1.1837e+03     45          8.3698e+02   7.6138e-05   5.0637e-04   3.1643e+00   3.1658e+00   2.0000e-01
  279  00:00:04   1.1879e+03     46  EP      8.3995e+02   7.3176e-05   4.8554e-04   3.1634e+00   3.1648e+00   2.0000e-01

 Run='SSMToolvonKarmanPlateIReps2.ep': Continue equilibria with varied omega at eps equal to 4.000000e-01.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          7.96e-13  1.08e+03    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   1.0799e+03      1  EP      7.6359e+02   4.4916e-03   1.1520e-02   5.2282e+00   6.0162e+00   4.0000e-01
   10  00:00:00   1.0759e+03      2          7.6073e+02   2.6889e-03   8.6110e-03   5.8365e+00   6.0770e+00   4.0000e-01
   20  00:00:00   1.0709e+03      3          7.5720e+02   1.5817e-03   6.2800e-03   6.0412e+00   6.1324e+00   4.0000e-01
   30  00:00:00   1.0659e+03      4          7.5367e+02   1.0811e-03   4.8468e-03   6.1205e+00   6.1669e+00   4.0000e-01
   40  00:00:00   1.0609e+03      5          7.5013e+02   8.1336e-04   3.9214e-03   6.1614e+00   6.1892e+00   4.0000e-01
   50  00:00:00   1.0559e+03      6          7.4660e+02   6.4976e-04   3.2849e-03   6.1861e+00   6.2045e+00   4.0000e-01
   60  00:00:00   1.0509e+03      7          7.4306e+02   5.4021e-04   2.8232e-03   6.2025e+00   6.2155e+00   4.0000e-01
   70  00:00:00   1.0459e+03      8          7.3953e+02   4.6198e-04   2.4740e-03   6.2142e+00   6.2239e+00   4.0000e-01
   80  00:00:00   1.0409e+03      9          7.3599e+02   4.0340e-04   2.2011e-03   6.2230e+00   6.2305e+00   4.0000e-01
   90  00:00:00   1.0359e+03     10          7.3245e+02   3.5794e-04   1.9822e-03   6.2298e+00   6.2357e+00   4.0000e-01
  100  00:00:00   1.0309e+03     11          7.2892e+02   3.2165e-04   1.8026e-03   6.2352e+00   6.2400e+00   4.0000e-01
  110  00:00:01   1.0260e+03     12  EP      7.2541e+02   2.9224e-04   1.6539e-03   6.2396e+00   6.2436e+00   4.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:01   1.0799e+03     13  EP      7.6359e+02   4.4916e-03   1.1520e-02   5.2282e+00   6.0162e+00   4.0000e-01
   10  00:00:01   1.0837e+03     14          7.6625e+02   7.2300e-03   1.4890e-02   3.9820e+00   5.8905e+00   4.0000e-01
   20  00:00:01   1.0885e+03     15          7.6968e+02   4.7207e-03   2.0124e-02   3.2418e+00   5.7533e+00   4.0000e-01
   30  00:00:01   1.0935e+03     16          7.7321e+02   3.0531e-03   2.5409e-02   3.0185e+00   5.6121e+00   4.0000e-01
   40  00:00:01   1.0985e+03     17          7.7674e+02   2.1797e-03   3.0124e-02   2.8804e+00   5.4617e+00   4.0000e-01
   50  00:00:01   1.1035e+03     18          7.8027e+02   1.5592e-03   3.4370e-02   2.7940e+00   5.2981e+00   4.0000e-01
   60  00:00:01   1.1085e+03     19          7.8380e+02   1.0332e-03   3.8211e-02   2.8058e+00   5.1005e+00   4.0000e-01
   70  00:00:01   1.1128e+03     20          7.8683e+02   6.5199e-04   4.1090e-02   3.2290e+00   4.8153e+00   4.0000e-01
   80  00:00:01   1.1134e+03     21          7.8727e+02   6.8186e-04   4.1286e-02   3.5554e+00   4.6626e+00   4.0000e-01
   81  00:00:01   1.1134e+03     22  FP      7.8727e+02   6.8306e-04   4.1283e-02   3.5585e+00   4.6608e+00   4.0000e-01
   81  00:00:01   1.1134e+03     23  SN      7.8727e+02   6.8306e-04   4.1283e-02   3.5585e+00   4.6608e+00   4.0000e-01
   90  00:00:01   1.1133e+03     24          7.8723e+02   7.1777e-04   4.1158e-02   3.6289e+00   4.6161e+00   4.0000e-01
  100  00:00:02   1.1127e+03     25          7.8680e+02   8.3903e-04   4.0451e-02   3.7460e+00   4.4995e+00   4.0000e-01
  110  00:00:02   1.1080e+03     26          7.8346e+02   1.2271e-03   3.5726e-02   3.7485e+00   4.1776e+00   4.0000e-01
  120  00:00:02   1.1030e+03     27          7.7993e+02   1.4222e-03   3.0158e-02   3.6295e+00   3.9518e+00   4.0000e-01
  130  00:00:02   1.0980e+03     28          7.7641e+02   1.4546e-03   2.2721e-02   3.4797e+00   3.7175e+00   4.0000e-01
  140  00:00:02   1.0957e+03     29          7.7478e+02   1.2507e-03   1.4850e-02   3.3569e+00   3.5052e+00   4.0000e-01
  141  00:00:02   1.0957e+03     30  SN      7.7478e+02   1.2449e-03   1.4711e-02   3.3551e+00   3.5016e+00   4.0000e-01
  141  00:00:02   1.0957e+03     31  FP      7.7478e+02   1.2449e-03   1.4710e-02   3.3551e+00   3.5016e+00   4.0000e-01
  150  00:00:02   1.0959e+03     32          7.7489e+02   1.1662e-03   1.2972e-02   3.3330e+00   3.4575e+00   4.0000e-01
  160  00:00:02   1.0999e+03     33          7.7774e+02   8.1970e-04   7.3877e-03   3.2658e+00   3.3194e+00   4.0000e-01
  170  00:00:02   1.1049e+03     34          7.8128e+02   6.4054e-04   5.2713e-03   3.2376e+00   3.2681e+00   4.0000e-01
  180  00:00:03   1.1099e+03     35          7.8482e+02   5.2992e-04   4.1389e-03   3.2208e+00   3.2408e+00   4.0000e-01
  190  00:00:03   1.1149e+03     36          7.8835e+02   4.5288e-04   3.4159e-03   3.2092e+00   3.2235e+00   4.0000e-01
  200  00:00:03   1.1199e+03     37          7.9189e+02   3.9573e-04   2.9107e-03   3.2006e+00   3.2113e+00   4.0000e-01
  210  00:00:03   1.1249e+03     38          7.9542e+02   3.5151e-04   2.5367e-03   3.1940e+00   3.2024e+00   4.0000e-01
  220  00:00:03   1.1299e+03     39          7.9896e+02   3.1623e-04   2.2483e-03   3.1888e+00   3.1954e+00   4.0000e-01
  230  00:00:03   1.1349e+03     40          8.0249e+02   2.8742e-04   2.0190e-03   3.1845e+00   3.1899e+00   4.0000e-01
  240  00:00:03   1.1399e+03     41          8.0603e+02   2.6344e-04   1.8322e-03   3.1809e+00   3.1855e+00   4.0000e-01
  250  00:00:03   1.1449e+03     42          8.0956e+02   2.4315e-04   1.6772e-03   3.1778e+00   3.1818e+00   4.0000e-01
  260  00:00:03   1.1499e+03     43          8.1310e+02   2.2578e-04   1.5463e-03   3.1753e+00   3.1786e+00   4.0000e-01
  270  00:00:03   1.1549e+03     44          8.1663e+02   2.1072e-04   1.4344e-03   3.1730e+00   3.1759e+00   4.0000e-01
  280  00:00:03   1.1599e+03     45          8.2017e+02   1.9755e-04   1.3377e-03   3.1710e+00   3.1736e+00   4.0000e-01
  290  00:00:04   1.1649e+03     46          8.2371e+02   1.8593e-04   1.2531e-03   3.1693e+00   3.1716e+00   4.0000e-01
  300  00:00:04   1.1699e+03     47  EP      8.2724e+02   1.7560e-04   1.1787e-03   3.1678e+00   3.1698e+00   4.0000e-01
Calculate FRC in physical domain at epsilon 2.000000e-01
the forcing frequency 7.2541e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the for...
```

![figure_24.png
](README_images/figure_24.png
)

![figure_25.png
](README_images/figure_25.png
)

![figure_26.png
](README_images/figure_26.png
)

![figure_27.png
](README_images/figure_27.png
)

![figure_28.png
](README_images/figure_28.png
)

![figure_29.png
](README_images/figure_29.png
)

![figure_30.png
](README_images/figure_30.png
)

![figure_31.png
](README_images/figure_31.png
)

# **Plot results**

```matlab:Code
fig = customFigure('subPlot',[2 1]);
subplot(211)
plotFRC(FRCSSMLearn, colors(colSSML,:), '- SSMLearn','freqscale',2*pi)
plotFRC(FRCSSMTool, colors(colSSMT,:), '- SSMTool','freqscale',2*pi)
xlabel('$\Omega$ [Hz]','interpreter','latex')
ylabel('amp($q_{A}$) [m]','interpreter','latex')
legend('off')
xlim([ceil(omegaSpan(1)/2/pi) floor(omegaSpan(2)/2/pi)-1])
ylim([0 0.0021])
subplot(212)
plotFRC(FRCSSMLearn, colors(colSSML,:), '- SSMLearn','outamp',2,'freqscale',2*pi)
plotFRC(FRCSSMTool, colors(colSSMT,:), '- SSMTool','outamp',2,'freqscale',2*pi)
xlabel('$\Omega$ [Hz]','interpreter','latex')
ylabel('amp($q_{B}$) [m]','interpreter','latex')
xlim([ceil(omegaSpan(1)/2/pi) floor(omegaSpan(2)/2/pi)-1])
ylim([0 0.00042])

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
subplot(2,1,2)
legend("Position", [0.53837,0.21482,0.44776,0.27273])
```

![figure_32.png
](README_images/figure_32.png
)

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
    18          21      0.000873128              1         0.0127  
    19          22      0.000545926              1         0.0156  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    20          23       0.00022258              1         0.0125  
    21          24      5.47161e-05              1        0.00521  
    22          25      2.02125e-05              1        0.00125  
    23          26      1.77684e-05              1       0.000775  
    24          27      1.75489e-05              1       0.000138  
    25          28      1.75277e-05              1       0.000224  
    26          29      1.75022e-05              1       9.35e-05  
    27          30      1.74774e-05              1       8.14e-05  
    28          31      1.73876e-05              1       5.23e-05  
    29          32       1.7376e-05              1       2.67e-05  
    30          33      1.73712e-05              1       2.53e-05  
    31          34      1.73677e-05              1       2.53e-05  
    32          35      1.73556e-05              1       5.55e-05  
    33          36       1.7332e-05              1       0.000106  
    34          37      1.72838e-05              1       0.000162  
    35          38      1.72194e-05              1        0.00017  
    36          39      1.71705e-05              1       0.000101  
    37          40      1.71531e-05              1       5.52e-05  
    38          41      1.71474e-05              1       5.28e-05  
    39          42      1.71402e-05              1       4.91e-05  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    40          43      1.71228e-05              1       9.05e-05  
    41          44      1.70838e-05              1       0.000142  
    42          45       1.7003e-05              1       0.000203  
    43          46       1.6875e-05              1       0.000261  
    44          47       1.6749e-05              1        0.00022  
    45          48      1.66872e-05              1         0.0001  
    46          49      1.66697e-05              1       5.71e-05  
    47          50      1.66602e-05              1       5.61e-05  
    48          51      1.66415e-05              1       9.99e-05  
    49          52      1.66012e-05              1        0.00019  
    50          53      1.65046e-05              1       0.000322  
    51          54      1.62752e-05              1       0.000507  
    52          55      1.57471e-05              1       0.000742  
    53          56      1.46874e-05              1       0.000947  
    54          57      1.31357e-05              1        0.00091  
    55          58      1.18881e-05              1       0.000508  
    56          59      1.14772e-05              1       0.000136  
    57          60      1.14286e-05              1       5.23e-05  
    58          61      1.14229e-05              1       5.04e-05  
    59          62      1.14151e-05              1       4.84e-05  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    60          63      1.13952e-05              1       6.05e-05  
    61          64      1.13518e-05              1       0.000106  
    62          65      1.12658e-05              1       0.000151  
    63          66      1.11459e-05              1       0.000154  
    64          67      1.10541e-05              1       9.73e-05  
    65          68      1.10244e-05              1       4.63e-05  
    66          69        1.102e-05              1       1.62e-05  
    67          70      1.10188e-05              1       1.62e-05  
    68          71      1.10168e-05              1        1.9e-05  
    69          72      1.10137e-05              1       4.32e-05  
    70          73      1.10075e-05              1       7.22e-05  
    71          74      1.09951e-05              1       0.000101  
    72          75        1.097e-05              1       0.000122  
    73          76      1.09221e-05              1       0.000171  
    74          77      1.08379e-05              1       0.000259  
    75          78       1.0699e-05              1       0.000347  
    76          79       1.0471e-05              1       0.000408  
    77          80      1.00948e-05              1       0.000441  
    78          81      9.53725e-06              1       0.000567  
    79          82      8.95372e-06              1       0.000476  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    80          83      8.63267e-06              1       0.000207  
    81          84      8.55637e-06              1       7.67e-05  
    82          85      8.54499e-06              1       7.01e-05  
    83          86      8.53869e-06              1       6.28e-05  
    84          87      8.52676e-06              1       4.53e-05  
    85          88      8.51202e-06              1       3.47e-05  
    86          89      8.50152e-06              1       2.37e-05  
    87          90      8.49846e-06              1       1.73e-05  
    88          91      8.49774e-06              1       1.73e-05  
    89          92      8.49697e-06              1       1.72e-05  
    90          93      8.49469e-06              1       1.74e-05  
    91          94      8.48908e-06              1       3.22e-05  
    92          95      8.47442e-06              1       5.55e-05  
    93          96      8.43847e-06              1       8.91e-05  
    94          97      8.35693e-06              1       0.000131  
    95          98      8.20753e-06              1       0.000159  
    96          99      8.03032e-06              1       0.000132  
    97         100       7.9304e-06              1       5.63e-05  
    98         101      7.90958e-06              1        2.2e-05  
    99         102      7.90773e-06              1          2e-05  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   100         103      7.90712e-06              1       1.82e-05  
   101         104      7.90555e-06              1        1.2e-05  
   102         105      7.90385e-06              1       7.59e-06  
   103         106       7.9027e-06              1       7.77e-06  
   104         107      7.90236e-06              1       7.97e-06  
   105         108      7.90222e-06              1       8.07e-06  
   106         109      7.90196e-06              1       8.17e-06  
   107         110      7.90127e-06              1       8.29e-06  
   108         111      7.89952e-06              1       9.73e-06  
   109         112      7.89519e-06              1       1.76e-05  
   110         113      7.88561e-06              1       2.51e-05  
   111         114      7.86882e-06              1       3.08e-05  
   112         115      7.85058e-06              1       2.35e-05  
   113         116      7.84158e-06              1       9.47e-06  
   114         117      7.83999e-06              1       5.42e-06  
   115         118      7.83987e-06              1       5.17e-06  
   116         119      7.83982e-06              1       4.99e-06  
   117         120      7.83965e-06              1       4.29e-06  
   118         121      7.83931e-06              1       3.37e-06  
   119         122      7.83869e-06              1       3.37e-06  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   120         123      7.83805e-06              1       3.37e-06  
   121         124      7.83771e-06              1       4.51e-06  
   122         125      7.83762e-06              1        4.4e-06  
   123         126      7.83758e-06              1       4.07e-06  
   124         127      7.83746e-06              1       3.47e-06  
   125         128      7.83719e-06              1       3.37e-06  
   126         129      7.83645e-06              1       3.37e-06  
   127         130      7.83457e-06              1       3.98e-06  
   128         131       7.8298e-06              1       6.44e-06  
   129         132      7.81839e-06              1       1.12e-05  
   130         133      7.79446e-06              1       1.75e-05  
   131         134      7.75733e-06              1          2e-05  
   132         135      7.72492e-06              1       1.36e-05  
   133         136      7.71317e-06              1       8.83e-06  
   134         137      7.71167e-06              1       8.44e-06  
   135         138      7.71152e-06              1       8.27e-06  
   136         139      7.71132e-06              1       8.11e-06  
   137         140      7.71072e-06              1       7.81e-06  
   138         141      7.70926e-06              1       7.34e-06  
   139         142      7.70555e-06              1       6.59e-06  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   140         143      7.69707e-06              1       5.53e-06  
   141         144      7.68078e-06              1       1.05e-05  
   142         145      7.65963e-06              1       2.03e-05  
   143         146      7.64555e-06              1       2.64e-05  
   144         147      7.64151e-06              1       2.64e-05  
   145         148      7.64051e-06              1       2.51e-05  
   146         149      7.63923e-06              1       2.37e-05  
   147         150      7.63571e-06              1       2.11e-05  
   148         151      7.62705e-06              1       1.69e-05  
   149         152      7.60557e-06              1       9.75e-06  
   150         153      7.55836e-06              1       9.27e-06  
   151         154      7.47608e-06              1       1.25e-05  
   152         155      7.38809e-06              1       1.68e-05  
   153         156      7.34567e-06              1       1.01e-05  
   154         157      7.33841e-06              1       7.72e-06  
   155         158       7.3379e-06              1       7.57e-06  
   156         159      7.33772e-06              1       7.53e-06  
   157         160      7.33698e-06              1       7.41e-06  
   158         161      7.33535e-06              1       7.26e-06  
   159         162      7.33093e-06              1          7e-06  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   160         163      7.32059e-06              1       7.27e-06  
   161         164      7.29852e-06              1       1.23e-05  
   162         165      7.26363e-06              1       1.78e-05  
   163         166      7.23151e-06              1       2.98e-05  
   164         167      7.21834e-06              1       3.43e-05  
   165         168      7.21556e-06              1       3.35e-05  
   166         169      7.21406e-06              1       3.24e-05  
   167         170      7.20993e-06              1       3.02e-05  
   168         171      7.20001e-06              1       2.67e-05  
   169         172      7.17407e-06              1       2.07e-05  
   170         173      7.11155e-06              1       1.36e-05  
   171         174      6.97464e-06              1       2.34e-05  
   172         175      6.74337e-06              1       2.97e-05  
   173         176      6.50775e-06              1       2.35e-05  
   174         177      6.40236e-06              1       1.44e-05  
   175         178      6.38595e-06              1       1.45e-05  
   176         179      6.38485e-06              1       1.42e-05  
   177         180       6.3844e-06              1       1.41e-05  
   178         181      6.38252e-06              1       1.38e-05  
   179         182      6.37836e-06              1       1.33e-05  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   180         183      6.36686e-06              1       1.24e-05  
   181         184      6.33832e-06              1       1.08e-05  
   182         185       6.2683e-06              1       1.54e-05  
   183         186      6.11691e-06              1       2.25e-05  
   184         187        5.863e-06              1        2.6e-05  
   185         188      5.60958e-06              1        1.9e-05  
   186         189      5.49985e-06              1       1.08e-05  
   187         190      5.48353e-06              1       1.05e-05  
   188         191      5.48259e-06              1       1.05e-05  
   189         192      5.48235e-06              1       1.05e-05  
   190         193      5.48129e-06              1       1.05e-05  
   191         194      5.47899e-06              1       1.04e-05  
   192         195      5.47255e-06              1       1.03e-05  
   193         196      5.45642e-06              1       1.02e-05  
   194         197      5.41554e-06              1       9.89e-06  
   195         198      5.32018e-06              1       9.46e-06  
   196         199      5.12999e-06              1       1.17e-05  
   197         200      4.86498e-06              1       2.01e-05  
   198         201      4.67353e-06              1       2.27e-05  
   199         202       4.6202e-06              1       1.97e-05  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   200         203      4.61475e-06              1       1.77e-05  
   201         204      4.61399e-06              1       1.72e-05  
   202         205      4.61233e-06              1       1.65e-05  
   203         206      4.60809e-06              1       1.55e-05  
   204         207      4.59694e-06              1       1.37e-05  
   205         208      4.56804e-06              1       1.15e-05  
   206         209       4.4938e-06              1       1.03e-05  
   207         210      4.30933e-06              1       9.37e-06  
   208         211      3.88643e-06              1       1.11e-05  
   209         212      3.08978e-06              1       2.15e-05  
   210         213      2.10838e-06              1       2.36e-05  
   211         214       1.5248e-06              1       1.34e-05  
   212         215      1.38949e-06              1       6.31e-06  
   213         216      1.37965e-06              1       6.26e-06  
   214         217       1.3794e-06              1       6.22e-06  
   215         218      1.37933e-06              1       6.21e-06  
   216         219       1.3789e-06              1       6.16e-06  
   217         220        1.378e-06              1       6.08e-06  
   218         221      1.37547e-06              1       5.91e-06  
   219         222      1.36925e-06              1       5.55e-06  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   220         223      1.35396e-06              1       4.78e-06  
   221         224      1.32115e-06              1       3.26e-06  
   222         225      1.26685e-06              1       1.84e-06  
   223         226      1.21403e-06              1       2.75e-06  
   224         227      1.19185e-06              1       3.13e-06  
   225         228      1.18868e-06              1       2.96e-06  
   226         229      1.18853e-06              1       2.87e-06  
   227         230       1.1885e-06              1       2.85e-06  
   228         231      1.18839e-06              1       2.79e-06  
   229         232      1.18816e-06              1        2.7e-06  
   230         233       1.1875e-06              1        2.5e-06  
   231         234      1.18587e-06              1       2.12e-06  
   232         235      1.18187e-06              1       1.92e-06  
   233         236       1.1733e-06              1        1.9e-06  
   234         237      1.15913e-06              1       2.49e-06  
   235         238      1.14536e-06              1       4.48e-06  
   236         239      1.13956e-06              1       4.82e-06  
   237         240      1.13871e-06              1       4.38e-06  
   238         241      1.13864e-06              1       4.19e-06  
   239         242      1.13858e-06              1       4.09e-06  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   240         243      1.13838e-06              1       3.86e-06  
   241         244      1.13793e-06              1       3.52e-06  
   242         245      1.13668e-06              1       2.94e-06  
   243         246      1.13352e-06              1       2.16e-06  
   244         247      1.12544e-06              1        2.2e-06  
   245         248      1.10606e-06              1       2.17e-06  
   246         249       1.0649e-06              1       4.08e-06  
   247         250      9.99402e-07              1       6.71e-06  
   248         251      9.39642e-07              1       8.27e-06  
   249         252      9.16764e-07              1       7.56e-06  
   250         253      9.13808e-07              1       6.65e-06  
   251         254      9.13648e-07              1       6.42e-06  
   252         255      9.13585e-07              1       6.36e-06  
   253         256      9.13304e-07              1       6.15e-06  
   254         257      9.12699e-07              1       5.84e-06  
   255         258      9.11041e-07              1        5.2e-06  
   256         259      9.07138e-07              1       4.03e-06  
   257         260      8.98647e-07              1       1.99e-06  
   258         261      8.84769e-07              1       2.01e-06  
   259         262      8.71369e-07              1       3.03e-06  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   260         263      8.65823e-07              1       2.98e-06  
   261         264      8.65032e-07              1       2.83e-06  
   262         265      8.64982e-07              1       2.84e-06  
   263         266       8.6496e-07              1       2.84e-06  
   264         267      8.64865e-07              1       2.84e-06  
   265         268      8.64654e-07              1       2.84e-06  
   266         269      8.64066e-07              1       2.83e-06  
   267         270      8.62578e-07              1       2.81e-06  
   268         271      8.58747e-07              1       2.75e-06  
   269         272      8.49396e-07              1       3.06e-06  
   270         273      8.28667e-07              1       5.53e-06  
   271         274      7.92473e-07              1        7.4e-06  
   272         275       7.5338e-07              1       6.42e-06  
   273         276      7.34302e-07              1       3.01e-06  
   274         277      7.31039e-07              1       3.28e-06  
   275         278      7.30851e-07              1       3.26e-06  
   276         279      7.30826e-07              1       3.24e-06  
   277         280      7.30741e-07              1        3.2e-06  
   278         281      7.30549e-07              1       3.12e-06  
   279         282      7.30027e-07              1       2.97e-06  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   280         283      7.28744e-07              1       2.64e-06  
   281         284      7.25711e-07              1       1.98e-06  
   282         285       7.1968e-07              1       1.32e-06  
   283         286      7.11231e-07              1       1.14e-06  
   284         287      7.05054e-07              1       1.73e-06  
   285         288      7.03295e-07              1       1.65e-06  
   286         289      7.03129e-07              1       1.58e-06  
   287         290      7.03118e-07              1       1.59e-06  
   288         291      7.03104e-07              1        1.6e-06  
   289         292       7.0306e-07              1       1.61e-06  
   290         293      7.02954e-07              1       1.62e-06  
   291         294       7.0267e-07              1       1.63e-06  
   292         295      7.01954e-07              1       1.62e-06  
   293         296      7.00189e-07              1       1.51e-06  
   294         297      6.96278e-07              1       1.17e-06  
   295         298      6.89393e-07              1       1.92e-06  
   296         299      6.81862e-07              1       1.96e-06  
   297         300      6.78113e-07              1        1.3e-06  
   298         301      6.77458e-07              1       1.35e-06  
   299         302      6.77419e-07              1       1.29e-06  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   300         303      6.77415e-07              1       1.27e-06  
   301         304      6.77399e-07              1       1.24e-06  
   302         305      6.77364e-07              1       1.19e-06  
   303         306      6.77268e-07              1        1.1e-06  
   304         307      6.77021e-07              1       9.62e-07  
   305         308      6.76377e-07              1       7.76e-07  
   306         309      6.74741e-07              1       7.64e-07  
   307         310      6.70753e-07              1       7.58e-07  
   308         311      6.62046e-07              1       1.35e-06  
   309         312      6.47278e-07              1       1.93e-06  
   310         313      6.32184e-07              1          2e-06  
   311         314      6.25388e-07              1       2.22e-06  
   312         315      6.24335e-07              1       2.27e-06  
   313         316      6.24281e-07              1       2.23e-06  
   314         317      6.24274e-07              1       2.22e-06  
   315         318      6.24248e-07              1       2.18e-06  
   316         319      6.24189e-07              1       2.13e-06  
   317         320       6.2403e-07              1       2.02e-06  
   318         321      6.23641e-07              1        1.8e-06  
   319         322      6.22728e-07              1       1.36e-06  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   320         323      6.20945e-07              1       5.85e-07  
   321         324      6.18534e-07              1       5.83e-07  
   322         325      6.16873e-07              1       8.79e-07  
   323         326      6.16435e-07              1       8.35e-07  
   324         327      6.16397e-07              1       8.56e-07  
   325         328      6.16395e-07              1       8.56e-07  
   326         329      6.16391e-07              1       8.56e-07  
   327         330       6.1638e-07              1       8.55e-07  
   328         331      6.16351e-07              1       8.53e-07  
   329         332      6.16276e-07              1        8.5e-07  
   330         333       6.1608e-07              1       8.43e-07  
   331         334      6.15572e-07              1       8.28e-07  
   332         335      6.14275e-07              1       8.69e-07  
   333         336      6.11089e-07              1       1.72e-06  
   334         337      6.03992e-07              1       2.79e-06  
   335         338      5.91417e-07              1       3.54e-06  
   336         339      5.77481e-07              1       2.97e-06  
   337         340       5.7041e-07              1       1.55e-06  
   338         341      5.69146e-07              1       1.58e-06  
   339         342      5.69074e-07              1       1.55e-06  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   340         343       5.6907e-07              1       1.54e-06  
   341         344      5.69061e-07              1       1.53e-06  
   342         345      5.69038e-07              1       1.51e-06  
   343         346      5.68979e-07              1       1.46e-06  
   344         347      5.68828e-07              1       1.37e-06  
   345         348      5.68453e-07              1       1.18e-06  
   346         349      5.67605e-07              1       8.04e-07  
   347         350      5.66044e-07              1       3.93e-07  
   348         351      5.64185e-07              1       5.37e-07  
   349         352      5.63133e-07              1       8.21e-07  
   350         353      5.62917e-07              1       7.16e-07  
   351         354      5.62903e-07              1       6.32e-07  
   352         355      5.62901e-07              1       6.15e-07  
   353         356      5.62897e-07              1       5.87e-07  
   354         357      5.62887e-07              1       5.52e-07  
   355         358      5.62861e-07              1       5.51e-07  
   356         359      5.62794e-07              1       5.49e-07  
   357         360       5.6262e-07              1       5.45e-07  
   358         361      5.62178e-07              1       5.37e-07  
   359         362      5.61101e-07              1       5.87e-07  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   360         363      5.58757e-07              1       1.18e-06  
   361         364       5.5481e-07              1       1.64e-06  
   362         365      5.50828e-07              1       1.45e-06  
   363         366       5.4907e-07              1       7.44e-07  
   364         367      5.48803e-07              1       7.32e-07  
   365         368       5.4879e-07              1       7.34e-07  
   366         369      5.48789e-07              1       7.31e-07  
   367         370      5.48783e-07              1       7.23e-07  
   368         371      5.48772e-07              1       7.11e-07  
   369         372      5.48739e-07              1       6.85e-07  
   370         373      5.48657e-07              1       6.32e-07  
   371         374      5.48448e-07              1       5.19e-07  
   372         375      5.47954e-07              1       4.13e-07  
   373         376      5.46935e-07              1       3.99e-07  
   374         377      5.45412e-07              1       8.14e-07  
   375         378      5.44167e-07              1       1.29e-06  
   376         379      5.43757e-07              1       1.33e-06  
   377         380      5.43711e-07              1       1.23e-06  
   378         381      5.43707e-07              1       1.19e-06  
   379         382      5.43704e-07              1       1.17e-06  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   380         383      5.43692e-07              1       1.12e-06  
   381         384      5.43665e-07              1       1.05e-06  
   382         385      5.43591e-07              1       9.32e-07  
   383         386      5.43401e-07              1       7.37e-07  
   384         387      5.42916e-07              1        6.1e-07  
   385         388      5.41737e-07              1       5.85e-07  
   386         389      5.39154e-07              1       7.42e-07  
   387         390      5.34761e-07              1       1.33e-06  
   388         391      5.30245e-07              1       1.26e-06  
   389         392      5.28196e-07              1       5.75e-07  
   390         393      5.27876e-07              1       3.22e-07  
   391         394      5.27861e-07              1       3.17e-07  
   392         395      5.27861e-07              1       3.16e-07  
   393         396       5.2786e-07              1       3.15e-07  
   394         397      5.27858e-07              1       3.12e-07  
   395         398      5.27853e-07              1       3.07e-07  
   396         399      5.27839e-07              1       2.99e-07  
   397         400      5.27803e-07              1       2.81e-07  
   398         401      5.27711e-07              1       2.56e-07  
   399         402      5.27482e-07              1       2.52e-07  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   400         403      5.26947e-07              1       2.46e-07  
   401         404      5.25902e-07              1       3.35e-07  
   402         405      5.24498e-07              1       7.59e-07  
   403         406      5.23536e-07              1       1.02e-06  
   404         407      5.23286e-07              1          1e-06  
   405         408      5.23264e-07              1       9.48e-07  
   406         409      5.23263e-07              1       9.34e-07  
   407         410      5.23261e-07              1       9.22e-07  
   408         411      5.23255e-07              1       8.99e-07  
   409         412      5.23241e-07              1       8.65e-07  
   410         413      5.23202e-07              1       8.07e-07  
   411         414      5.23102e-07              1       7.13e-07  
   412         415      5.22842e-07              1        5.6e-07  
   413         416       5.2218e-07              1       3.12e-07  
   414         417      5.20561e-07              1       2.37e-07  
   415         418      5.16998e-07              1       5.84e-07  
   416         419      5.10846e-07              1       1.03e-06  
   417         420      5.04351e-07              1       9.58e-07  
   418         421      5.01286e-07              1       5.82e-07  
   419         422      5.00786e-07              1       5.15e-07  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   420         423      5.00761e-07              1       4.95e-07  
   421         424       5.0076e-07              1       4.93e-07  
   422         425      5.00759e-07              1       4.91e-07  
   423         426      5.00756e-07              1       4.88e-07  
   424         427      5.00748e-07              1       4.83e-07  
   425         428      5.00726e-07              1       4.74e-07  
   426         429      5.00672e-07              1       4.59e-07  
   427         430      5.00536e-07              1       4.34e-07  
   428         431      5.00227e-07              1       3.98e-07  
   429         432      4.99653e-07              1       3.57e-07  
   430         433      4.98954e-07              1       3.38e-07  
   431         434      4.98544e-07              1        4.3e-07  
   432         435      4.98456e-07              1       3.99e-07  
   433         436       4.9845e-07              1       4.01e-07  
   434         437       4.9845e-07              1       4.02e-07  
   435         438      4.98449e-07              1       4.03e-07  
   436         439      4.98447e-07              1       4.04e-07  
   437         440      4.98441e-07              1       4.06e-07  
   438         441      4.98426e-07              1       4.08e-07  
   439         442      4.98387e-07              1        4.1e-07  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   440         443      4.98288e-07              1       4.08e-07  
   441         444      4.98046e-07              1       3.92e-07  
   442         445      4.97507e-07              1        3.4e-07  
   443         446      4.96558e-07              1       4.06e-07  
   444         447      4.95519e-07              1       4.13e-07  
   445         448      4.95003e-07              1        2.1e-07  
   446         449      4.94913e-07              1       7.55e-08  
   447         450      4.94908e-07              1       7.53e-08  
   448         451      4.94908e-07              1       7.53e-08  
   449         452      4.94908e-07              1       7.53e-08  
   450         453      4.94908e-07              1       7.52e-08  
   451         454      4.94907e-07              1       7.52e-08  
   452         455      4.94906e-07              1       7.51e-08  
   453         456      4.94904e-07              1       7.49e-08  
   454         457      4.94897e-07              1       7.45e-08  
   455         458      4.94882e-07              1        7.4e-08  
   456         459      4.94847e-07              1        7.3e-08  
   457         460      4.94788e-07              1       9.62e-08  
   458         461      4.94726e-07              1       1.81e-07  
   459         462      4.94696e-07              1       2.13e-07  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   460         463      4.94691e-07              1       2.04e-07  
   461         464      4.94691e-07              1       1.97e-07  
   462         465      4.94691e-07              1       1.95e-07  
   463         466      4.94691e-07              1        1.9e-07  
   464         467       4.9469e-07              1       1.84e-07  
   465         468      4.94689e-07              1       1.72e-07  
   466         469      4.94684e-07              1       1.54e-07  
   467         470      4.94674e-07              1       1.23e-07  
   468         471      4.94647e-07              1       7.27e-08  
   469         472      4.94581e-07              1       7.21e-08  
   470         473      4.94436e-07              1       1.22e-07  
   471         474      4.94191e-07              1       2.32e-07  
   472         475      4.93939e-07              1       2.44e-07  
   473         476      4.93825e-07              1       1.47e-07  
   474         477      4.93808e-07              1       1.19e-07  
   475         478      4.93807e-07              1       1.22e-07  
   476         479      4.93807e-07              1       1.22e-07  
   477         480      4.93807e-07              1       1.22e-07  
   478         481      4.93807e-07              1       1.21e-07  
   479         482      4.93806e-07              1        1.2e-07  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   480         483      4.93805e-07              1       1.17e-07  
   481         484      4.93802e-07              1        1.1e-07  
   482         485      4.93796e-07              1       9.58e-08  
   483         486       4.9378e-07              1       6.57e-08  
   484         487      4.93749e-07              1       5.47e-08  
   485         488       4.9371e-07              1       5.94e-08  
   486         489      4.93685e-07              1       9.85e-08  
   487         490      4.93679e-07              1       9.45e-08  
   488         491      4.93678e-07              1       8.61e-08  
   489         492      4.93678e-07              1       8.42e-08  
   490         493      4.93678e-07              1       8.28e-08  
   491         494      4.93678e-07              1       7.99e-08  
   492         495      4.93678e-07              1       7.56e-08  
   493         496      4.93677e-07              1       6.84e-08  
   494         497      4.93676e-07              1       5.68e-08  
   495         498      4.93671e-07              1       5.03e-08  
   496         499      4.93661e-07              1        4.9e-08  
   497         500      4.93633e-07              1       4.84e-08  
   498         501      4.93572e-07              1       1.05e-07  
   499         502      4.93456e-07              1       1.69e-07  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   500         503      4.93313e-07              1       1.73e-07  
   501         504      4.93228e-07              1       1.22e-07  
   502         505      4.93209e-07              1       1.26e-07  
   503         506      4.93208e-07              1       1.22e-07  
   504         507      4.93208e-07              1       1.21e-07  
   505         508      4.93208e-07              1        1.2e-07  
   506         509      4.93208e-07              1       1.19e-07  
   507         510      4.93208e-07              1       1.17e-07  
   508         511      4.93207e-07              1       1.13e-07  
   509         512      4.93205e-07              1       1.08e-07  
   510         513        4.932e-07              1       9.84e-08  
   511         514      4.93187e-07              1       8.31e-08  
   512         515      4.93154e-07              1       5.79e-08  
   513         516      4.93074e-07              1       3.55e-08  
   514         517      4.92905e-07              1       3.65e-08  
   515         518      4.92633e-07              1       8.47e-08  
   516         519       4.9238e-07              1       8.22e-08  
   517         520      4.92282e-07              1       3.61e-08  
   518         521      4.92269e-07              1       3.01e-08  
   519         522      4.92269e-07              1       2.93e-08  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   520         523      4.92269e-07              1       2.92e-08  
   521         524      4.92269e-07              1       2.92e-08  
   522         525      4.92269e-07              1       2.91e-08  
   523         526      4.92269e-07              1       2.89e-08  
   524         527      4.92268e-07              1       2.85e-08  
   525         528      4.92268e-07              1       2.78e-08  
   526         529      4.92267e-07              1       2.62e-08  
   527         530      4.92263e-07              1       2.43e-08  
   528         531      4.92255e-07              1       2.44e-08  
   529         532      4.92237e-07              1       2.44e-08  
   530         533      4.92202e-07              1       3.93e-08  
   531         534      4.92159e-07              1       8.36e-08  
   532         535      4.92132e-07              1       1.08e-07  
   533         536      4.92126e-07              1       1.07e-07  
   534         537      4.92126e-07              1       1.02e-07  
   535         538      4.92126e-07              1       1.01e-07  
   536         539      4.92126e-07              1          1e-07  
   537         540      4.92126e-07              1       9.84e-08  
   538         541      4.92125e-07              1        9.6e-08  
   539         542      4.92125e-07              1       9.19e-08  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   540         543      4.92123e-07              1       8.52e-08  
   541         544      4.92119e-07              1       7.42e-08  
   542         545       4.9211e-07              1       5.61e-08  
   543         546      4.92084e-07              1       2.83e-08  
   544         547      4.92023e-07              1       2.83e-08  
   545         548      4.91891e-07              1       8.29e-08  
   546         549      4.91676e-07              1       1.38e-07  
   547         550       4.9147e-07              1       1.32e-07  
   548         551      4.91386e-07              1       1.21e-07  
   549         552      4.91375e-07              1       1.28e-07  
   550         553      4.91374e-07              1       1.27e-07  
   551         554      4.91374e-07              1       1.27e-07  
   552         555      4.91374e-07              1       1.26e-07  
   553         556      4.91374e-07              1       1.25e-07  
   554         557      4.91374e-07              1       1.23e-07  
   555         558      4.91372e-07              1       1.21e-07  
   556         559      4.91369e-07              1       1.16e-07  
   557         560      4.91361e-07              1       1.07e-07  
   558         561      4.91341e-07              1       9.06e-08  
   559         562      4.91293e-07              1       5.98e-08  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   560         563      4.91191e-07              1       3.74e-08  
   561         564      4.91032e-07              1       6.05e-08  
   562         565      4.90894e-07              1       1.02e-07  
   563         566      4.90844e-07              1       9.46e-08  
   564         567      4.90839e-07              1       7.86e-08  
   565         568      4.90838e-07              1       7.45e-08  
   566         569      4.90838e-07              1        7.4e-08  
   567         570      4.90838e-07              1       7.29e-08  
   568         571      4.90838e-07              1       7.15e-08  
   569         572      4.90838e-07              1       7.01e-08  
   570         573      4.90836e-07              1       6.93e-08  
   571         574      4.90833e-07              1        6.8e-08  
   572         575      4.90825e-07              1       6.57e-08  
   573         576      4.90804e-07              1       6.21e-08  
   574         577      4.90752e-07              1       5.61e-08  
   575         578      4.90628e-07              1       4.72e-08  
   576         579       4.9038e-07              1       1.21e-07  
   577         580      4.90031e-07              1       2.12e-07  
   578         581      4.89771e-07              1       2.61e-07  
   579         582      4.89696e-07              1       2.49e-07  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   580         583      4.89688e-07              1       2.33e-07  
   581         584      4.89688e-07              1       2.29e-07  
   582         585      4.89688e-07              1       2.29e-07  
   583         586      4.89688e-07              1       2.25e-07  
   584         587      4.89686e-07              1       2.21e-07  
   585         588      4.89683e-07              1       2.14e-07  
   586         589      4.89675e-07              1       2.02e-07  
   587         590      4.89655e-07              1        1.8e-07  
   588         591      4.89603e-07              1       1.44e-07  
   589         592      4.89476e-07              1       8.11e-08  
   590         593      4.89207e-07              1        5.6e-08  
   591         594      4.88772e-07              1       1.24e-07  
   592         595      4.88366e-07              1       1.68e-07  
   593         596      4.88205e-07              1       1.82e-07  
   594         597      4.88184e-07              1       1.93e-07  
   595         598      4.88183e-07              1       1.93e-07  
   596         599      4.88183e-07              1       1.92e-07  
   597         600      4.88183e-07              1       1.92e-07  
   598         601      4.88183e-07              1       1.91e-07  
   599         602      4.88181e-07              1       1.89e-07  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   600         603      4.88179e-07              1       1.86e-07  
   601         604      4.88171e-07              1       1.81e-07  
   602         605      4.88152e-07              1       1.72e-07  
   603         606      4.88102e-07              1       1.57e-07  
   604         607       4.8798e-07              1       1.29e-07  
   605         608        4.877e-07              1       8.03e-08  
   606         609      4.87173e-07              1       1.08e-07  
   607         610      4.86521e-07              1       1.15e-07  
   608         611      4.86127e-07              1       7.62e-08  
   609         612      4.86041e-07              1       7.64e-08  
   610         613      4.86035e-07              1       7.33e-08  
   611         614      4.86035e-07              1       7.26e-08  
   612         615      4.86035e-07              1       7.25e-08  
   613         616      4.86035e-07              1       7.22e-08  
   614         617      4.86034e-07              1       7.18e-08  
   615         618      4.86034e-07              1       7.09e-08  
   616         619      4.86032e-07              1       6.93e-08  
   617         620      4.86029e-07              1       6.63e-08  
   618         621      4.86019e-07              1       6.03e-08  
   619         622      4.85996e-07              1        4.8e-08  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   620         623      4.85941e-07              1       3.98e-08  
   621         624      4.85832e-07              1       3.53e-08  
   622         625      4.85685e-07              1       7.77e-08  
   623         626      4.85581e-07              1        1.1e-07  
   624         627      4.85554e-07              1       1.05e-07  
   625         628      4.85551e-07              1       9.69e-08  
   626         629      4.85551e-07              1       9.51e-08  
   627         630      4.85551e-07              1       9.48e-08  
   628         631      4.85551e-07              1       9.37e-08  
   629         632      4.85551e-07              1       9.23e-08  
   630         633       4.8555e-07              1       8.97e-08  
   631         634      4.85549e-07              1       8.56e-08  
   632         635      4.85546e-07              1       7.86e-08  
   633         636      4.85537e-07              1       6.69e-08  
   634         637      4.85515e-07              1       4.75e-08  
   635         638      4.85463e-07              1        3.3e-08  
   636         639      4.85356e-07              1       2.59e-08  
   637         640      4.85198e-07              1       6.14e-08  
   638         641      4.85072e-07              1       5.78e-08  
   639         642      4.85032e-07              1       2.79e-08  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   640         643      4.85028e-07              1       2.86e-08  
   641         644      4.85028e-07              1       2.84e-08  
   642         645      4.85028e-07              1       2.84e-08  
   643         646      4.85027e-07              1       2.83e-08  
   644         647      4.85027e-07              1       2.82e-08  
   645         648      4.85027e-07              1       2.81e-08  
   646         649      4.85027e-07              1       2.78e-08  
   647         650      4.85027e-07              1       2.74e-08  
   648         651      4.85026e-07              1       2.66e-08  
   649         652      4.85022e-07              1       2.53e-08  
   650         653      4.85014e-07              1       2.28e-08  
   651         654      4.84994e-07              1       2.26e-08  
   652         655      4.84945e-07              1       2.42e-08  
   653         656      4.84842e-07              1       3.36e-08  
   654         657      4.84683e-07              1       3.59e-08  
   655         658      4.84545e-07              1        4.3e-08  
   656         659      4.84496e-07              1       4.15e-08  
   657         660      4.84491e-07              1       3.69e-08  
   658         661       4.8449e-07              1       3.57e-08  
   659         662       4.8449e-07              1       3.56e-08  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   660         663       4.8449e-07              1       3.55e-08  
   661         664       4.8449e-07              1       3.53e-08  
   662         665       4.8449e-07              1       3.51e-08  
   663         666       4.8449e-07              1       3.46e-08  
   664         667       4.8449e-07              1        3.4e-08  
   665         668      4.84488e-07              1       3.28e-08  
   666         669      4.84485e-07              1        3.1e-08  
   667         670      4.84478e-07              1        2.8e-08  
   668         671      4.84458e-07              1        2.3e-08  
   669         672      4.84409e-07              1       2.15e-08  
   670         673      4.84302e-07              1       3.13e-08  
   671         674      4.84116e-07              1       3.69e-08  
   672         675      4.83919e-07              1       3.75e-08  
   673         676      4.83825e-07              1       3.56e-08  
   674         677      4.83809e-07              1       2.88e-08  
   675         678      4.83808e-07              1       2.66e-08  
   676         679      4.83808e-07              1       2.64e-08  
   677         680      4.83808e-07              1       2.64e-08  
   678         681      4.83808e-07              1       2.63e-08  
   679         682      4.83808e-07              1       2.62e-08  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   680         683      4.83808e-07              1        2.6e-08  
   681         684      4.83808e-07              1       2.57e-08  
   682         685      4.83807e-07              1       2.52e-08  
   683         686      4.83806e-07              1       2.41e-08  
   684         687      4.83802e-07              1       2.21e-08  
   685         688      4.83792e-07              1       1.85e-08  
   686         689       4.8377e-07              1       2.56e-08  
   687         690       4.8373e-07              1       3.23e-08  
   688         691       4.8368e-07              1       3.93e-08  
   689         692      4.83652e-07              1       4.93e-08  
   690         693      4.83646e-07              1       4.82e-08  
   691         694      4.83645e-07              1       4.63e-08  
   692         695      4.83645e-07              1       4.59e-08  
   693         696      4.83645e-07              1       4.58e-08  
   694         697      4.83645e-07              1       4.55e-08  
   695         698      4.83645e-07              1        4.5e-08  
   696         699      4.83645e-07              1       4.42e-08  
   697         700      4.83645e-07              1       4.28e-08  
   698         701      4.83644e-07              1       4.04e-08  
   699         702      4.83642e-07              1       3.61e-08  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   700         703      4.83636e-07              1       2.79e-08  
   701         704      4.83622e-07              1       2.36e-08  
   702         705      4.83594e-07              1       3.34e-08  
   703         706      4.83548e-07              1       4.63e-08  
   704         707      4.83507e-07              1       6.75e-08  
   705         708      4.83491e-07              1       6.38e-08  
   706         709      4.83489e-07              1       5.54e-08  
   707         710      4.83489e-07              1       5.32e-08  
   708         711      4.83489e-07              1       5.29e-08  
   709         712      4.83489e-07              1       5.24e-08  
   710         713      4.83489e-07              1       5.16e-08  
   711         714      4.83489e-07              1       5.04e-08  
   712         715      4.83489e-07              1       4.84e-08  
   713         716      4.83488e-07              1       4.51e-08  
   714         717      4.83485e-07              1       3.99e-08  
   715         718      4.83479e-07              1       3.13e-08  
   716         719      4.83462e-07              1       2.43e-08  
   717         720       4.8342e-07              1       3.95e-08  
   718         721      4.83311e-07              1       6.44e-08  
   719         722      4.83052e-07              1       9.94e-08  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   720         723      4.82509e-07              1       1.46e-07  
   721         724       4.8167e-07              1       1.71e-07  
   722         725      4.80944e-07              1       1.18e-07  
   723         726      4.80685e-07              1       3.81e-08  
   724         727      4.80655e-07              1       3.82e-08  
   725         728      4.80654e-07              1       3.78e-08  
   726         729      4.80654e-07              1       3.78e-08  
   727         730      4.80654e-07              1       3.77e-08  
   728         731      4.80654e-07              1       3.77e-08  
   729         732      4.80654e-07              1       3.76e-08  
   730         733      4.80654e-07              1       3.74e-08  
   731         734      4.80653e-07              1       3.72e-08  
   732         735      4.80652e-07              1       3.67e-08  
   733         736      4.80649e-07              1       3.58e-08  
   734         737       4.8064e-07              1       3.43e-08  
   735         738      4.80618e-07              1       3.43e-08  
   736         739      4.80563e-07              1       5.65e-08  
   737         740      4.80434e-07              1       8.67e-08  
   738         741      4.80178e-07              1       1.15e-07  
   739         742      4.79822e-07              1       1.12e-07  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   740         743      4.79565e-07              1       7.79e-08  
   741         744      4.79493e-07              1       7.64e-08  
   742         745      4.79487e-07              1       7.11e-08  
   743         746      4.79487e-07              1       6.99e-08  
   744         747      4.79487e-07              1       6.98e-08  
   745         748      4.79487e-07              1       6.96e-08  
   746         749      4.79486e-07              1       6.92e-08  
   747         750      4.79486e-07              1       6.87e-08  
   748         751      4.79486e-07              1       6.78e-08  
   749         752      4.79484e-07              1       6.64e-08  
   750         753      4.79481e-07              1       6.41e-08  
   751         754      4.79472e-07              1       6.04e-08  
   752         755      4.79448e-07              1       5.45e-08  
   753         756      4.79387e-07              1       6.36e-08  
   754         757      4.79236e-07              1       1.03e-07  
   755         758      4.78891e-07              1       1.53e-07  
   756         759      4.78243e-07              1       1.93e-07  
   757         760      4.77442e-07              1       1.92e-07  
   758         761       4.7696e-07              1       2.39e-07  
   759         762      4.76854e-07              1       2.35e-07  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   760         763      4.76847e-07              1       2.26e-07  
   761         764      4.76847e-07              1       2.25e-07  
   762         765      4.76847e-07              1       2.24e-07  
   763         766      4.76846e-07              1       2.23e-07  
   764         767      4.76846e-07              1       2.21e-07  
   765         768      4.76843e-07              1       2.17e-07  
   766         769      4.76837e-07              1       2.11e-07  
   767         770      4.76821e-07              1       2.01e-07  
   768         771       4.7678e-07              1       1.84e-07  
   769         772      4.76673e-07              1       1.55e-07  
   770         773      4.76414e-07              1       1.06e-07  
   771         774      4.75836e-07              1       1.26e-07  
   772         775       4.7482e-07              1       1.51e-07  
   773         776       4.7371e-07              1       1.21e-07  
   774         777      4.73158e-07              1       8.84e-08  
   775         778      4.73062e-07              1       6.21e-08  
   776         779      4.73057e-07              1       6.27e-08  
   777         780      4.73057e-07              1       6.26e-08  
   778         781      4.73057e-07              1       6.25e-08  
   779         782      4.73057e-07              1       6.24e-08  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   780         783      4.73057e-07              1       6.22e-08  
   781         784      4.73056e-07              1       6.19e-08  
   782         785      4.73055e-07              1       6.13e-08  
   783         786      4.73052e-07              1       6.02e-08  
   784         787      4.73042e-07              1       5.82e-08  
   785         788      4.73019e-07              1       6.05e-08  
   786         789      4.72962e-07              1       9.83e-08  
   787         790      4.72837e-07              1       1.44e-07  
   788         791      4.72622e-07              1        1.7e-07  
   789         792      4.72397e-07              1       1.32e-07  
   790         793      4.72292e-07              1       5.07e-08  
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
the forcing frequency 3.7616e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.7708e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.7788e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.7857e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.7919e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.7977e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8046e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8112e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8165e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8208e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8244e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8275e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8302e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8326e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8348e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8368e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8386e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8403e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8419e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8434e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8448e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8462e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8475e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8488e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8501e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8513e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8525e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8537e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8546e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8551e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8553e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8553e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8556e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8560e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8569e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8581e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8593e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8605e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8618e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8631e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8644e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8656e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8665e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8671e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8675e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8679e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8682e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8684e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8687e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8689e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8691e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8692e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8694e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8696e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8697e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8699e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8701e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8702e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8704e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8705e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8706e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8708e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8709e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8710e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8711e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8712e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8713e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8714e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8715e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8716e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8718e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8720e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8722e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8725e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8728e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8732e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8737e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8742e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8748e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8756e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8764e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8774e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8785e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8799e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8814e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8832e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8854e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8879e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8910e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8947e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8993e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.9052e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.9130e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.9236e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.9389e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.9543e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.9695e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.9890e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 4.0167e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 4.0770e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 4.0866e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
Calculate FRC in physical domain at epsilon 7.000000e-01
the forcing frequency 3.7589e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.7592e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.7689e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.7796e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.7880e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.7949e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8007e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8056e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8098e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8135e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8168e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8197e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8224e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8249e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8272e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8293e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8313e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8331e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8349e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8366e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8382e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8398e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8413e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8428e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8442e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8456e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8470e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8483e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8497e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8510e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8523e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8536e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8546e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8551e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8553e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8553e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8556e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8561e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8570e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8583e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8596e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8609e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8622e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8635e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8648e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8662e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8675e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8688e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8702e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8716e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8729e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8743e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8757e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8772e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8786e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8801e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8816e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8831e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8846e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8861e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8877e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8893e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8910e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8926e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8941e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8949e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8954e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8958e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8961e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8963e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8964e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8966e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8967e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8968e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8969e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8970e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8971e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8972e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8973e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8973e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8974e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8974e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8974e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8974e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8974e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8974e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8974e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8973e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8971e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8969e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8965e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8960e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8953e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8947e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8940e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8933e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8926e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8919e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8913e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8906e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8900e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8894e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8888e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8882e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8876e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8871e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8866e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8861e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8856e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8852e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8848e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8845e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8841e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8839e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8836e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8834e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8832e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8831e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8830e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8830e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8830e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8830e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8830e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8831e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8832e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8834e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8836e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8839e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8844e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8848e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8854e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8861e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8869e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8878e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8889e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8901e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8915e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8931e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8950e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8971e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.8996e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.9025e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.9058e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.9098e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.9146e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.9203e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.9274e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.9363e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.9478e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.9631e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 3.9845e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 4.0165e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 4.0689e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
the forcing frequency 4.0866e+00 is nearly resonant with the eigenvalue -1.0000e-02 + i3.8553e+00
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

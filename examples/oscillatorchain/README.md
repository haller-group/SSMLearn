This is a preview of the livescript `oscillator4D.mlx` (see also the other two present in the folder).

# Identifying a 4D SSM for an oscillator chain

This is an example of how to reconstruct the slow 4D SSM of a mechanical system using synthetic measurements of the full state space. In this example, we consider a damped oscillator chain with an additional nonlinear spring attached to the leftmost mass. The measurements for this example are transients occurring from exact initial conditions on the slow 4D SSM.

![image_0.png
](README_images/image_0.png
)

```matlab:Code
clearvars
close all
```

# Example setup

The $N$-degree of freedom dynamical system is of the form

$$
{M\ddot{q} }+{C\dot{q} }+{Kq}+f(q,{\dot{q} })=0
$$

where $f=\mathcal{O}(|q|^2 ,|{\dot{q} }|^2 ,|q||{\dot{q} }|)$ represents the nonlinearities and $M$, $C$, and $K$ are the $n\times m$ mass, stiffness, and damping matrices, respectively.

We rewrite the system in first-order form as

$$
{\dot{y} }=Ay+G(y)=F(y)
$$

with

$y=\left\lbrack \begin{array}{c}
q\\
\dot{q} 
\end{array}\right\rbrack ,~~A=\left\lbrack \begin{array}{cc}
0 & I\\
-M^{-1} K & -M^{-1} C
\end{array}\right\rbrack ,~~G(y)=\left\lbrack \begin{array}{c}
0\\
-M^{-1} f(y)
\end{array}\right\rbrack$.

```matlab:Code
n         = 5;      % number of masses
mass      = 1;
mass1     = 1.5;
stiffness = 1;
dampingM   = 0.002;
dampingK   = 0.005;

vMass = ones(1,n);%fliplr(linspace(0.5,1,n)); %
M = mass*diag(vMass); M(1,1) = mass1;
K = stiffness*(2*eye(n) - diag(ones(1,n-1),-1) - diag(ones(1,n-1),1));
K(end,end) = 1;
C = dampingM*M + dampingK*K;

% nonlinearities: f1 = 0.33*q1dot^2 + 3*q1^3 + 0.7*q1^2*q1dot + 0.5*q1dot^3
F2 = sptensor(zeros(n, 2*n, 2*n)); F3 = sptensor(zeros(n, 2*n, 2*n, 2*n));
F2(1,n+1,n+1) = 0.33;    % q1dot^2
F3(1,1,1,1) = 3;         % q1^3
F3(1,1,1,n+1) = 0.7;     % q1^2*q1dot 
F3(1,n+1,n+1,n+1) = 0.5; % q1dot^3
fnl = {F2*mass1, F3*mass1};

[F, lambda, E] = functionFromTensors(M, C, K, fnl);
```

# Generation of Synthetic Data

Having set up the dynamics of the problem, we now move on to generate synthetic data, which will be used to fit a parametrisation of the manifold. We will divide the data into a training set, for model fitting, and a test set, for validation.

```matlab:Code
nTraj = 6;
indTest = [1];
indTrain = setdiff(1:nTraj, indTest);
```

We now set the initial conditions for the trajectories, which are picked exactly on the slow 4D SSM of the systems, computed using SSMTool.

```matlab:Code
ICRadius = 0.3;
IC = getSSMIC(M, C, K, fnl, nTraj, ICRadius, 4, 1);
```

```text:Output
Due to high-dimensionality, we compute only the first 5 eigenvalues with the smallest magnitude. These would also be used to compute the spectral quotients
Assuming a proportional damping hypthesis with symmetric matrices
modal damping ratio for 1 mode is 4.245669e-03
modal damping ratio for 2 mode is 3.236388e-03
modal damping ratio for 3 mode is 3.848172e-03
modal damping ratio for 4 mode is 4.639624e-03
modal damping ratio for 5 mode is 5.267721e-03

 The first 10 eigenvalues are given as 
   -0.0011996 +    0.28254i
   -0.0011996 -    0.28254i
   -0.0025406 +      0.785i
   -0.0025406 -      0.785i
   -0.0046493 +     1.2082i
   -0.0046493 -     1.2082i
   -0.0074555 +     1.6069i
   -0.0074555 -     1.6069i
   -0.0099883 +     1.8961i
   -0.0099883 -     1.8961i

(near) outer resonance detected for the following combination of master eigenvalues
     0     4     3     0
     7     0     0     1
     4     0     0     3
     0     7     1     0
     3     0     1     0
     3     0     2     1
     4     1     1     0
     3     0     3     2
     4     1     2     1
     5     2     1     0
     0     3     0     1
     0     3     1     2
     1     4     0     1
     0     3     2     3
     1     4     1     2
     2     5     0     1
     4     0     1     0
     4     0     2     1
     5     1     1     0
     0     4     0     1
     0     4     1     2
     1     5     0     1

These are in resonance with the follwing eigenvalues of the slave subspace
   -0.0046493 +     1.2082i
   -0.0046493 +     1.2082i
   -0.0046493 -     1.2082i
   -0.0046493 -     1.2082i
   -0.0074555 +     1.6069i
   -0.0074555 +     1.6069i
   -0.0074555 +     1.6069i
   -0.0074555 +     1.6069i
   -0.0074555 +     1.6069i
   -0.0074555 +     1.6069i
   -0.0074555 -     1.6069i
   -0.0074555 -     1.6069i
   -0.0074555 -     1.6069i
   -0.0074555 -     1.6069i
   -0.0074555 -     1.6069i
   -0.0074555 -     1.6069i
   -0.0099883 +     1.8961i
   -0.0099883 +     1.8961i
   -0.0099883 +     1.8961i
   -0.0099883 -     1.8961i
   -0.0099883 -     1.8961i
   -0.0099883 -     1.8961i

sigma_out = 8
(near) inner resonance detected for the following combination of master eigenvalues
     1     0     1     1
     2     1     0     0
     1     0     2     2
     2     1     1     1
     3     2     0     0
     1     0     3     3
     2     1     2     2
     3     2     1     1
     4     3     0     0
     0     1     1     1
     1     2     0     0
     0     1     2     2
     1     2     1     1
     2     3     0     0
     0     1     3     3
     1     2     2     2
     2     3     1     1
     3     4     0     0
     0     0     2     1
     1     1     1     0
     0     0     3     2
     1     1     2     1
     2     2     1     0
     0     0     4     3
     1     1     3     2
     2     2     2     1
     3     3     1     0
     0     0     1     2
     1     1     0     1
     0     0     2     3
     1     1     1     2
     2     2     0     1
     0     0     3     4
     1     1     2     3
     2     2     1     2
     3     3     0     1

These are in resonance with the follwing eigenvalues of the master subspace
   -0.0011996 +    0.28254i
   -0.0011996 +    0.28254i
   -0.0011996 +    0.28254i
   -0.0011996 +    0.28254i
   -0.0011996 +    0.28254i
   -0.0011996 +    0.28254i
   -0.0011996 +    0.28254i
   -0.0011996 +    0.28254i
   -0.0011996 +    0.28254i
   -0.0011996 -    0.28254i
   -0.0011996 -    0.28254i
   -0.0011996 -    0.28254i
   -0.0011996 -    0.28254i
   -0.0011996 -    0.28254i
   -0.0011996 -    0.28254i
   -0.0011996 -    0.28254i
   -0.0011996 -    0.28254i
   -0.0011996 -    0.28254i
   -0.0025406 +      0.785i
   -0.0025406 +      0.785i
   -0.0025406 +      0.785i
   -0.0025406 +      0.785i
   -0.0025406 +      0.785i
   -0.0025406 +      0.785i
   -0.0025406 +      0.785i
   -0.0025406 +      0.785i
   -0.0025406 +      0.785i
   -0.0025406 -      0.785i
   -0.0025406 -      0.785i
   -0.0025406 -      0.785i
   -0.0025406 -      0.785i
   -0.0025406 -      0.785i
   -0.0025406 -      0.785i
   -0.0025406 -      0.785i
   -0.0025406 -      0.785i
   -0.0025406 -      0.785i

sigma_in = 8
Due to (near) outer resonance, the exisitence of the manifold is questionable and the underlying computation may suffer.
Attempting manifold computation
Manifold computation time at order 2 = 00:00:00
Estimated memory usage at order  2 = 2.82E-02 MB
Manifold computation time at order 3 = 00:00:00
Estimated memory usage at order  3 = 5.44E-02 MB
Manifold computation time at order 4 = 00:00:00
Estimated memory usage at order  4 = 7.99E-02 MB
Manifold computation time at order 5 = 00:00:00
Estimated memory usage at order  5 = 1.70E-01 MB
Manifold computation time at order 6 = 00:00:00
Estimated memory usage at order  6 = 2.52E-01 MB
Manifold computation time at order 7 = 00:00:00
Estimated memory usage at order  7 = 4.81E-01 MB
Manifold computation time at order 8 = 00:00:01
Estimated memory usage at order  8 = 6.89E-01 MB
Manifold computation time at order 9 = 00:00:01
Estimated memory usage at order  9 = 1.18E+00 MB
```

The data is generated by time-integration of the right-hand side of the system. In this case we are observing the full state space, so the `observable` function is the identity.

```matlab:Code
observable = @(x) x;
slowTime = 2*pi/abs(lambda(1));
nPeriods = 140; 
nDataPerPeriod = 50*2;
endTime = nPeriods * slowTime;
nSamp = nDataPerPeriod * nPeriods + 1;
dt = endTime/(nSamp-1);
tic
yData = integrateTrajectories(F, endTime, IC, nSamp, observable);
```

```text:Output
simulating trajectory 1 of 6...
simulating trajectory 2 of 6...
simulating trajectory 3 of 6...
simulating trajectory 4 of 6...
simulating trajectory 5 of 6...
simulating trajectory 6 of 6...
```

```matlab:Code
toc
```

```text:Output
Elapsed time is 3.034720 seconds.
```

# Plot generated trajectory data

```matlab:Code
customFigure(); colororder(winter(nTraj));
for iTraj = 1:size(yData,1)
    plot3(yData{iTraj,2}(1,:), yData{iTraj,2}(n+1,:), yData{iTraj,2}(2*n,:),'Linewidth',1,'Displayname',['Trajectory ' num2str(iTraj)])
end
xlabel('$q_1 \, [$m$]$','Interpreter','latex'); 
ylabel('$\dot{q}_1 \, [$m/s$]$','Interpreter','latex'); 
zlabel(['$\dot{q}_{' num2str(n) '} \, [$m/s$]$'],'Interpreter','latex'); 
view(3)
```

![figure_0.png
](README_images/figure_0.png
)

```matlab:Code

customFigure(); colororder(winter(nTraj));
for iTraj = 1:size(yData,1)
    plot(yData{iTraj,1}, yData{iTraj,2}(n,:),'Linewidth',1,'Displayname',['Trajectory ' num2str(iTraj)])
end
xlim([0 endTime])
xlabel('$t \, [s]$','Interpreter','latex'); 
ylabel('$q_1 \, [$m$]$','Interpreter','latex'); 
legend('location','NE')
```

![figure_1.png
](README_images/figure_1.png
)

# Datadriven manifold fitting

The measured trajectories are assumed to lie very close to a four-dimensional manifold that is tangent at the origin to the eigenspace corresponding to the two slowest pair of eigenvalues. We now want to fit a polynomial of order $M$ to the data points to approximate the manifold. Here we use a graph style parametrization, meaning that the manifold is parametrized using coordinates in the computed eigenspace $V_e$. This excludes the possibility of folds on the manifold (however, if necessary we may circumvent this problem by increasing the dimensionality of the observed space). 

We seek the $2n\times m_M$ polynomial coefficient matrix $H$ of the manifold on the form

$y=V_e V_e^{\top } y+H\phi_{m,2:M} (V_e^{\top } y)$,

where the function $\phi_{m,2:M} (q)$ computes a vector of all $m_M$ monomials from orders 2 up to $M$ of an $m$-vector $q$. This optimization problem amounts to minimizing a cost function computed from the measured data $y_k$,

$$
C_r (V_e ,H)=\sum_{k=1}^N ||y_k -V_e V_e^{\top } y_k -H\phi_{m,2:M} (V_e^{\top } y_k )||
$$

to find $H$ and the $2n\times 2$ eigenspace matrix $V_e$, under the constraints that $V_e^{\top } H=0$ and $V_e^{\top } V_e =I$. The minimization is performed with the Matlab routine `fmincon` in `IMGeometry`. Since we know the exact eigenspace, we set this as input to the code.

```matlab:Code
SSMOrder = 3; SSMDim = 4;
[IMInfo, SSMChart, SSMFunction] = IMGeometry(yData(indTrain,:), SSMDim, SSMOrder,'V_e',[real(E(:,1)) imag(E(:,1)) real(E(:,2)) imag(E(:,2))]);

% Get reduced-order coordinates
etaData = projectTrajectories(IMInfo, yData);

% Plot static reconstruction
yStaticRec = yData(indTest,:);
yStaticRec {1,2} = IMInfo.parametrization.map(etaData{indTest,2});
plotTrajectories(yData(indTest,:), yStaticRec, 'm', fix(n/2), {'Test set', 'Static prediction'})
```

![figure_2.png
](README_images/figure_2.png
)

# Reduced order model

We compute a model for the reduced dynamics with the truncated training data projected onto the manifold. The function `IMDynamicsFlow` fits a polynomial map

$$
\dot{\eta} =R(\eta )=W_r \phi (\eta )
$$

where $\phi (\eta )$ again computes a vector of all monomials of $\eta$, and  $W_r$ is a matrix of polynomial coefficients. 

We are also specifying that we want the reduced dynamics in normal form, so that the map $N$ will fulfill

$$
\dot{z} =N(z)=Dz+W_n \phi (z)
$$

with $D$ a diagonal matrix and $W_n$ coefficients for the near-resonant nonlinear terms, after a near-identity change of coordinates

$$
z=T^{-1} (\eta )=\eta +W_t \phi (\eta )
$$

with $W_t$ containing the coefficients for the nonlinear terms of the coordinate change.

If we know $D$ from the model, we can enforce it in the data-driven dynamics identification. To do so, we fix the linear part of $R$. These dynamics act on coordinates being orthogonal projections on the eigenspace, which are not modal coordinates where the linear part of the dynamics is indeed the diagonal matrix $D$. However, these two coordinate systems are linked by a linear transformation. More precisely, we define

$P=V_e^{\top } E$ ,

where the columns of $E$ are the eigenvectors related to the SSM. So, we then have 

$DR(0)=PDP^{-1}$ .

We expect the polar normal form of the dynamics for this 4D nonresonant SSM case to take the form

${\dot{\rho} }_1 =c_1 (\rho_1 ,\rho_2 )\rho_1$ ,

${\dot{\rho} }_2 =c_2 (\rho_1 ,\rho_2 )\rho_2$ ,

${\dot{\theta} }_1 =\omega_1 (\rho_1 ,\rho_2 )$ ,

${\dot{\theta} }_2 =\omega_2 (\rho_1 ,\rho_2 )$ ,

which is identified from data using the function `IMDynamicsFlow`. In case of internal resonance, the polar vector field above also shows dependence on phases.

```matlab:Code
P = transpose(IMInfo.parametrization.tangentSpaceAtOrigin)*[E(:,1) E(:,2) E(:,1+n) E(:,2+n)];
DR0 = (P * diag(lambda([1 2 1+n 2+n])) ) / P ;
ROMOrder = 3; % Order for the normal form dynamics
RDInfo = IMDynamicsFlow(etaData(indTrain,:), 'R_PolyOrd', 1,'R_coeff', DR0, 'N_PolyOrd', ROMOrder, 'style', 'normalform');
```

```text:Output
 Done. 
Estimation of the reduced dynamics in normal form...
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
     0           1      1.16126e-05                      4.74e-05
     1           5      8.37699e-06            820        2.2e-05  
     2           6      6.41647e-06              1       1.11e-05  
     3           7      5.79993e-06              1       9.31e-06  
     4           8      5.13654e-06              1       8.71e-06  
     5           9      4.78871e-06              1       9.75e-06  
     6          10        4.424e-06              1        6.3e-06  
     7          11      4.17623e-06              1       5.81e-06  
     8          12      3.98872e-06              1       5.58e-06  
     9          13      3.81404e-06              1       6.16e-06  
    10          14      3.61171e-06              1       6.02e-06  
    11          15      3.44841e-06              1       4.28e-06  
    12          16      3.35985e-06              1       2.22e-06  
    13          17      3.30435e-06              1       2.32e-06  
    14          18      3.23291e-06              1        3.6e-06  
    15          19      3.14108e-06              1       4.39e-06  
    16          20      3.06459e-06              1       3.13e-06  
    17          21      3.02467e-06              1       1.72e-06  
    18          22      3.00141e-06              1       1.24e-06  
    19          23      2.97513e-06              1       2.23e-06  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    20          24      2.94407e-06              1       2.39e-06  
    21          25      2.91775e-06              1       2.02e-06  
    22          26      2.89825e-06              1       1.92e-06  
    23          27      2.87915e-06              1       1.49e-06  
    24          28      2.85671e-06              1       2.12e-06  
    25          29      2.83704e-06              1       1.57e-06  
    26          30       2.8252e-06              1       1.18e-06  
    27          31      2.81746e-06              1       1.02e-06  
    28          32      2.80865e-06              1       1.27e-06  
    29          33      2.79854e-06              1       1.38e-06  
    30          34      2.79124e-06              1       7.44e-07  
    31          35      2.78783e-06              1       4.09e-07  
    32          36      2.78564e-06              1       4.23e-07  
    33          37      2.78238e-06              1       7.68e-07  
    34          38      2.77718e-06              1       8.96e-07  
    35          39      2.77053e-06              1       1.15e-06  
    36          40      2.76451e-06              1       9.77e-07  
    37          41      2.76029e-06              1       4.33e-07  
    38          42       2.7571e-06              1       4.89e-07  
    39          43      2.75378e-06              1          7e-07  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    40          44       2.7498e-06              1       1.02e-06  
    41          45      2.74546e-06              1       9.57e-07  
    42          46      2.74172e-06              1       5.46e-07  
    43          47       2.7392e-06              1       3.68e-07  
    44          48      2.73742e-06              1       4.06e-07  
    45          49      2.73548e-06              1       5.85e-07  
    46          50      2.73279e-06              1       5.97e-07  
    47          51      2.72962e-06              1       6.47e-07  
    48          52      2.72697e-06              1       4.83e-07  
    49          53      2.72522e-06              1       4.08e-07  
    50          54      2.72384e-06              1       3.04e-07  
    51          55      2.72239e-06              1       4.23e-07  
    52          56      2.72112e-06              1       3.02e-07  
    53          57      2.72024e-06              1       2.68e-07  
    54          58      2.71953e-06              1       2.55e-07  
    55          59      2.71861e-06              1       3.59e-07  
    56          60      2.71741e-06              1       4.09e-07  
    57          61      2.71629e-06              1       2.54e-07  
    58          62      2.71557e-06              1       1.54e-07  
    59          63      2.71508e-06              1       1.48e-07  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    60          64      2.71448e-06              1       2.53e-07  
    61          65       2.7136e-06              1       2.92e-07  
    62          66      2.71259e-06              1       3.58e-07  
    63          67      2.71176e-06              1       2.56e-07  
    64          68      2.71119e-06              1       1.89e-07  
    65          69      2.71069e-06              1       1.68e-07  
    66          70      2.71013e-06              1       2.85e-07  
    67          71      2.70952e-06              1       3.25e-07  
    68          72      2.70902e-06              1       2.14e-07  
    69          73      2.70868e-06              1        1.3e-07  
    70          74      2.70844e-06              1       1.07e-07  
    71          75      2.70815e-06              1       2.14e-07  
    72          76      2.70771e-06              1       2.75e-07  
    73          77       2.7071e-06              1        2.7e-07  
    74          78      2.70654e-06              1        2.3e-07  
    75          79      2.70622e-06              1        9.2e-08  
    76          80      2.70606e-06              1       8.95e-08  
    77          81       2.7059e-06              1       1.29e-07  
    78          82      2.70564e-06              1       2.06e-07  
    79          83      2.70533e-06              1       1.99e-07  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    80          84      2.70512e-06              1       9.51e-08  
    81          85      2.70502e-06              1       6.85e-08  
    82          86      2.70493e-06              1       7.75e-08  
    83          87       2.7048e-06              1       1.41e-07  
    84          88       2.7046e-06              1       1.64e-07  
    85          89       2.7044e-06              1        1.1e-07  
    86          90      2.70425e-06              1       1.16e-07  
    87          91      2.70412e-06              1       1.09e-07  
    88          92      2.70392e-06              1       1.63e-07  
    89          93      2.70362e-06              1       2.15e-07  
    90          94      2.70332e-06              1       1.68e-07  
    91          95      2.70315e-06              1        6.8e-08  
    92          96      2.70306e-06              1       6.31e-08  
    93          97      2.70297e-06              1          1e-07  
    94          98      2.70283e-06              1        1.5e-07  
    95          99      2.70267e-06              1       1.34e-07  
    96         100      2.70258e-06              1        5.6e-08  
    97         101      2.70254e-06              1       3.84e-08  
    98         102      2.70252e-06              1       4.47e-08  
    99         103      2.70247e-06              1       7.91e-08  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   100         104      2.70241e-06              1       8.57e-08  
   101         105      2.70237e-06              1       4.99e-08  
   102         106      2.70234e-06              1       4.56e-08  
   103         107      2.70232e-06              1        4.2e-08  
   104         108      2.70229e-06              1       5.49e-08  
   105         109      2.70227e-06              1       5.66e-08  
   106         110      2.70224e-06              1       2.98e-08  
   107         111      2.70223e-06              1       3.13e-08  
   108         112      2.70222e-06              1       3.14e-08  
   109         113      2.70221e-06              1        4.5e-08  
   110         114      2.70218e-06              1       4.98e-08  
   111         115      2.70217e-06              1       2.98e-08  
   112         116      2.70216e-06              1       1.91e-08  
   113         117      2.70215e-06              1       1.86e-08  
   114         118      2.70214e-06              1       3.07e-08  
   115         119      2.70213e-06              1       3.53e-08  
   116         120      2.70212e-06              1       3.44e-08  
   117         121       2.7021e-06              1       3.42e-08  
   118         122       2.7021e-06              1       2.37e-08  
   119         123      2.70209e-06              1       2.61e-08  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   120         124      2.70208e-06              1       2.54e-08  
   121         125      2.70207e-06              1       2.78e-08  
   122         126      2.70206e-06              1       3.16e-08  
   123         127      2.70205e-06              1       2.69e-08  
   124         128      2.70204e-06              1       2.34e-08  
   125         129      2.70204e-06              1       2.77e-08  
   126         130      2.70203e-06              1        1.8e-08  
   127         131      2.70203e-06              1        8.3e-09  
   128         132      2.70203e-06              1       8.76e-09  
   129         133      2.70203e-06              1       1.91e-08  
   130         134      2.70202e-06              1       2.41e-08  
   131         135      2.70202e-06              1       2.03e-08  
   132         136      2.70201e-06              1       1.49e-08  
   133         137      2.70201e-06              1       1.57e-08  
   134         138        2.702e-06              1       2.33e-08  
   135         139        2.702e-06              1       2.82e-08  
   136         140      2.70199e-06              1       2.12e-08  
   137         141      2.70198e-06              1       1.66e-08  
   138         142      2.70198e-06              1       7.37e-09  
   139         143      2.70198e-06              1       9.38e-09  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   140         144      2.70198e-06              1       1.12e-08  
   141         145      2.70197e-06              1       1.74e-08  
   142         146      2.70197e-06              1       1.63e-08  
   143         147      2.70197e-06              1       1.15e-08  
   144         148      2.70197e-06              1       1.42e-08  
   145         149      2.70196e-06              1       1.58e-08  
   146         150      2.70196e-06              1       2.44e-08  
   147         151      2.70195e-06              1       2.93e-08  
   148         152      2.70194e-06              1       2.85e-08  
   149         153      2.70194e-06              1       1.62e-08  
   150         154      2.70193e-06              1       1.54e-08  
   151         155      2.70193e-06              1       1.54e-08  
   152         156      2.70192e-06              1       3.14e-08  
   153         157      2.70191e-06              1       4.31e-08  
   154         158       2.7019e-06              1       4.01e-08  
   155         159      2.70189e-06              1       2.06e-08  
   156         160      2.70189e-06              1       1.92e-08  
   157         161      2.70188e-06              1       2.74e-08  
   158         162      2.70187e-06              1       4.49e-08  
   159         163      2.70185e-06              1       5.07e-08  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   160         164      2.70184e-06              1       3.68e-08  
   161         165      2.70182e-06              1        2.2e-08  
   162         166      2.70182e-06              1       1.62e-08  
   163         167      2.70182e-06              1       1.21e-08  
   164         168      2.70181e-06              1       1.55e-08  
   165         169      2.70181e-06              1        1.9e-08  
   166         170       2.7018e-06              1       2.16e-08  
   167         171       2.7018e-06              1       2.32e-08  
   168         172       2.7018e-06              1       1.87e-08  
   169         173      2.70179e-06              1       1.82e-08  
   170         174      2.70179e-06              1       1.76e-08  
   171         175      2.70179e-06              1       1.09e-08  
   172         176      2.70179e-06              1       1.27e-08  
   173         177      2.70178e-06              1       1.32e-08  
   174         178      2.70178e-06              1       2.37e-08  
   175         179      2.70177e-06              1       2.73e-08  
   176         180      2.70177e-06              1       1.79e-08  
   177         181      2.70176e-06              1       1.32e-08  
   178         182      2.70176e-06              1       1.29e-08  
   179         183      2.70176e-06              1       2.32e-08  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   180         184      2.70175e-06              1       3.46e-08  
   181         185      2.70174e-06              1       3.41e-08  
   182         186      2.70173e-06              1       1.76e-08  
   183         187      2.70173e-06              1       1.44e-08  
   184         188      2.70173e-06              1       1.34e-08  
   185         189      2.70172e-06              1       2.18e-08  
   186         190      2.70172e-06              1       3.18e-08  
   187         191      2.70171e-06              1       3.09e-08  
   188         192       2.7017e-06              1       1.58e-08  
   189         193       2.7017e-06              1       8.27e-09  
   190         194       2.7017e-06              1       6.22e-09  
   191         195       2.7017e-06              1       9.53e-09  
   192         196      2.70169e-06              1       1.15e-08  
   193         197      2.70169e-06              1       1.12e-08  
   194         198      2.70169e-06              1       1.36e-08  
   195         199      2.70169e-06              1       1.26e-08  
   196         200      2.70169e-06              1       9.17e-09  
   197         201      2.70169e-06              1       1.33e-08  
   198         202      2.70168e-06              1       1.58e-08  
   199         203      2.70168e-06              1       1.08e-08  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   200         204      2.70168e-06              1       6.54e-09  
   201         205      2.70168e-06              1       6.41e-09  
   202         206      2.70168e-06              1       1.48e-08  
   203         207      2.70168e-06              1       1.89e-08  
   204         208      2.70167e-06              1       1.47e-08  
   205         209      2.70167e-06              1       6.01e-09  
   206         210      2.70167e-06              1       5.83e-09  
   207         211      2.70167e-06              1       1.08e-08  
   208         212      2.70167e-06              1       1.44e-08  
   209         213      2.70167e-06              1       1.26e-08  
   210         214      2.70167e-06              1       6.28e-09  
   211         215      2.70167e-06              1       5.69e-09  
   212         216      2.70167e-06              1       5.61e-09  
   213         217      2.70167e-06              1       6.48e-09  
   214         218      2.70166e-06              1       6.51e-09  
   215         219      2.70166e-06              1       5.54e-09  
   216         220      2.70166e-06              1       5.79e-09  
   217         221      2.70166e-06              1       6.08e-09  
   218         222      2.70166e-06              1       6.66e-09  
   219         223      2.70166e-06              1       1.01e-08  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   220         224      2.70166e-06              1       1.11e-08  
   221         225      2.70166e-06              1       6.85e-09  
   222         226      2.70166e-06              1       6.69e-09  
   223         227      2.70166e-06              1       5.29e-09  
   224         228      2.70166e-06              1       8.87e-09  
   225         229      2.70165e-06              1       1.21e-08  
   226         230      2.70165e-06              1       1.04e-08  
   227         231      2.70165e-06              1        7.4e-09  
   228         232      2.70165e-06              1          6e-09  
   229         233      2.70165e-06              1       6.32e-09  
   230         234      2.70165e-06              1       1.11e-08  
   231         235      2.70165e-06              1       1.31e-08  
   232         236      2.70164e-06              1       9.53e-09  
   233         237      2.70164e-06              1       8.39e-09  
   234         238      2.70164e-06              1       6.06e-09  
   235         239      2.70164e-06              1       5.91e-09  
   236         240      2.70164e-06              1       7.88e-09  
   237         241      2.70164e-06              1       7.14e-09  
   238         242      2.70164e-06              1       7.73e-09  
   239         243      2.70164e-06              1       7.46e-09  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   240         244      2.70164e-06              1       7.71e-09  
   241         245      2.70164e-06              1       8.44e-09  
   242         246      2.70164e-06              1       7.86e-09  
   243         247      2.70163e-06              1       5.06e-09  
   244         248      2.70163e-06              1       4.86e-09  
   245         249      2.70163e-06              1       5.31e-09  
   246         250      2.70163e-06              1       8.65e-09  
   247         251      2.70163e-06              1       1.22e-08  
   248         252      2.70163e-06              1       1.13e-08  
   249         253      2.70163e-06              1       5.38e-09  
   250         254      2.70163e-06              1       3.12e-09  
   251         255      2.70163e-06              1       2.98e-09  
   252         256      2.70163e-06              1       3.79e-09  
   253         257      2.70163e-06              1       6.42e-09  
   254         258      2.70162e-06              1       8.14e-09  
   255         259      2.70162e-06              1       6.61e-09  
   256         260      2.70162e-06              1       4.69e-09  
   257         261      2.70162e-06              1       4.47e-09  
   258         262      2.70162e-06              1       3.85e-09  
   259         263      2.70162e-06              1       4.93e-09  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   260         264      2.70162e-06              1       6.55e-09  
   261         265      2.70162e-06              1       5.66e-09  
   262         266      2.70162e-06              1       5.46e-09  
   263         267      2.70162e-06              1       5.38e-09  
   264         268      2.70162e-06              1       4.74e-09  
   265         269      2.70162e-06              1       5.61e-09  
   266         270      2.70162e-06              1       8.43e-09  
   267         271      2.70161e-06              1       9.27e-09  
   268         272      2.70161e-06              1       5.98e-09  
   269         273      2.70161e-06              1       3.03e-09  
   270         274      2.70161e-06              1       2.83e-09  
   271         275      2.70161e-06              1       2.67e-09  
   272         276      2.70161e-06              1       3.18e-09  
   273         277      2.70161e-06              1       4.34e-09  
   274         278      2.70161e-06              1       4.11e-09  
   275         279      2.70161e-06              1       5.86e-09  
   276         280      2.70161e-06              1       7.48e-09  
   277         281      2.70161e-06              1       7.93e-09  
   278         282      2.70161e-06              1       7.67e-09  
   279         283       2.7016e-06              1       6.17e-09  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   280         284       2.7016e-06              1       5.44e-09  
   281         285       2.7016e-06              1       3.99e-09  
   282         286       2.7016e-06              1       3.11e-09  
   283         287       2.7016e-06              1       2.54e-09  
   284         288       2.7016e-06              1       2.41e-09  
   285         289       2.7016e-06              1       2.63e-09  
   286         290       2.7016e-06              1       5.02e-09  
   287         291       2.7016e-06              1        7.7e-09  
   288         292       2.7016e-06              1       8.09e-09  
   289         293       2.7016e-06              1       7.02e-09  
   290         294       2.7016e-06              1       4.84e-09  
   291         295       2.7016e-06              1       3.93e-09  
   292         296       2.7016e-06              1       2.98e-09  
   293         297       2.7016e-06              1        2.9e-09  
   294         298       2.7016e-06              1       2.29e-09  
   295         299       2.7016e-06              1       2.15e-09  
   296         300       2.7016e-06              1       1.86e-09  
   297         301       2.7016e-06              1       1.71e-09  
   298         302       2.7016e-06              1       2.78e-09  
   299         303       2.7016e-06              1       2.32e-09  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   300         304       2.7016e-06              1       2.22e-09  
   301         305       2.7016e-06              1       2.13e-09  
   302         306       2.7016e-06              1       1.72e-09  
   303         307       2.7016e-06              1       1.24e-09  
   304         308       2.7016e-06              1       1.48e-09  
   305         309       2.7016e-06              1       1.16e-09  
   306         310       2.7016e-06              1        1.2e-09  
   307         311       2.7016e-06              1       1.05e-09  
   308         312       2.7016e-06              1       7.02e-10  
   309         313       2.7016e-06              1       6.81e-10  
   310         314       2.7016e-06              1       4.99e-10  
   311         315       2.7016e-06              1       4.44e-10  
   312         316       2.7016e-06              1       4.42e-10  
   313         317       2.7016e-06              1       4.37e-10  
   314         318       2.7016e-06              1       4.19e-10  
   315         319       2.7016e-06              1       5.26e-10  
   316         320       2.7016e-06              1       5.44e-10  
   317         321       2.7016e-06              1       5.99e-10  
   318         322       2.7016e-06              1       6.06e-10  
   319         323       2.7016e-06              1       5.24e-10  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   320         324       2.7016e-06              1       3.51e-10  
   321         325       2.7016e-06              1       3.33e-10  
   322         326       2.7016e-06              1       2.03e-10  
   323         327       2.7016e-06              1       1.56e-10  
   324         328       2.7016e-06              1       1.48e-10  
   325         329       2.7016e-06              1       1.36e-10  
   326         330       2.7016e-06              1       1.97e-10  
   327         331       2.7016e-06              1       2.02e-10  
   328         332       2.7016e-06              1       1.34e-10  
   329         333       2.7016e-06              1       7.98e-11  
   330         334       2.7016e-06              1       5.85e-11  
   331         335       2.7016e-06              1       5.51e-11  

Local minimum possible.

fminunc stopped because the size of the current step is less than
the value of the step size tolerance.

<stopping criteria details>
Plotting figure with the polar normal form equations ... Done. 
```

![figure_3.png
](README_images/figure_3.png
)

We transform the truncated initial condition of our test trajectory according to the obtained change of coordinates, and integrate our reduced order evolution rule to predict the development of the trajectory.

```matlab:Code
[yRec, etaRec, zRec] = advect(IMInfo, RDInfo, yData);
```

# Evaluation of reduced dynamics

The error RRMSE is computed as the average distance of the predicted trajectory to the measured one in the full state space. We also plot the testing trajectory

```matlab:Code
normedTrajDist = computeTrajectoryErrors(yRec, yData);
RRMSE = mean(normedTrajDist(indTest))
```

```text:Output
RRMSE = 
     0.026549

```

```matlab:Code

% Plot reduced coordinates
plotReducedCoordinates(etaData(indTest,:), etaRec(indTest,:))
view(3)
customFigure();

% Plot normal form coordinates
T_1 = RDInfo.inverseTransformation.map;
zDataTest = T_1(etaData{indTest,2});
zRecTest = zRec{indTest,2};
plot3(real(zDataTest(1,:)),imag(zDataTest(1,:)),real(zDataTest(2,:)),'Linewidth',2)
plot3(real(zRecTest(1,:)),imag(zRecTest(1,:)),real(zRecTest(2,:)),'--','Linewidth',2)
xlabel('$\rho_1\cos\theta_1$','Interpreter','latex')
ylabel('$\rho_1\sin\theta_1$','Interpreter','latex')
zlabel('$\rho_2\cos\theta_2$','Interpreter','latex')
legend('Test data','Prediction')
view(3)
```

![figure_4.png
](README_images/figure_4.png
)

```matlab:Code

% Plot physical coordinates
plotTrajectories(yData(indTest,:), yRec(indTest,:), 'm', n, {'Test set', 'Prediction'})
```

![figure_5.png
](README_images/figure_5.png
)

# Backbone curves and frequency responses

We can now study the reduced order model. We extract backbones curves and forced response curves (FRCs). For the latter, our 

$\dot{z} =N(z)+f_{red} e^{i\Omega t}$,

where $\Omega$ is the forcing frequency and $f_{red}$ is the forcing amplitude vector. For FRC extraction, we use the numerical continuation core coco. In this multimodal case, the backbone curves are those of the uncoupled oscillator limit, i.e., the backbones curves of a certain mode are computed setting to zero the amplitudes of the other modes.

```matlab:Code
% Define the physical amplitude metric
reconstructedEigenvalues = RDInfo.eigenvaluesLinPartFlow; ndofSSM = SSMDim/2; 
amplitudeFunction = @(x) x(n,:); maxRho = max(abs(zRec{1,2}(1:2,:)),[],2);

% Compute backbone curves in the training range of data
BBCInfo = backboneCurves(IMInfo, RDInfo, amplitudeFunction,maxRho, 'norm');
```

![figure_6.png
](README_images/figure_6.png
)

![figure_7.png
](README_images/figure_7.png
)

```matlab:Code
BBSInfo = backboneSurfaces(RDInfo, maxRho, 'norm');
```

![figure_8.png
](README_images/figure_8.png
)

![figure_9.png
](README_images/figure_9.png
)

```matlab:Code
% Compute frequency response curves
omegaSpan = [0.95 1.05].*transpose(abs(reconstructedEigenvalues(1:ndofSSM))); 
forcingSpan = [1 2.5 5; [.5 1.25 2.5]*3]*1e-4; % Assumes forcings on different dofs to have the same phase
[IMInfoF,RDInfoF] = forcedSSMROM(IMInfo,RDInfo,'nForcingFrequencies',1);
```

```text:Output
Forced SSM reduced-order model assumes external forcing only along the tangent (modal) subspace at the origin. 
```

```matlab:Code
FRCDataCOCO = continuationFRCpo(IMInfoF,RDInfoF,forcingSpan,omegaSpan,amplitudeFunction);
```

```text:Output
Frequency sweep for the set of forcing values number 1 ...

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          1.48e-04  4.68e+01    0.0    0.0    0.0
   1   1  1.00e+00  1.06e-04  7.47e-13  4.68e+01    0.0    0.0    0.0
   2   1  1.00e+00  3.66e-11  1.89e-17  4.68e+01    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE     po.period            T
    0  00:00:00   4.6817e+01      1  EP      2.3409e+01   2.3409e+01
   10  00:00:00   4.4736e+01      2          2.2367e+01   2.2367e+01
   20  00:00:01   4.4568e+01      3          2.2283e+01   2.2283e+01
   30  00:00:02   4.4473e+01      4          2.2235e+01   2.2235e+01
   40  00:00:02   4.4388e+01      5          2.2193e+01   2.2193e+01
   50  00:00:03   4.4267e+01      6          2.2133e+01   2.2133e+01
   60  00:00:03   4.3971e+01      7          2.1986e+01   2.1986e+01
   70  00:00:04   4.0637e+01      8          2.0319e+01   2.0319e+01
   80  00:00:04   3.5637e+01      9          1.7819e+01   1.7819e+01
   90  00:00:05   3.0637e+01     10          1.5319e+01   1.5319e+01
  100  00:00:05   2.5637e+01     11          1.2819e+01   1.2819e+01
  110  00:00:06   2.0637e+01     12          1.0319e+01   1.0319e+01
  120  00:00:07   1.6233e+01     13          8.1164e+00   8.1164e+00
  130  00:00:07   1.6066e+01     14          8.0324e+00   8.0324e+00
  140  00:00:08   1.6024e+01     15          8.0105e+00   8.0105e+00
  150  00:00:08   1.5996e+01     16          7.9962e+00   7.9962e+00
  160  00:00:09   1.5972e+01     17          7.9848e+00   7.9848e+00
  170  00:00:09   1.5944e+01     18          7.9711e+00   7.9711e+00
  180  00:00:10   1.5881e+01     19          7.9405e+00   7.9405e+00
  190  00:00:10   1.5631e+01     20          7.8155e+00   7.8155e+00
  193  00:00:11   1.5246e+01     21  EP      7.6229e+00   7.6229e+00
Frequency sweep for the set of forcing values number 2 ...

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          6.19e-04  4.68e+01    0.0    0.0    0.0
   1   1  1.00e+00  1.43e-04  4.39e-12  4.68e+01    0.0    0.0    0.0
   2   1  1.00e+00  2.22e-10  4.67e-17  4.68e+01    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE     po.period            T
    0  00:00:00   4.6817e+01      1  EP      2.3409e+01   2.3409e+01
   10  00:00:00   4.4908e+01      2          2.2453e+01   2.2453e+01
   20  00:00:01   4.4632e+01      3          2.2311e+01   2.2311e+01
   30  00:00:01   4.4492e+01      4          2.2236e+01   2.2236e+01
   40  00:00:02   4.4400e+01      5          2.2190e+01   2.2190e+01
   50  00:00:03   4.4326e+01      6          2.2154e+01   2.2154e+01
   60  00:00:03   4.4218e+01      7          2.2105e+01   2.2105e+01
   70  00:00:04   4.3958e+01      8          2.1978e+01   2.1978e+01
   80  00:00:05   4.2754e+01      9          2.1377e+01   2.1377e+01
   90  00:00:05   3.7755e+01     10          1.8877e+01   1.8877e+01
  100  00:00:06   3.2755e+01     11          1.6377e+01   1.6377e+01
  110  00:00:06   2.7755e+01     12          1.3877e+01   1.3877e+01
  120  00:00:07   2.2755e+01     13          1.1377e+01   1.1377e+01
  130  00:00:07   1.7755e+01     14          8.8773e+00   8.8773e+00
  140  00:00:08   1.6166e+01     15          8.0821e+00   8.0821e+00
  150  00:00:08   1.6054e+01     16          8.0235e+00   8.0235e+00
  160  00:00:09   1.5993e+01     17          7.9907e+00   7.9907e+00
  170  00:00:10   1.5956e+01     18          7.9686e+00   7.9686e+00
  180  00:00:10   1.5930e+01     19          7.9570e+00   7.9570e+00
  190  00:00:11   1.5910e+01     20          7.9501e+00   7.9501e+00
  200  00:00:12   1.5855e+01     21          7.9262e+00   7.9262e+00
  210  00:00:13   1.5696e+01     22          7.8479e+00   7.8479e+00
  216  00:00:13   1.5246e+01     23  EP      7.6229e+00   7.6229e+00
Frequency sweep for the set of forcing values number 3 ...

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          4.25e-04  4.68e+01    0.0    0.0    0.0
   1   1  1.00e+00  2.64e-04  2.45e-11  4.68e+01    0.0    0.0    0.0
   2   1  1.00e+00  1.08e-09  8.14e-17  4.68e+01    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE     po.period            T
    0  00:00:00   4.6818e+01      1  EP      2.3409e+01   2.3409e+01
   10  00:00:00   4.5019e+01      2          2.2505e+01   2.2505e+01
   20  00:00:01   4.4613e+01      3          2.2288e+01   2.2288e+01
   30  00:00:02   4.4392e+01      4          2.2161e+01   2.2161e+01
   40  00:00:03   4.4259e+01      5          2.2087e+01   2.2087e+01
   50  00:00:04   4.4199e+01      6          2.2062e+01   2.2062e+01
   60  00:00:05   4.4160e+01      7          2.2058e+01   2.2058e+01
   70  00:00:06   4.3990e+01      8          2.1988e+01   2.1988e+01
   80  00:00:07   4.3472e+01      9          2.1734e+01   2.1734e+01
   90  00:00:08   4.0301e+01     10          2.0150e+01   2.0150e+01
  100  00:00:09   3.5301e+01     11          1.7651e+01   1.7651e+01
  110  00:00:09   3.0301e+01     12          1.5151e+01   1.5151e+01
  120  00:00:10   2.5301e+01     13          1.2651e+01   1.2651e+01
  130  00:00:10   2.0301e+01     14          1.0151e+01   1.0151e+01
  140  00:00:11   1.6332e+01     15          8.1647e+00   8.1647e+00
  150  00:00:12   1.6097e+01     16          8.0430e+00   8.0430e+00
  160  00:00:12   1.5992e+01     17          7.9813e+00   7.9813e+00
  170  00:00:13   1.5914e+01     18          7.9362e+00   7.9362e+00
  180  00:00:13   1.5865e+01     19          7.9114e+00   7.9114e+00
  187  00:00:14   1.5853e+01     20  FP      7.9079e+00   7.9079e+00
  187  00:00:14   1.5853e+01     21  SN      7.9079e+00   7.9079e+00
  190  00:00:14   1.5850e+01     22          7.9088e+00   7.9088e+00
  197  00:00:15   1.5844e+01     23  FP      7.9109e+00   7.9109e+00
  197  00:00:15   1.5844e+01     24  SN      7.9109e+00   7.9109e+00
  200  00:00:16   1.5837e+01     25          7.9098e+00   7.9098e+00
  210  00:00:16   1.5743e+01     26          7.8695e+00   7.8695e+00
  220  00:00:17   1.5445e+01     27          7.7220e+00   7.7220e+00
  223  00:00:17   1.5246e+01     28  EP      7.6229e+00   7.6229e+00
```

```matlab:Code

freqNorm = BBCInfo.frequency(1,1);
% Plot backbone curves and FRCs
customFigure('subPlot',[1 ndofSSM]); colors = colororder;
for idof = 1:ndofSSM
    subplot(1,ndofSSM,idof)
    plot(BBCInfo.frequency(idof,:)/freqNorm, BBCInfo.amplitude(idof,:),'k',...
        'Linewidth',1,'DisplayName', 'Backbones')
    plotFRC(FRCDataCOCO, colors(1,:),', 4DSSM','freqscale',freqNorm)
    ylim([0 max(max(BBCInfo.amplitude))])
    xlim([0.96 1.04].*transpose(abs(reconstructedEigenvalues(idof)))/freqNorm); 
    legend('off')
    xlabel('$\Omega/\omega(0)$','Interpreter','latex'); 
    ylabel('$q_5 \, [$m$]$','Interpreter','latex');
end
legend('location','north')
```

![figure_10.png
](README_images/figure_10.png
)

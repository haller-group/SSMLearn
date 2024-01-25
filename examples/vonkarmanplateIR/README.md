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

```matlab:Code

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
Mode shape
  -0.00014705
  -0.00014181
   2.7274e-18
    -0.069108
  -0.00078624
   8.0752e-18
     -0.12208
  -0.00081778
   8.9702e-18
     -0.14761
   -0.0005186
  -2.5233e-17
     -0.14319
   0.00013401
  -6.0444e-19
     -0.11506
   0.00088987
  -3.3621e-18
    -0.075457
    0.0015077
   -5.682e-18
    -0.037929
    0.0016542
  -1.9181e-18
    -0.012192
     0.001569
  -2.0465e-17
   -0.0011337
   0.00037124
   1.3782e-17
    0.0002299
    0.0002299
  -3.0398e-19
  -0.00078592
    -0.069069
   5.1454e-18
   2.2008e-19
  -1.3482e-18
    0.0065286
    -0.061123
    -0.061154
  -1.8717e-17
  -1.2506e-18
  -9.1481e-19
     0.011494
     -0.10659
    -0.038893
   6.3512e-18
  -1.8541e-18
  -1.2021e-18
     0.013834
     -0.12686
   -0.0095195
    8.948e-18
  -3.4444e-18
   3.6607e-18
     0.013321
     -0.12005
     0.018009
   4.1543e-19
  -2.4455e-18
   1.8959e-18
     0.010575
    -0.092426
     0.035879
   1.0221e-18
  -8.9804e-19
  -3.9333e-20
     0.006786
    -0.055759
     0.040218
  -7.2178e-18
  -1.1865e-18
   7.1516e-19
     0.003255
    -0.022719
     0.032259
   1.9007e-19
  -1.5622e-18
   1.1834e-18
   0.00091105
   -0.0021867
     0.017601
    5.628e-18
   9.8776e-19
   1.3308e-18
   2.8511e-15
    0.0038327
    0.0038327
   3.0229e-18
   0.00037124
   -0.0011337
  -1.2367e-17
   -0.0008192
     -0.12197
  -8.7252e-18
   2.8361e-19
   -3.203e-18
      0.01149
    -0.038839
      -0.1066
  -6.2259e-18
   1.3117e-18
  -2.1906e-19
     0.020075
    -0.065316
    -0.065381
   2.0597e-17
  -3.3194e-19
  -1.6257e-18
      0.02384
    -0.072937
     -0.01181
  -1.1004e-17
  -2.4291e-18
    1.177e-18
     0.022427
    -0.060971
     0.037223
   8.1485e-21
  -9.0425e-19
   1.7833e-18
     0.017062
    -0.035336
     0.067368
   1.5289e-17
   -1.933e-18
   6.6837e-19
     0.010027
     -0.00621
     0.071873
  -4.3927e-18
  -1.7847e-18
   1.3739e-18
    0.0037698
     0.015883
     0.053713
  -4.2469e-18
   8.2615e-19
  -7.2801e-19
   8.7955e-15
     0.024053
     0.024053
   9.5616e-18
   8.5574e-19
   5.5559e-19
  -0.00091105
     0.017601
   -0.0021867
  -1.5315e-18
     0.001569
    -0.012192
   6.6198e-18
  -0.00052062
     -0.14743
    9.681e-18
    9.342e-19
  -2.8897e-18
     0.013825
   -0.0094745
     -0.12684
   6.8843e-18
    1.286e-18
  -2.6989e-19
     0.023833
    -0.011745
    -0.072994
  -2.4486e-17
   1.5552e-18
    -4.53e-18
     0.027584
   -0.0039053
   -0.0039752
   1.6636e-18
   1.6831e-18
  -8.8143e-19
     0.024727
     0.013346
     0.057426
   1.0111e-17
   8.6074e-19
   9.4864e-19
     0.017039
     0.034938
     0.092135
  -8.5164e-18
  -1.9918e-19
  -1.0587e-19
    0.0076855
     0.053236
     0.091657
   1.3472e-17
  -6.3734e-19
   2.3578e-18
   1.2306e-14
     0.060895
     0.060895
   2.3369e-17
   2.6172e-18
   5.2015e-19
   -0.0037698
     0.053713
     0.015883
   -1.054e-17
   2.5468e-18
  -6.7818e-19
    -0.003255
     0.032259
    -0.022719
  -9.8757e-18
    0.0016542
    -0.037929
  -2.4131e-18
     0.000131
     -0.14299
  -1.2665e-17
  -1.8551e-18
  -2.2935e-18
      0.01331
     0.018026
     -0.12001
  -6.2941e-18
  -3.1406e-18
  -1.9626e-19
     0.022416
     0.037256
    -0.061013
  -2.8427e-18
  -1.6706e-18
  -3.8717e-19
     0.024722
     0.057468
      0.01328
   2.1139e-17
   4.1873e-19
   9.4545e-19
     0.020018
     0.076821
     0.076776
   9.9855e-18
   1.8723e-18
   1.4228e-18
     0.010523
      0.09172
      0.10795
  -2.5157e-18
   1.9436e-18
   1.9156e-18
   9.3695e-15
     0.097842
     0.097842
   5.7914e-18
   2.5297e-18
   1.7861e-18
   -0.0076855
     0.091657
     0.053236
  -2.9834e-18
  -4.2831e-20
   2.0266e-19
    -0.010027
     0.071873
     -0.00621
  -2.0798e-17
   6.1908e-19
   9.3858e-19
    -0.006786
     0.040218
    -0.055759
  -1.0684e-18
    0.0015077
    -0.075457
  -7.4458e-19
   0.00088636
     -0.11487
   9.5898e-19
  -3.1511e-19
  -2.7799e-18
     0.010563
     0.035866
    -0.092377
   -1.927e-18
  -3.4006e-19
  -2.2813e-18
      0.01705
     0.067359
    -0.035367
  -1.8384e-19
   -8.981e-19
  -9.8439e-19
     0.017031
     0.092132
     0.034879
   9.0097e-18
  -4.8489e-19
   3.9527e-19
     0.010521
      0.10795
     0.091678
  -1.7725e-18
   6.5075e-19
   -2.147e-18
   1.0311e-15
      0.11339
      0.11339
    -1.44e-17
  -1.0065e-19
    2.028e-18
    -0.010523
      0.10795
      0.09172
   1.6626e-17
   5.6793e-19
   3.2685e-18
    -0.017039
     0.092135
     0.034938
   3.1518e-18
  -1.8002e-18
    1.683e-19
    -0.017062
     0.067368
    -0.035336
  -2.0053e-17
  -1.8087e-18
   1.9438e-18
    -0.010575
     0.035879
    -0.092426
   -2.691e-18
   0.00088987
     -0.11506
   8.0074e-18
    0.0015041
    -0.075315
  -3.2347e-17
  -2.1995e-18
   -2.423e-18
    0.0067768
     0.040183
     -0.05572
   -9.316e-18
  -1.2968e-18
  -2.3287e-18
     0.010017
      0.07183
   -0.0062382
   7.7168e-18
   1.0549e-18
  -2.2094e-18
    0.0076805
     0.091615
     0.053181
   2.4024e-17
   3.4111e-19
   -2.285e-18
  -7.1417e-15
     0.097802
     0.097802
   -1.561e-17
   8.3253e-19
  -4.0109e-18
    -0.010521
     0.091678
      0.10795
   9.9637e-18
   8.7799e-19
   1.8019e-18
    -0.020018
     0.076776
     0.076821
    2.169e-17
  -6.0335e-19
   2.7046e-18
    -0.024727
     0.057426
     0.013346
   -3.497e-17
  -3.6295e-18
  -1.2511e-19
    -0.022427
     0.037223
    -0.060971
  -2.1999e-18
   -1.155e-18
   2.4959e-18
    -0.013321
     0.018009
     -0.12005
   1.9453e-17
   0.00013401
     -0.14319
  -1.0149e-17
    0.0016509
    -0.037848
   1.0733e-18
   1.7428e-18
  -1.5839e-18
    0.0032497
     0.032218
      -0.0227
    1.857e-17
   3.0022e-18
     -2.8e-19
    0.0037652
     0.053661
     0.015851
   7.6128e-18
    1.239e-18
   2.0473e-18
  -1.0232e-14
     0.060842
     0.060842
  -2.3691e-17
   1.4974e-19
  -3.7433e-18
   -0.0076805
     0.053181
     0.091615
  -2.1591e-17
   4.9906e-19
  -2.7227e-18
    -0.017031
     0.034879
     0.092132
    1.115e-17
  -5.1673e-19
   1.3432e-18
    -0.024722
      0.01328
     0.057468
  -9.3424e-18
  -3.7245e-18
   5.2483e-19
    -0.027584
   -0.0039752
   -0.0039053
  -4.9576e-18
  -2.1869e-18
    1.571e-18
     -0.02384
     -0.01181
    -0.072937
   3.1546e-17
   3.0567e-18
   1.3516e-18
    -0.013834
   -0.0095195
     -0.12686
  -5.1027e-18
   -0.0005186
     -0.14761
  -1.6096e-17
    0.0015659
    -0.012162
   1.3162e-18
   -6.208e-19
   9.9114e-19
   0.00090928
     0.017572
   -0.0021876
   6.5966e-18
   5.2564e-19
   6.5516e-19
  -7.4559e-15
     0.024018
     0.024018
  -5.2283e-18
  -1.2017e-18
   1.6459e-18
   -0.0037652
     0.015851
     0.053661
  -2.2408e-17
  -1.2016e-18
  -1.2038e-18
    -0.010017
   -0.0062382
      0.07183
   1.0184e-17
   -2.411e-19
  -1.3797e-18
     -0.01705
    -0.035367
     0.067359
  -3.0277e-18
  -3.0044e-19
   3.8971e-19
    -0.022416
    -0.061013
     0.037256
  -4.5282e-18
  -2.9115e-18
   5.0394e-19
    -0.023833
    -0.072994
    -0.011745
   1.2538e-17
  -2.2965e-18
   1.6467e-18
    -0.020075
    -0.065381
    -0.065316
   1.0008e-17
   7.5361e-19
   1.3853e-18
    -0.011494
    -0.038893
     -0.10659
  -1.7088e-17
  -0.00081778
     -0.12208
   1.1401e-17
   0.00037016
   -0.0011306
   5.4378e-18
  -8.5866e-19
   1.1883e-18
   -2.414e-15
    0.0038237
    0.0038237
   4.5486e-19
  -3.1661e-18
   2.5015e-18
  -0.00090928
   -0.0021876
     0.017572
  -7.2492e-18
  -2.6039e-18
   1.1324e-18
   -0.0032497
      -0.0227
     0.032218
   2.0196e-18
  -1.9324e-18
   2.0595e-18
   -0.0067768
     -0.05572
     0.040183
   1.4305e-17
  -1.3128e-18
   8.9638e-19
    -0.010563
    -0.092377
     0.035866
  -2.6544e-18
  -1.3822e-19
   6.3449e-19
     -0.01331
     -0.12001
     0.018026
  -9.0213e-18
  -5.2404e-19
  -3.5435e-19
    -0.013825
     -0.12684
   -0.0094745
   5.5163e-18
  -1.8879e-18
   7.4014e-19
     -0.01149
      -0.1066
    -0.038839
   -5.688e-18
  -2.2833e-18
   2.2794e-19
   -0.0065286
    -0.061154
    -0.061123
  -8.7474e-18
  -0.00078624
    -0.069108
   2.2523e-17
   0.00022941
   0.00022941
    2.367e-18
   -0.0011306
   0.00037016
   1.9951e-18
    -0.012162
    0.0015659
   1.0638e-17
    -0.037848
    0.0016509
   6.7918e-18
    -0.075315
    0.0015041
  -6.1356e-18
     -0.11487
   0.00088636
   9.5252e-18
     -0.14299
     0.000131
  -9.1894e-18
     -0.14743
  -0.00052062
    7.412e-18
     -0.12197
   -0.0008192
   7.4862e-18
    -0.069069
  -0.00078592
   4.6513e-18
  -0.00014181
  -0.00014705
  -2.4456e-17

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
Mode shape
   0.00086981
  -0.00086881
   1.1961e-19
    0.0021813
  -0.00082061
   2.6942e-19
     0.013476
   -0.0015382
   2.5018e-19
     0.038394
   -0.0014675
  -7.6041e-19
     0.075069
   -0.0010419
  -4.3955e-20
      0.11411
  -0.00035574
  -1.6904e-19
      0.14218
   0.00027782
  -1.6181e-19
       0.1469
   0.00065313
   -7.516e-20
      0.12173
   0.00055752
  -5.3815e-19
     0.069069
   0.00029876
   4.0785e-19
   0.00021654
  -0.00021654
  -1.1683e-20
   0.00082087
   -0.0021986
   1.5146e-19
   5.1067e-21
  -4.2387e-20
   1.0145e-06
    -0.003085
    0.0030639
   -5.012e-19
  -3.9857e-20
  -2.7017e-20
  -0.00094763
    0.0030592
      0.01504
   1.5878e-19
  -6.3305e-20
  -3.8177e-20
   -0.0032802
     0.023484
     0.029219
   2.3562e-19
  -1.0228e-19
   1.1033e-19
   -0.0067896
     0.056362
     0.037822
   3.6029e-21
  -8.3362e-20
   6.1709e-20
    -0.010567
     0.092836
     0.034856
   1.3954e-20
  -3.4326e-20
   1.8148e-21
    -0.013318
      0.12023
     0.018405
  -2.1644e-19
  -3.2076e-20
   2.1023e-20
     -0.01384
      0.12671
    -0.008221
   3.8584e-21
  -3.8343e-20
   3.3967e-20
    -0.011503
      0.10601
    -0.037486
   1.6194e-19
   2.6953e-20
   3.9286e-20
   -0.0065354
     0.060191
    -0.060191
   8.5012e-20
  -0.00029876
    -0.069069
  -3.6224e-19
    0.0015381
    -0.013505
  -2.6911e-19
   7.4791e-21
  -9.2071e-20
   0.00095163
    -0.015056
   -0.0031069
  -1.8383e-19
   3.3832e-20
   -1.261e-20
   4.5059e-06
    -0.022508
      0.02246
   5.6679e-19
  -1.8578e-20
  -4.7603e-20
   -0.0037912
    -0.015045
     0.051929
  -3.0084e-19
  -8.4889e-20
   3.2907e-20
    -0.010072
    0.0065408
     0.070725
   2.9751e-20
  -2.5975e-20
   5.4386e-20
    -0.017136
     0.035336
     0.067284
   4.1837e-19
  -6.5176e-20
   2.2967e-20
    -0.022524
     0.060714
     0.038053
  -1.1322e-19
  -4.9712e-20
     4.19e-20
    -0.023934
     0.072386
    -0.010614
   -1.054e-19
   2.1777e-20
  -1.3799e-20
    -0.020133
     0.064383
    -0.064383
    2.757e-19
    2.552e-20
   1.8338e-20
    -0.011503
     0.037486
     -0.10601
  -5.9868e-20
  -0.00055752
     -0.12173
   1.5746e-19
    0.0014667
    -0.038414
    2.586e-19
    2.675e-20
  -8.4669e-20
    0.0032849
    -0.029216
    -0.023553
   1.6743e-19
   3.4049e-20
  -1.6104e-20
     0.003804
    -0.051952
     0.014959
  -6.4406e-19
   4.0445e-20
   -1.284e-19
   1.0208e-05
    -0.060419
     0.060357
   6.0691e-20
     4.29e-20
  -2.5672e-20
   -0.0077284
    -0.053714
     0.091914
   3.1818e-19
    1.837e-20
   3.4747e-20
     -0.01715
    -0.035939
      0.09319
  -1.9848e-19
  -1.1422e-20
   3.1449e-21
    -0.024889
    -0.014514
     0.058861
   3.9721e-19
  -1.9869e-20
   7.8913e-20
    -0.027746
    0.0027569
   -0.0027569
   6.4607e-19
    7.712e-20
   1.9968e-20
    -0.023934
     0.010614
    -0.072386
  -2.8015e-19
   7.3728e-20
  -1.3583e-20
     -0.01384
     0.008221
     -0.12671
  -3.1871e-19
  -0.00065313
      -0.1469
  -8.0225e-20
    0.0010403
    -0.075054
   -3.984e-19
  -5.1455e-20
  -7.6027e-20
    0.0067929
    -0.037794
    -0.056436
  -1.9361e-19
  -8.2421e-20
   -1.536e-20
     0.010085
    -0.070705
   -0.0066547
  -7.9555e-20
   -5.316e-20
  -1.7225e-20
    0.0077526
    -0.091923
     0.053613
   6.4481e-19
   9.8001e-21
   2.2353e-20
   1.5818e-05
    -0.099213
     0.099168
   3.5172e-19
    5.111e-20
   4.3047e-20
    -0.010596
    -0.093628
      0.10989
  -2.2548e-20
   5.2021e-20
   6.1285e-20
    -0.020168
     -0.07868
      0.07868
   1.7147e-19
   6.6961e-20
   6.0824e-20
    -0.024889
    -0.058861
     0.014514
  -8.4302e-20
  -1.8988e-21
   1.5839e-20
    -0.022524
    -0.038053
    -0.060714
  -6.4657e-19
   1.6842e-20
   2.9953e-20
    -0.013318
    -0.018405
     -0.12023
  -6.7887e-20
  -0.00027782
     -0.14218
  -5.8918e-20
   0.00035334
     -0.11403
   8.4996e-22
  -1.0212e-20
  -7.6366e-20
     0.010567
    -0.034807
    -0.092894
   -8.012e-20
  -1.0345e-20
  -6.7261e-20
     0.017145
    -0.067222
     -0.03546
   2.7089e-23
  -2.5529e-20
   -3.427e-20
     0.017172
    -0.093135
     0.035812
   2.9401e-19
  -1.3235e-20
   7.1646e-21
     0.010628
     -0.10985
     0.093549
  -2.0552e-20
   1.7359e-20
  -5.4654e-20
   1.8178e-05
     -0.11564
      0.11564
  -3.7415e-19
  -5.7109e-21
   6.6554e-20
    -0.010596
     -0.10989
     0.093628
   4.4534e-19
   1.2411e-20
   9.0014e-20
     -0.01715
     -0.09319
     0.035939
    7.481e-20
  -5.9098e-20
   9.7971e-21
    -0.017136
    -0.067284
    -0.035336
  -5.6706e-19
  -5.8674e-20
   5.1682e-20
    -0.010567
    -0.034856
    -0.092836
  -1.1654e-19
   0.00035574
     -0.11411
   2.0059e-19
  -0.00028095
     -0.14204
  -9.4785e-19
  -6.3508e-20
  -6.8971e-20
     0.013313
     -0.01835
     -0.12025
  -2.6563e-19
  -3.4844e-20
  -8.1332e-20
     0.022525
    -0.037966
    -0.060824
   2.0921e-19
   3.0686e-20
   -5.968e-20
     0.024903
    -0.058758
     0.014382
   7.1545e-19
    1.117e-20
  -6.8765e-20
     0.020194
     -0.07858
      0.07858
  -4.2842e-19
   2.2643e-20
  -1.0623e-19
     0.010628
    -0.093549
      0.10985
   3.1717e-19
   2.1801e-20
   5.9307e-20
   1.5818e-05
    -0.099168
     0.099213
   6.0009e-19
  -2.5836e-20
   7.4385e-20
   -0.0077284
    -0.091914
     0.053714
  -1.0004e-18
  -1.0129e-19
   2.6245e-21
    -0.010072
    -0.070725
   -0.0065408
   -6.175e-20
  -3.5229e-20
   7.4708e-20
   -0.0067896
    -0.037822
    -0.056362
   5.1383e-19
    0.0010419
    -0.075069
  -3.0881e-19
   -0.0006561
     -0.14671
   1.2736e-20
   4.8574e-20
  -4.0838e-20
      0.01383
    0.0082612
     -0.12669
   5.2382e-19
   8.3815e-20
  -1.3181e-20
     0.023927
     0.010694
    -0.072462
   2.0978e-19
   3.4265e-20
   5.1271e-20
     0.027749
    0.0028716
   -0.0028716
  -6.4716e-19
   5.0259e-21
  -1.1447e-19
     0.024903
    -0.014382
     0.058758
  -6.0879e-19
   1.3199e-20
  -7.8932e-20
     0.017172
    -0.035812
     0.093135
   3.2537e-19
  -1.9442e-20
   4.4374e-20
    0.0077526
    -0.053613
     0.091923
  -2.3994e-19
  -1.1074e-19
   3.0253e-20
   1.0208e-05
    -0.060357
     0.060419
  -1.2808e-19
  -7.0103e-20
   4.9271e-20
   -0.0037912
    -0.051929
     0.015045
   8.8263e-19
   8.0813e-20
   3.1961e-20
   -0.0032802
    -0.029219
    -0.023484
  -1.6379e-19
    0.0014675
    -0.038394
  -4.3582e-19
  -0.00056145
     -0.12153
    3.124e-20
  -1.8477e-20
   2.6758e-20
      0.01149
     0.037489
     -0.10596
    1.932e-19
   1.2891e-20
   1.6523e-20
      0.02012
     0.064418
    -0.064418
  -1.4978e-19
  -3.3822e-20
   4.0705e-20
     0.023927
     0.072462
    -0.010694
  -6.0518e-19
  -3.3314e-20
  -3.6737e-20
     0.022525
     0.060824
     0.037966
   2.7329e-19
  -1.0394e-20
  -3.8213e-20
     0.017145
      0.03546
     0.067222
  -7.6033e-20
  -9.8093e-21
   1.2664e-20
     0.010085
    0.0066547
     0.070705
  -1.1247e-19
  -8.7293e-20
   2.1613e-20
     0.003804
    -0.014959
     0.051952
   3.4479e-19
  -6.0857e-20
   5.4461e-20
   4.5059e-06
     -0.02246
     0.022508
   3.0222e-19
    1.745e-20
   4.3472e-20
  -0.00094763
     -0.01504
   -0.0030592
  -4.8274e-19
    0.0015382
    -0.013476
   3.1519e-19
  -0.00029636
    -0.068934
   1.4502e-19
  -2.6821e-20
   3.1677e-20
    0.0065258
     0.060139
    -0.060139
   1.4396e-20
  -1.0014e-19
   7.2845e-20
      0.01149
      0.10596
    -0.037489
  -1.9587e-19
  -7.0128e-20
   3.0338e-20
      0.01383
      0.12669
   -0.0082612
   5.7633e-20
  -5.9712e-20
   5.8188e-20
     0.013313
      0.12025
      0.01835
   4.1473e-19
  -3.7219e-20
   2.2691e-20
     0.010567
     0.092894
     0.034807
  -6.2914e-20
  -5.6882e-21
   1.9487e-20
    0.0067929
     0.056436
     0.037794
  -2.2528e-19
  -1.8195e-20
  -8.0842e-21
    0.0032849
     0.023553
     0.029216
   1.7909e-19
  -5.4308e-20
   2.1127e-20
   0.00095163
    0.0031069
     0.015056
  -1.4947e-19
  -6.1059e-20
    1.224e-20
   1.0145e-06
   -0.0030639
     0.003085
  -2.3766e-19
   0.00082061
   -0.0021813
   5.8727e-19
   0.00021479
  -0.00021479
    7.003e-20
     0.068934
   0.00029636
   8.1622e-20
      0.12153
   0.00056145
   3.3927e-19
      0.14671
    0.0006561
   2.1676e-19
      0.14204
   0.00028095
  -1.4453e-19
      0.11403
  -0.00035334
   2.6474e-19
     0.075054
   -0.0010403
  -2.6582e-19
     0.038414
   -0.0014667
   1.9108e-19
     0.013505
   -0.0015381
   2.0892e-19
    0.0021986
  -0.00082087
   1.6469e-19
   0.00086881
  -0.00086981
  -6.6662e-19

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
    12          17          28.3356              1           4.34  
    13          18          28.2852              1           3.84  
    14          19          28.2487              1           12.2  
    15          21          28.2288       0.296525           2.66  
    16          22          28.2144              1           1.17  
    17          23          28.2062              1          0.609  
    18          24          28.2028              1           3.53  
    19          25          28.1992              1           6.37  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    20          26          28.1827              1           1.57  
    21          27          28.1705              1           2.86  
    22          28          28.1636              1           2.04  
    23          29          28.1588              1          0.861  
    24          30          28.1518              1           1.82  
    25          31           28.142              1           3.43  
    26          32          28.1326              1           3.11  
    27          33          28.1277              1           1.23  
    28          34          28.1258              1          0.441  
    29          35          28.1245              1          0.963  
    30          36           28.123              1           1.26  
    31          37          28.1216              1           0.84  
    32          38          28.1207              1          0.347  
    33          39          28.1201              1          0.435  
    34          40          28.1193              1          0.809  
    35          41          28.1181              1          0.887  
    36          42          28.1167              1          0.674  
    37          43          28.1153              1          0.493  
    38          44          28.1138              1          0.623  
    39          45          28.1118              1          0.995  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    40          46          28.1088              1           1.13  
    41          47          28.1047              1           1.52  
    42          48          28.1006              1           1.26  
    43          49          28.0978              1          0.532  
    44          50          28.0963              1          0.479  
    45          51          28.0952              1          0.621  
    46          52          28.0937              1          0.825  
    47          53           28.092              1          0.663  
    48          54          28.0906              1          0.607  
    49          55          28.0896              1          0.395  
    50          56          28.0885              1          0.539  
    51          57          28.0872              1          0.684  
    52          58          28.0857              1          0.868  
    53          59          28.0844              1          0.581  
    54          60          28.0834              1          0.457  
    55          61          28.0826              1          0.427  
    56          62          28.0817              1          0.703  
    57          63          28.0809              1          0.573  
    58          64          28.0804              1          0.207  
    59          65            28.08              1          0.263  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    60          66          28.0792              1          0.766  
    61          67           28.078              1           1.21  
    62          68          28.0763              1           1.22  
    63          69           28.075              1          0.641  
    64          70          28.0743              1          0.446  
    65          71          28.0737              1          0.548  
    66          72          28.0727              1              1  
    67          73          28.0713              1           1.17  
    68          74          28.0698              1          0.798  
    69          75          28.0688              1           0.28  
    70          76          28.0683              1          0.294  
    71          77           28.068              1          0.433  
    72          78          28.0677              1          0.309  
    73          79          28.0675              1          0.123  
    74          80          28.0674              1          0.146  
    75          81          28.0673              1          0.295  
    76          82           28.067              1          0.368  
    77          83          28.0666              1          0.353  
    78          84          28.0661              1          0.377  
    79          85          28.0655              1          0.258  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    80          86           28.065              1           0.37  
    81          87          28.0646              1          0.313  
    82          88          28.0642              1          0.291  
    83          89          28.0639              1          0.221  
    84          90          28.0636              1          0.305  
    85          91          28.0631              1          0.351  
    86          92          28.0625              1          0.494  
    87          93          28.0615              1          0.702  
    88          94          28.0602              1           0.67  
    89          95           28.059              1          0.505  
    90          96          28.0581              1          0.353  
    91          97          28.0573              1          0.405  
    92          98          28.0564              1          0.589  
    93          99          28.0553              1           0.59  
    94         100          28.0546              1           0.47  
    95         101          28.0543              1          0.158  
    96         102          28.0542              1          0.129  
    97         103           28.054              1          0.339  
    98         104          28.0536              1           0.56  
    99         105          28.0532              1          0.571  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   100         106          28.0528              1            0.3  
   101         107          28.0526              1          0.121  
   102         108          28.0525              1          0.256  
   103         109          28.0522              1          0.506  
   104         110          28.0517              1          0.667  
   105         111          28.0511              1          0.561  
   106         112          28.0507              1          0.241  
   107         113          28.0503              1          0.193  
   108         114            28.05              1          0.448  
   109         115          28.0494              1           0.66  
   110         116          28.0485              1           0.66  
   111         117          28.0473              1          0.574  
   112         118          28.0461              1          0.668  
   113         119          28.0446              1          0.669  
   114         120          28.0429              1          0.921  
   115         121          28.0414              1          0.703  
   116         122          28.0405              1          0.231  
   117         123            28.04              1          0.244  
   118         124          28.0392              1          0.591  
   119         125          28.0378              1          0.983  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   120         126          28.0351              1           1.23  
   121         127          28.0314              1           1.12  
   122         128          28.0283              1          0.686  
   123         129          28.0264              1          0.472  
   124         130          28.0252              1          0.484  
   125         131          28.0241              1          0.582  
   126         132          28.0231              1          0.489  
   127         133          28.0225              1          0.223  
   128         134          28.0221              1          0.224  
   129         135          28.0216              1          0.408  
   130         136          28.0208              1          0.686  
   131         137          28.0196              1          0.756  
   132         138          28.0184              1          0.479  
   133         139          28.0174              1          0.318  
   134         140          28.0162              1          0.524  
   135         141          28.0141              1           1.16  
   136         142          28.0101              1           1.79  
   137         143          28.0043              1           1.97  
   138         144          27.9989              1           1.27  
   139         145          27.9959              1           0.35  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   140         146          27.9943              1          0.454  
   141         147          27.9921              1           1.05  
   142         148          27.9879              1           1.62  
   143         149          27.9812              1           1.79  
   144         150          27.9738              1            1.2  
   145         151           27.968              1          0.853  
   146         152          27.9632              1          0.807  
   147         153          27.9563              1           1.61  
   148         154          27.9442              1           2.28  
   149         155          27.9272              1           2.44  
   150         156          27.9119              1           1.67  
   151         157          27.9033              1          0.889  
   152         158          27.8975              1          0.826  
   153         159          27.8892              1           1.69  
   154         160          27.8754              1           2.47  
   155         161          27.8584              1            2.3  
   156         162           27.845              1           1.05  
   157         163          27.8365              1          0.903  
   158         164          27.8283              1           1.38  
   159         165          27.8157              1           2.23  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   160         166          27.7988              1            2.3  
   161         167          27.7829              1           1.65  
   162         168          27.7727              1           1.46  
   163         169          27.7654              1           1.11  
   164         170          27.7575              1           1.38  
   165         171          27.7484              1           1.67  
   166         172          27.7403              1            1.2  
   167         173          27.7336              1           1.19  
   168         174           27.726              1           1.19  
   169         175          27.7138              1           2.24  
   170         176          27.6948              1           3.27  
   171         177          27.6728              1           3.09  
   172         178          27.6562              1           1.55  
   173         179          27.6469              1          0.936  
   174         180          27.6399              1           1.26  
   175         181          27.6304              1           2.09  
   176         182          27.6183              1           2.19  
   177         183          27.6075              1           1.29  
   178         184          27.6012              1            0.7  
   179         185          27.5973              1          0.635  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   180         186          27.5923              1           1.19  
   181         187          27.5832              1           1.64  
   182         188          27.5694              1           2.32  
   183         189          27.5549              1           1.99  
   184         190           27.545              1           1.07  
   185         191          27.5383              1            1.1  
   186         192            27.53              1           1.62  
   187         193          27.5166              1           2.59  
   188         194          27.4997              1            2.6  
   189         195          27.4859              1           1.35  
   190         196          27.4782              1          0.959  
   191         197          27.4731              1           0.99  
   192         198          27.4679              1           1.33  
   193         199          27.4623              1          0.999  
   194         200          27.4562              1           1.29  
   195         201          27.4477              1           1.45  
   196         202          27.4329              1           2.51  
   197         203          27.4077              1            3.7  
   198         204          27.3729              1           3.67  
   199         205          27.3363              1            1.9  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   200         206          27.3004              1           2.16  
   201         207          27.2555              1           3.63  
   202         208          27.1903              1           5.94  
   203         209           27.117              1           5.82  
   204         210          27.0698              1           2.89  
   205         211          27.0518              1          0.886  
   206         212          27.0417              1            1.4  
   207         213          27.0258              1           2.75  
   208         214          27.0002              1           3.46  
   209         215           26.967              1           2.75  
   210         216          26.9337              1           3.11  
   211         217          26.9022              1            2.6  
   212         218          26.8708              1           3.05  
   213         219          26.8431              1           2.84  
   214         220          26.8251              1           1.43  
   215         221           26.813              1            1.7  
   216         222          26.7976              1           2.25  
   217         223          26.7696              1           4.05  
   218         224          26.7266              1            4.9  
   219         225          26.6825              1           3.56  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   220         226          26.6535              1           2.18  
   221         227          26.6327              1            2.1  
   222         228          26.6024              1           3.77  
   223         229          26.5461              1           6.04  
   224         230          26.4653              1           6.62  
   225         231           26.397              1           4.11  
   226         232          26.3666              1           1.45  
   227         233          26.3531              1            1.4  
   228         234          26.3351              1           2.68  
   229         235          26.2991              1           4.31  
   230         236          26.2423              1           4.93  
   231         237           26.182              1           3.41  
   232         238          26.1434              1           1.99  
   233         239          26.1224              1           1.52  
   234         240          26.1042              1           2.09  
   235         241          26.0845              1           1.98  
   236         242          26.0684              1           1.24  
   237         243          26.0562              1           1.07  
   238         244          26.0423              1           1.59  
   239         245          26.0186              1           2.58  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   240         246          25.9824              1           3.34  
   241         247          25.9447              1            3.5  
   242         248          25.9232              1           2.06  
   243         249          25.9152              1          0.804  
   244         250          25.9102              1          0.785  
   245         251          25.9022              1           1.89  
   246         252          25.8875              1           3.31  
   247         253          25.8598              1           4.55  
   248         254          25.8144              1            4.8  
   249         255          25.7575              1           3.19  
   250         256          25.7094              1           2.64  
   251         257          25.6807              1            1.8  
   252         258          25.6637              1           2.09  
   253         259          25.6511              1           1.39  
   254         260          25.6419              1          0.816  
   255         261          25.6347              1          0.985  
   256         262          25.6258              1           1.61  
   257         263          25.6102              1           2.13  
   258         264          25.5845              1           2.58  
   259         265          25.5511              1           2.42  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   260         266          25.5186              1           1.38  
   261         267          25.4877              1           1.74  
   262         268          25.4485              1           2.63  
   263         269          25.3905              1           3.97  
   264         270          25.3231              1           3.79  
   265         271          25.2735              1           2.07  
   266         272          25.2464              1           1.75  
   267         273          25.2244              1           1.74  
   268         274          25.1933              1           2.72  
   269         275          25.1525              1           2.62  
   270         276           25.108              1           2.73  
   271         277          25.0602              1           3.36  
   272         278           24.997              1           3.37  
   273         279          24.9066              1           4.68  
   274         280          24.8038              1           4.29  
   275         281          24.7246              1           2.68  
   276         282          24.6751              1           2.59  
   277         283          24.6334              1           2.43  
   278         284          24.5867              1           3.02  
   279         285          24.5445              1           2.01  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   280         286          24.5145              1           1.55  
   281         287          24.4896              1           1.36  
   282         288           24.459              1           2.46  
   283         289          24.4225              1           2.49  
   284         290          24.3884              1           1.57  
   285         291          24.3583              1           2.11  
   286         292          24.3225              1           2.22  
   287         293           24.273              1           3.36  
   288         294          24.2223              1           2.99  
   289         295          24.1931              1           1.23  
   290         296          24.1818              1          0.979  
   291         297          24.1738              1          0.965  
   292         298          24.1611              1           1.59  
   293         299          24.1453              1           1.61  
   294         300          24.1328              1            1.3  
   295         301           24.126              1          0.499  
   296         302          24.1207              1          0.658  
   297         303          24.1111              1           1.55  
   298         304          24.0905              1           3.29  
   299         305          24.0492              1           5.39  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   300         306          23.9802              1           6.78  
   301         307          23.8951              1           5.77  
   302         308          23.8188              1           3.77  
   303         309          23.7518              1           4.01  
   304         310          23.6643              1           6.23  
   305         311          23.5156              1           10.5  
   306         312          23.2971              1           12.3  
   307         313          23.0895              1            8.8  
   308         314          22.9862              1           2.62  
   309         315          22.9514              1           1.22  
   310         316          22.9273              1           2.71  
   311         317          22.8885              1           3.85  
   312         318          22.8391              1           3.53  
   313         319          22.7969              1           1.78  
   314         320          22.7744              1           1.16  
   315         321          22.7597              1           1.32  
   316         322          22.7371              1           1.99  
   317         323          22.6901              1           3.22  
   318         324          22.5999              1            4.8  
   319         325          22.4641              1           5.15  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   320         326          22.3268              1           3.23  
   321         327            22.23              1           2.87  
   322         328          22.1553              1           2.62  
   323         329          22.0614              1           3.89  
   324         330          21.9309              1           4.35  
   325         331          21.7885              1           3.91  
   326         332          21.6631              1           4.78  
   327         333          21.5403              1           4.89  
   328         334          21.3799              1           4.71  
   329         335          21.1788              1           6.39  
   330         336          21.0069              1           5.56  
   331         337           20.915              1           2.49  
   332         338          20.8664              1           1.83  
   333         339          20.8162              1            2.9  
   334         340          20.7481              1           4.91  
   335         341          20.6759              1           4.94  
   336         342          20.6194              1           2.73  
   337         343          20.5785              1            2.8  
   338         344          20.5356              1           2.94  
   339         345           20.471              1           5.57  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   340         346          20.3723              1           7.14  
   341         347          20.2505              1           6.12  
   342         348          20.1463              1           3.59  
   343         349          20.0772              1           3.87  
   344         350          20.0187              1           3.47  
   345         351          19.9381              1           4.97  
   346         352          19.8285              1            4.8  
   347         353           19.727              1           3.02  
   348         354          19.6699              1           1.81  
   349         355          19.6378              1           1.66  
   350         356          19.5962              1           2.23  
   351         357          19.5156              1            4.1  
   352         358          19.3886              1           5.15  
   353         359          19.2605              1           3.93  
   354         360          19.1918              1           1.77  
   355         361          19.1633              1           1.28  
   356         362          19.1362              1           1.77  
   357         363          19.0845              1           3.02  
   358         364          19.0034              1           3.51  
   359         365          18.9063              1            4.8  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   360         366          18.8194              1           5.48  
   361         367          18.7405              1           4.75  
   362         368          18.6544              1           3.08  
   363         369           18.571              1            2.9  
   364         370          18.5202              1           1.53  
   365         371          18.4967              1           1.58  
   366         372          18.4757              1           1.41  
   367         373          18.4351              1           2.56  
   368         374          18.3619              1           3.53  
   369         375          18.2613              1           3.32  
   370         376          18.1733              1           1.65  
   371         377          18.1154              1            1.9  
   372         378          18.0636              1           1.86  
   373         379          17.9878              1           2.96  
   374         380          17.8919              1           3.28  
   375         381           17.818              1            2.5  
   376         382          17.7853              1           1.06  
   377         383          17.7683              1           1.06  
   378         384          17.7449              1           1.66  
   379         385           17.706              1           3.05  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   380         386          17.6576              1           3.49  
   381         387          17.6194              1           2.32  
   382         388           17.597              1           1.08  
   383         389          17.5785              1           1.11  
   384         390          17.5507              1            2.8  
   385         391          17.5094              1           3.95  
   386         392          17.4658              1           3.44  
   387         393          17.4389              1           1.48  
   388         394          17.4274              1          0.504  
   389         395          17.4193              1           1.14  
   390         396          17.4067              1           1.93  
   391         397          17.3865              1           2.26  
   392         398          17.3616              1           1.96  
   393         399          17.3412              1          0.856  
   394         400          17.3273              1           1.04  
   395         401          17.3127              1           1.41  
   396         402          17.2866              1           2.27  
   397         403          17.2409              1           3.01  
   398         404          17.1813              1           2.75  
   399         405          17.1339              1           1.31  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   400         406          17.1071              1           1.16  
   401         407          17.0835              1           1.33  
   402         408          17.0396              1           2.67  
   403         409          16.9546              1           4.05  
   404         410          16.8226              1            4.6  
   405         411          16.6889              1           3.27  
   406         412          16.5995              1           3.37  
   407         413          16.5314              1            3.4  
   408         414          16.4268              1           3.84  
   409         415           16.225              1           6.71  
   410         416          15.9097              1           8.46  
   411         417          15.6046              1           6.73  
   412         418          15.4559              1           2.69  
   413         419          15.4083              1           2.08  
   414         420          15.3761              1           1.61  
   415         421          15.3249              1           3.11  
   416         422          15.2651              1           3.45  
   417         423          15.2178              1           2.16  
   418         424          15.1909              1           2.05  
   419         425          15.1689              1           1.81  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   420         426          15.1359              1           2.79  
   421         427          15.0825              1           4.19  
   422         428          15.0142              1           4.53  
   423         429          14.9571              1           2.86  
   424         430          14.9246              1           2.02  
   425         431          14.9027              1           1.75  
   426         432          14.8745              1           2.65  
   427         433          14.8326              1           3.45  
   428         434          14.7862              1           2.81  
   429         435          14.7533              1           1.22  
   430         436          14.7332              1           1.32  
   431         437          14.7106              1           1.82  
   432         438          14.6653              1           3.27  
   433         439          14.5761              1           4.59  
   434         440          14.4349              1            4.8  
   435         441          14.2911              1           2.98  
   436         442          14.2002              1           2.01  
   437         443          14.1417              1           2.16  
   438         444          14.0627              1           3.16  
   439         445          13.9179              1           4.53  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   440         446          13.7143              1           4.55  
   441         447          13.5452              1           2.91  
   442         448          13.4704              1           1.79  
   443         449          13.4368              1           1.77  
   444         450          13.3952              1           1.87  
   445         451          13.3264              1            3.2  
   446         452           13.253              1           3.33  
   447         453          13.2129              1           1.82  
   448         454           13.201              1          0.588  
   449         455          13.1945              1          0.578  
   450         456          13.1807              1           1.94  
   451         457          13.1508              1           3.91  
   452         458            13.09              1           6.23  
   453         459          12.9992              1           7.25  
   454         460          12.9196              1           5.12  
   455         461           12.887              1           1.66  
   456         462          12.8785              1          0.671  
   457         463          12.8726              1           1.22  
   458         464          12.8592              1           2.64  
   459         465          12.8341              1           4.04  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   460         466          12.7952              1           4.55  
   461         467          12.7608              1           3.08  
   462         468          12.7459              1          0.871  
   463         469          12.7418              1          0.537  
   464         470          12.7389              1          0.795  
   465         471          12.7317              1           1.57  
   466         472          12.7159              1           2.49  
   467         473          12.6826              1            3.4  
   468         474          12.6305              1           3.53  
   469         475          12.5808              1            2.2  
   470         476          12.5566              1            1.3  
   471         477          12.5472              1           1.04  
   472         478          12.5368              1           1.11  
   473         479          12.5116              1           1.98  
   474         480          12.4553              1              3  
   475         481          12.3458              1            3.8  
   476         482          12.2012              1           3.41  
   477         483          12.0996              1           1.89  
   478         484          12.0663              1           1.27  
   479         485          12.0551              1          0.893  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   480         486          12.0387              1          0.842  
   481         487          12.0011              1            1.4  
   482         488          11.9308              1           1.84  
   483         489          11.8401              1           1.66  
   484         490          11.7805              1           1.06  
   485         491          11.7626              1           1.09  
   486         492           11.757              1            1.1  
   487         493          11.7483              1            1.1  
   488         494          11.7258              1           1.09  
   489         495          11.6722              1           1.37  
   490         496          11.5554              1           2.06  
   491         497          11.3605              1           2.63  
   492         498           11.161              1           2.19  
   493         499          11.0627              1           1.44  
   494         500          11.0345              1           1.76  
   495         501          11.0158              1           1.84  
   496         502          10.9714              1           1.85  
   497         503          10.8809              1           1.66  
   498         504          10.7282              1           1.83  
   499         505          10.5783              1           1.48  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   500         506           10.511              1          0.585  
   501         507          10.4972              1          0.576  
   502         508           10.492              1          0.554  
   503         509          10.4797              1          0.519  
   504         510          10.4519              1          0.673  
   505         511          10.3922              1           1.05  
   506         512           10.299              1           1.21  
   507         513          10.2141              1           0.84  
   508         514          10.1796              1          0.486  
   509         515          10.1721              1          0.493  
   510         516          10.1677              1          0.486  
   511         517          10.1559              1          0.461  
   512         518          10.1296              1          0.561  
   513         519          10.0721              1          0.839  
   514         520          9.98161              1          0.944  
   515         521          9.89743              1          0.643  
   516         522          9.86277              1          0.479  
   517         523          9.85558              1          0.449  
   518         524          9.85195              1           0.42  
   519         525          9.84229              1          0.418  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   520         526          9.82009              1          0.415  
   521         527          9.76832              1          0.628  
   522         528           9.6731              1          0.792  
   523         529          9.55531              1          0.819  
   524         530          9.48191              1           1.01  
   525         531          9.46008              1          0.962  
   526         532          9.45183              1          0.888  
   527         533          9.43691              1          0.764  
   528         534          9.40149              1           0.62  
   529         535          9.32555              1          0.563  
   530         536          9.20064              1          0.686  
   531         537          9.07865              1          0.886  
   532         538          9.02466              1          0.875  
   533         539          9.01334              1          0.741  
   534         540          9.00864              1          0.681  
   535         541           8.9975              1          0.555  
   536         542          8.97432              1          0.552  
   537         543          8.93227              1          0.479  
   538         544          8.88559              1          0.468  
   539         545          8.86031              1          0.456  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   540         546          8.85407              1          0.375  
   541         547          8.85181              1           0.38  
   542         548          8.84718              1          0.389  
   543         549          8.83586              1          0.401  
   544         550          8.80876              1          0.415  
   545         551           8.7528              1          0.673  
   546         552          8.66754              1          0.926  
   547         553          8.59419              1          0.932  
   548         554          8.56584              1           1.01  
   549         555          8.55913              1           1.13  
   550         556          8.55407              1           1.14  
   551         557            8.541              1            1.1  
   552         558          8.51497              1          0.907  
   553         559          8.47214              1          0.446  
   554         560          8.43247              1           0.52  
   555         561          8.41575              1          0.477  
   556         562          8.41229              1            0.5  
   557         563           8.4106              1          0.473  
   558         564          8.40618              1          0.406  
   559         565          8.39605              1            0.3  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   560         566          8.37272              1          0.409  
   561         567          8.33099              1          0.616  
   562         568          8.28241              1          0.589  
   563         569          8.25542              1          0.446  
   564         570          8.24908              1           0.34  
   565         571          8.24752              1          0.281  
   566         572          8.24505              1          0.236  
   567         573          8.23883              1          0.206  
   568         574          8.22508              1          0.318  
   569         575          8.19992              1          0.432  
   570         576          8.17037              1          0.523  
   571         577          8.15317              1          0.465  
   572         578          8.14853              1          0.436  
   573         579          8.14686              1          0.413  
   574         580          8.14369              1           0.37  
   575         581          8.13609              1          0.288  
   576         582           8.1192              1          0.303  
   577         583           8.0894              1           0.45  
   578         584           8.0563              1          0.525  
   579         585          8.03892              1           0.44  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   580         586          8.03496              1          0.401  
   581         587          8.03381              1          0.363  
   582         588          8.03166              1          0.304  
   583         589          8.02645              1          0.204  
   584         590          8.01489              1          0.243  
   585         591          7.99439              1          0.415  
   586         592          7.97133              1          0.568  
   587         593          7.95874              1          0.605  
   588         594           7.9554              1          0.518  
   589         595          7.95397              1          0.452  
   590         596          7.95087              1          0.345  
   591         597          7.94337              1          0.255  
   592         598           7.9253              1          0.367  
   593         599          7.88799              1          0.605  
   594         600          7.83095              1           1.08  
   595         601          7.78176              1           1.22  
   596         602          7.76279              1          0.987  
   597         603          7.75845              1          0.788  
   598         604          7.75535              1          0.655  
   599         605          7.74717              1           0.48  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   600         606          7.73025              1           0.45  
   601         607          7.70003              1            0.4  
   602         608          7.66775              1          0.616  
   603         609          7.65154              1          0.453  
   604         610          7.64822              1            0.3  
   605         611          7.64742              1          0.301  
   606         612            7.646              1            0.3  
   607         613          7.64232              1          0.296  
   608         614          7.63328              1          0.366  
   609         615           7.6126              1           0.73  
   610         616           7.5739              1           1.11  
   611         617          7.52546              1           1.19  
   612         618           7.4944              1          0.831  
   613         619          7.48483              1           0.73  
   614         620           7.4811              1          0.678  
   615         621          7.47445              1          0.594  
   616         622          7.45875              1          0.528  
   617         623           7.4257              1          0.642  
   618         624          7.37352              1          0.897  
   619         625           7.3265              1          0.708  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   620         626          7.30845              1          0.264  
   621         627          7.30568              1          0.189  
   622         628          7.30503              1          0.186  
   623         629          7.30355              1          0.181  
   624         630          7.29991              1          0.195  
   625         631          7.29071              1          0.363  
   626         632          7.26969              1          0.582  
   627         633          7.22956              1          0.768  
   628         634          7.17814              1           0.71  
   629         635          7.14479              1          0.419  
   630         636           7.1357              1          0.443  
   631         637          7.13372              1          0.429  
   632         638          7.13142              1          0.409  
   633         639           7.1251              1           0.36  
   634         640          7.11072              1          0.333  
   635         641          7.08047              1          0.483  
   636         642          7.03503              1          0.519  
   637         643          6.99648              1          0.449  
   638         644          6.98237              1          0.461  
   639         645           6.9796              1          0.414  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   640         646          6.97796              1           0.38  
   641         647          6.97329              1          0.311  
   642         648          6.96259              1          0.238  
   643         649          6.93739              1          0.338  
   644         650          6.89061              1          0.434  
   645         651           6.8321              1          0.386  
   646         652          6.79629              1          0.289  
   647         653          6.78754              1          0.292  
   648         654          6.78615              1          0.311  
   649         655          6.78489              1          0.315  
   650         656          6.78109              1          0.316  
   651         657          6.77244              1          0.305  
   652         658          6.75276              1           0.26  
   653         659          6.71875              1            0.3  
   654         660          6.68169              1          0.306  
   655         661          6.66328              1           0.31  
   656         662          6.65972              1          0.271  
   657         663          6.65912              1          0.251  
   658         664          6.65832              1          0.225  
   659         665          6.65618              1          0.167  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   660         666           6.6516              1          0.112  
   661         667          6.64326              1           0.14  
   662         668           6.6337              1          0.348  
   663         669          6.62817              1          0.432  
   664         670          6.62652              1          0.405  
   665         671          6.62567              1          0.369  
   666         672           6.6238              1          0.308  
   667         673           6.6192              1          0.207  
   668         674          6.60768              1          0.179  
   669         675          6.58158              1          0.243  
   670         676          6.53276              1          0.512  
   671         677          6.47245              1          0.661  
   672         678          6.43521              1          0.509  
   673         679          6.42524              1           0.53  
   674         680          6.42257              1          0.534  
   675         681           6.4186              1          0.522  
   676         682          6.40831              1           0.48  
   677         683          6.38458              1          0.369  
   678         684          6.33627              1          0.538  
   679         685          6.26681              1          0.669  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   680         686          6.21248              1          0.629  
   681         687          6.19441              1          0.704  
   682         688          6.19093              1          0.663  
   683         689          6.18848              1          0.621  
   684         690          6.18147              1          0.527  
   685         691          6.16614              1          0.372  
   686         692          6.13408              1          0.419  
   687         693          6.08807              1          0.481  
   688         694          6.05182              1          0.322  
   689         695          6.04027              1          0.185  
   690         696          6.03869              1          0.166  
   691         697          6.03815              1          0.167  
   692         698          6.03654              1          0.168  
   693         699          6.03285              1          0.166  
   694         700          6.02346              1          0.163  
   695         701           6.0031              1          0.313  
   696         702          5.96767              1          0.425  
   697         703           5.9301              1          0.363  
   698         704          5.91192              1          0.224  
   699         705          5.90854              1          0.228  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   700         706          5.90799              1           0.22  
   701         707          5.90723              1           0.21  
   702         708          5.90508              1          0.191  
   703         709          5.89978              1          0.158  
   704         710          5.88643              1          0.213  
   705         711          5.85615              1           0.33  
   706         712          5.79937              1          0.421  
   707         713          5.72939              1           0.37  
   708         714          5.68733              1           0.22  
   709         715          5.67774              1          0.198  
   710         716          5.67668              1          0.195  
   711         717          5.67611              1          0.191  
   712         718          5.67418              1          0.178  
   713         719          5.66994              1          0.155  
   714         720          5.65974              1          0.138  
   715         721          5.64088              1          0.126  
   716         722          5.61706              1           0.17  
   717         723          5.60234              1           0.25  
   718         724           5.5987              1          0.244  
   719         725          5.59811              1          0.227  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   720         726          5.59759              1          0.212  
   721         727          5.59599              1          0.183  
   722         728          5.59224              1          0.153  
   723         729          5.58304              1          0.156  
   724         730          5.56417              1          0.152  
   725         731          5.53534              1          0.197  
   726         732          5.51102              1          0.205  
   727         733          5.50244              1          0.198  
   728         734          5.50122              1          0.212  
   729         735          5.50089              1           0.21  
   730         736          5.50006              1          0.205  
   731         737          5.49807              1          0.196  
   732         738          5.49301              1          0.176  
   733         739          5.48152              1          0.136  
   734         740          5.45973              1          0.219  
   735         741          5.43237              1          0.253  
   736         742          5.41531              1          0.193  
   737         743          5.41111              1          0.174  
   738         744          5.41048              1          0.176  
   739         745          5.40998              1          0.175  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   740         746           5.4084              1          0.173  
   741         747          5.40474              1           0.17  
   742         748          5.39577              1          0.162  
   743         749          5.37758              1          0.148  
   744         750          5.35038              1          0.221  
   745         751          5.32823              1          0.207  
   746         752           5.3207              1          0.182  
   747         753          5.31957              1          0.195  
   748         754          5.31915              1          0.196  
   749         755          5.31794              1          0.196  
   750         756          5.31512              1          0.192  
   751         757          5.30794              1          0.179  
   752         758          5.29213              1          0.145  
   753         759           5.2638              1          0.207  
   754         760          5.23206              1          0.219  
   755         761          5.21548              1          0.154  
   756         762          5.21223              1          0.158  
   757         763          5.21182              1          0.157  
   758         764          5.21146              1          0.155  
   759         765           5.2103              1          0.151  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   760         766           5.2076              1          0.145  
   761         767          5.20094              1          0.135  
   762         768           5.1872              1          0.156  
   763         769          5.16584              1          0.262  
   764         770           5.1472              1          0.302  
   765         771          5.14014              1          0.261  
   766         772          5.13886              1          0.223  
   767         773          5.13827              1          0.205  
   768         774          5.13663              1          0.171  
   769         775          5.13281              1          0.145  
   770         776          5.12324              1          0.148  
   771         777          5.10286              1          0.146  
   772         778          5.06905              1          0.202  
   773         779          5.03616              1          0.215  
   774         780          5.02222              1          0.122  
   775         781          5.02007              1          0.122  
   776         782          5.01981              1          0.121  
   777         783           5.0195              1          0.119  
   778         784          5.01856              1          0.116  
   779         785          5.01634              1          0.109  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   780         786          5.01095              1         0.0963  
   781         787          4.99997              1          0.162  
   782         788           4.9836              1          0.219  
   783         789          4.97032              1          0.202  
   784         790          4.96587              1           0.14  
   785         791          4.96526              1          0.109  
   786         792          4.96508              1          0.101  
   787         793          4.96459              1         0.0985  
   788         794          4.96344              1          0.092  
   789         795          4.96053              1         0.0743  
   790         796          4.95417              1         0.0743  
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
   50  00:00:01   1.0560e+03      6          7.4668e+02   1.5593e-03   7.9396e-03   1.4749e+00   4.6352e+00   2.0000e-01
   60  00:00:01   1.0510e+03      7          7.4315e+02   1.2924e-03   6.8079e-03   1.4917e+00   4.6465e+00   2.0000e-01
   70  00:00:01   1.0460e+03      8          7.3961e+02   1.1032e-03   5.9575e-03   1.5037e+00   4.6550e+00   2.0000e-01
   80  00:00:01   1.0410e+03      9          7.3608e+02   9.6224e-04   5.2954e-03   1.5125e+00   4.6616e+00   2.0000e-01
   90  00:00:01   1.0360e+03     10          7.3254e+02   8.5312e-04   4.7655e-03   1.5194e+00   4.6669e+00   2.0000e-01
  100  00:00:01   1.0310e+03     11          7.2900e+02   7.6620e-04   4.3318e-03   1.5249e+00   4.6713e+00   2.0000e-01
  110  00:00:01   1.0260e+03     12          7.2547e+02   6.9532e-04   3.9703e-03   1.5293e+00   4.6749e+00   2.0000e-01
  111  00:00:01   1.0259e+03     13  EP      7.2541e+02   6.9430e-04   3.9651e-03   1.5294e+00   4.6749e+00   2.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:01   1.0799e+03     14  EP      7.6359e+02   1.4414e-02   3.3645e-02   1.8524e-01   4.3780e+00   2.0000e-01
   10  00:00:01   1.0837e+03     15          7.6626e+02   1.0511e-02   5.1933e-02  -9.8992e-01   4.1644e+00   2.0000e-01
   20  00:00:02   1.0886e+03     16          7.6974e+02   5.2264e-03   8.2415e-02  -1.3523e+00   3.7411e+00   2.0000e-01
   30  00:00:02   1.0924e+03     17          7.7246e+02   3.6455e-03   9.9255e-02  -1.1458e+00   3.0464e+00   2.0000e-01
   38  00:00:02   1.0925e+03     18  FP      7.7250e+02   3.7252e-03   9.7795e-02  -1.1258e+00   2.9454e+00   2.0000e-01
   38  00:00:02   1.0925e+03     19  SN      7.7250e+02   3.7252e-03   9.7795e-02  -1.1258e+00   2.9454e+00   2.0000e-01
   40  00:00:02   1.0925e+03     20          7.7250e+02   3.7513e-03   9.7222e-02  -1.1227e+00   2.9178e+00   2.0000e-01
   50  00:00:02   1.0923e+03     21          7.7237e+02   3.9150e-03   9.2331e-02  -1.1249e+00   2.7541e+00   2.0000e-01
   60  00:00:02   1.0913e+03     22          7.7165e+02   3.8887e-03   6.2201e-02  -1.2665e+00   2.2444e+00   2.0000e-01
   62  00:00:02   1.0913e+03     23  SN      7.7165e+02   3.8715e-03   6.1349e-02  -1.2703e+00   2.2336e+00   2.0000e-01
   62  00:00:02   1.0913e+03     24  FP      7.7165e+02   3.8715e-03   6.1349e-02  -1.2703e+00   2.2336e+00   2.0000e-01
   70  00:00:03   1.0913e+03     25          7.7167e+02   3.7719e-03   5.6954e-02  -1.2893e+00   2.1789e+00   2.0000e-01
   80  00:00:03   1.0917e+03     26          7.7196e+02   3.4473e-03   4.6120e-02  -1.3325e+00   2.0522e+00   2.0000e-01
   90  00:00:03   1.0962e+03     27          7.7513e+02   2.3209e-03   2.2918e-02  -1.4215e+00   1.8041e+00   2.0000e-01
  100  00:00:03   1.1012e+03     28          7.7866e+02   1.7620e-03   1.5366e-02  -1.4576e+00   1.7271e+00   2.0000e-01
  110  00:00:03   1.1062e+03     29          7.8220e+02   1.4250e-03   1.1600e-02  -1.4790e+00   1.6892e+00   2.0000e-01
  120  00:00:03   1.1112e+03     30          7.8573e+02   1.1972e-03   9.3226e-03  -1.4934e+00   1.6663e+00   2.0000e-01
  130  00:00:03   1.1162e+03     31          7.8927e+02   1.0324e-03   7.7940e-03  -1.5038e+00   1.6509e+00   2.0000e-01
  140  00:00:03   1.1212e+03     32          7.9280e+02   9.0757e-04   6.6964e-03  -1.5116e+00   1.6399e+00   2.0000e-01
  150  00:00:04   1.1262e+03     33          7.9634e+02   8.0969e-04   5.8698e-03  -1.5178e+00   1.6317e+00   2.0000e-01
  160  00:00:04   1.1312e+03     34          7.9988e+02   7.3087e-04   5.2249e-03  -1.5227e+00   1.6252e+00   2.0000e-01
  170  00:00:04   1.1362e+03     35          8.0341e+02   6.6604e-04   4.7076e-03  -1.5268e+00   1.6200e+00   2.0000e-01
  180  00:00:04   1.1412e+03     36          8.0695e+02   6.1177e-04   4.2835e-03  -1.5302e+00   1.6158e+00   2.0000e-01
  190  00:00:04   1.1462e+03     37          8.1048e+02   5.6568e-04   3.9295e-03  -1.5331e+00   1.6123e+00   2.0000e-01
  200  00:00:04   1.1512e+03     38          8.1402e+02   5.2605e-04   3.6295e-03  -1.5356e+00   1.6093e+00   2.0000e-01
  210  00:00:04   1.1562e+03     39          8.1755e+02   4.9161e-04   3.3721e-03  -1.5377e+00   1.6067e+00   2.0000e-01
  220  00:00:04   1.1612e+03     40          8.2109e+02   4.6140e-04   3.1487e-03  -1.5396e+00   1.6045e+00   2.0000e-01
  230  00:00:05   1.1662e+03     41          8.2462e+02   4.3468e-04   2.9531e-03  -1.5413e+00   1.6025e+00   2.0000e-01
  240  00:00:05   1.1712e+03     42          8.2816e+02   4.1089e-04   2.7804e-03  -1.5428e+00   1.6008e+00   2.0000e-01
  250  00:00:05   1.1762e+03     43          8.3170e+02   3.8957e-04   2.6268e-03  -1.5442e+00   1.5992e+00   2.0000e-01
  260  00:00:05   1.1812e+03     44          8.3523e+02   3.7035e-04   2.4892e-03  -1.5454e+00   1.5979e+00   2.0000e-01
  270  00:00:05   1.1862e+03     45          8.3877e+02   3.5294e-04   2.3654e-03  -1.5465e+00   1.5966e+00   2.0000e-01
  274  00:00:05   1.1879e+03     46  EP      8.3995e+02   3.4747e-04   2.3266e-03  -1.5468e+00   1.5962e+00   2.0000e-01
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
   50  00:00:02   1.1035e+03     18          7.8027e+02   5.9079e-03   1.7080e-01  -1.7789e+00   3.6597e+00   4.0000e-01
   60  00:00:02   1.1085e+03     19          7.8379e+02   3.9690e-03   1.8897e-01  -1.5871e+00   3.4050e+00   4.0000e-01
   70  00:00:02   1.1108e+03     20          7.8542e+02   3.7724e-03   1.9523e-01  -1.1608e+00   3.1086e+00   4.0000e-01
   74  00:00:02   1.1108e+03     21  SN      7.8542e+02   3.8231e-03   1.9509e-01  -1.1398e+00   3.0900e+00   4.0000e-01
   74  00:00:02   1.1108e+03     22  FP      7.8542e+02   3.8231e-03   1.9509e-01  -1.1398e+00   3.0900e+00   4.0000e-01
   80  00:00:02   1.1107e+03     23          7.8541e+02   3.9099e-03   1.9475e-01  -1.1105e+00   3.0619e+00   4.0000e-01
   90  00:00:02   1.1105e+03     24          7.8524e+02   4.2011e-03   1.9307e-01  -1.0470e+00   2.9843e+00   4.0000e-01
  100  00:00:02   1.1067e+03     25          7.8255e+02   5.8058e-03   1.7344e-01  -1.0032e+00   2.6513e+00   4.0000e-01
  110  00:00:02   1.1017e+03     26          7.7903e+02   6.7861e-03   1.4495e-01  -1.1109e+00   2.3933e+00   4.0000e-01
  120  00:00:03   1.0967e+03     27          7.7551e+02   6.8438e-03   1.0411e-01  -1.2635e+00   2.1225e+00   4.0000e-01
  129  00:00:03   1.0953e+03     28  FP      7.7450e+02   6.0074e-03   7.3141e-02  -1.3541e+00   1.9478e+00   4.0000e-01
  129  00:00:03   1.0953e+03     29  SN      7.7450e+02   6.0074e-03   7.3141e-02  -1.3541e+00   1.9478e+00   4.0000e-01
  130  00:00:03   1.0953e+03     30          7.7450e+02   5.9407e-03   7.1465e-02  -1.3583e+00   1.9388e+00   4.0000e-01
  140  00:00:03   1.0956e+03     31          7.7474e+02   5.4466e-03   6.0484e-02  -1.3852e+00   1.8804e+00   4.0000e-01
  150  00:00:03   1.1003e+03     32          7.7806e+02   3.7724e-03   3.3882e-02  -1.4491e+00   1.7432e+00   4.0000e-01
  160  00:00:03   1.1053e+03     33          7.8159e+02   2.9790e-03   2.4591e-02  -1.4748e+00   1.6962e+00   4.0000e-01
  170  00:00:03   1.1103e+03     34          7.8513e+02   2.4760e-03   1.9443e-02  -1.4908e+00   1.6703e+00   4.0000e-01
  180  00:00:04   1.1153e+03     35          7.8866e+02   2.1219e-03   1.6109e-02  -1.5020e+00   1.6536e+00   4.0000e-01
  190  00:00:04   1.1203e+03     36          7.9220e+02   1.8576e-03   1.3762e-02  -1.5103e+00   1.6418e+00   4.0000e-01
  200  00:00:04   1.1253e+03     37          7.9573e+02   1.6523e-03   1.2015e-02  -1.5167e+00   1.6331e+00   4.0000e-01
  210  00:00:04   1.1303e+03     38          7.9927e+02   1.4881e-03   1.0664e-02  -1.5219e+00   1.6263e+00   4.0000e-01
  220  00:00:04   1.1353e+03     39          8.0280e+02   1.3537e-03   9.5868e-03  -1.5261e+00   1.6209e+00   4.0000e-01
  230  00:00:04   1.1403e+03     40          8.0634e+02   1.2416e-03   8.7076e-03  -1.5296e+00   1.6165e+00   4.0000e-01
  240  00:00:04   1.1453e+03     41          8.0988e+02   1.1467e-03   7.9763e-03  -1.5326e+00   1.6128e+00   4.0000e-01
  250  00:00:04   1.1503e+03     42          8.1341e+02   1.0653e-03   7.3585e-03  -1.5352e+00   1.6098e+00   4.0000e-01
  260  00:00:05   1.1553e+03     43          8.1695e+02   9.9467e-04   6.8296e-03  -1.5374e+00   1.6071e+00   4.0000e-01
  270  00:00:05   1.1603e+03     44          8.2048e+02   9.3285e-04   6.3716e-03  -1.5393e+00   1.6048e+00   4.0000e-01
  280  00:00:05   1.1653e+03     45          8.2402e+02   8.7826e-04   5.9713e-03  -1.5410e+00   1.6028e+00   4.0000e-01
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
Estimated memory usage at order  2 = 3.09E+00 MB
Manifold computation time at order 3 = 00:00:00
Estimated memory usage at order  3 = 6.37E+00 MB
Manifold computation time at order 4 = 00:00:01
Estimated memory usage at order  4 = 9.22E+00 MB
Manifold computation time at order 5 = 00:00:04
Estimated memory usage at order  5 = 2.01E+01 MB

 Run='SSMToolvonKarmanPlateIReps.ep': Continue equilibria with varied epsilon.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          2.20e-01  7.64e+02    0.0    0.0    0.0
   1   1  1.00e+00  5.69e-02  1.39e-03  7.64e+02    0.0    0.0    0.0
   2   1  1.00e+00  4.30e-04  8.75e-08  7.64e+02    0.0    0.0    0.0
   3   1  1.00e+00  7.84e-08  1.78e-15  7.64e+02    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE           eps         rho1         rho2          th1          th2           om
    0  00:00:00   7.6367e+02      1  EP      2.0000e-01   2.9965e-03   6.9388e-03   4.9180e+00   5.9505e+00   7.6359e+02
    1  00:00:00   7.6367e+02      2  UZ      2.0000e-01   2.9965e-03   6.9388e-03   4.9180e+00   5.9505e+00   7.6359e+02
    1  00:00:00   7.6367e+02      3  EP      1.8000e-01   2.7578e-03   6.3613e-03   4.8855e+00   5.9430e+00   7.6359e+02

 STEP      TIME        ||U||  LABEL  TYPE           eps         rho1         rho2          th1          th2           om
    0  00:00:00   7.6367e+02      4  EP      2.0000e-01   2.9965e-03   6.9388e-03   4.9180e+00   5.9505e+00   7.6359e+02
    4  00:00:00   7.6368e+02      5  UZ      4.0000e-01   4.4921e-03   1.1521e-02   5.2281e+00   6.0162e+00   7.6359e+02
    5  00:00:00   7.6368e+02      6  EP      4.4000e-01   4.6550e-03   1.2247e-02   5.2798e+00   6.0262e+00   7.6359e+02

 Run='SSMToolvonKarmanPlateIReps1.ep': Continue equilibria with varied omega at eps equal to 2.000000e-01.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          1.78e-15  1.08e+03    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   1.0799e+03      1  EP      7.6359e+02   2.9965e-03   6.9388e-03   4.9180e+00   5.9505e+00   2.0000e-01
   10  00:00:00   1.0760e+03      2          7.6080e+02   1.5957e-03   4.7220e-03   5.7727e+00   6.0558e+00   2.0000e-01
   20  00:00:00   1.0710e+03      3          7.5728e+02   8.3992e-04   3.2666e-03   6.0289e+00   6.1262e+00   2.0000e-01
   30  00:00:00   1.0660e+03      4          7.5375e+02   5.5565e-04   2.4732e-03   6.1165e+00   6.1645e+00   2.0000e-01
   40  00:00:00   1.0610e+03      5          7.5021e+02   4.1304e-04   1.9849e-03   6.1597e+00   6.1880e+00   2.0000e-01
   50  00:00:00   1.0560e+03      6          7.4668e+02   3.2816e-04   1.6561e-03   6.1852e+00   6.2038e+00   2.0000e-01
   60  00:00:00   1.0510e+03      7          7.4314e+02   2.7205e-04   1.4202e-03   6.2020e+00   6.2151e+00   2.0000e-01
   70  00:00:00   1.0460e+03      8          7.3960e+02   2.3225e-04   1.2429e-03   6.2139e+00   6.2236e+00   2.0000e-01
   80  00:00:00   1.0410e+03      9          7.3607e+02   2.0258e-04   1.1048e-03   6.2228e+00   6.2303e+00   2.0000e-01
   90  00:00:01   1.0360e+03     10          7.3253e+02   1.7962e-04   9.9429e-04   6.2296e+00   6.2356e+00   2.0000e-01
  100  00:00:01   1.0310e+03     11          7.2900e+02   1.6132e-04   9.0382e-04   6.2351e+00   6.2399e+00   2.0000e-01
  110  00:00:01   1.0260e+03     12          7.2546e+02   1.4640e-04   8.2843e-04   6.2395e+00   6.2435e+00   2.0000e-01
  111  00:00:01   1.0260e+03     13  EP      7.2541e+02   1.4621e-04   8.2746e-04   6.2396e+00   6.2436e+00   2.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:01   1.0799e+03     14  EP      7.6359e+02   2.9965e-03   6.9388e-03   4.9180e+00   5.9505e+00   2.0000e-01
   10  00:00:01   1.0837e+03     15          7.6625e+02   2.3054e-03   1.0583e-02   3.7272e+00   5.7460e+00   2.0000e-01
   20  00:00:01   1.0886e+03     16          7.6974e+02   1.1753e-03   1.6622e-02   3.3194e+00   5.3572e+00   2.0000e-01
   30  00:00:01   1.0931e+03     17          7.7288e+02   7.2478e-04   2.0836e-02   3.5185e+00   4.6847e+00   2.0000e-01
   40  00:00:01   1.0932e+03     18  FP      7.7299e+02   7.4701e-04   2.0518e-02   3.5751e+00   4.5331e+00   2.0000e-01
   40  00:00:01   1.0932e+03     19  SN      7.7299e+02   7.4701e-04   2.0518e-02   3.5751e+00   4.5331e+00   2.0000e-01
   40  00:00:01   1.0932e+03     20          7.7299e+02   7.4876e-04   2.0489e-02   3.5770e+00   4.5252e+00   2.0000e-01
   50  00:00:02   1.0931e+03     21          7.7293e+02   7.7779e-04   1.9901e-02   3.5939e+00   4.4081e+00   2.0000e-01
   60  00:00:02   1.0917e+03     22          7.7191e+02   8.3009e-04   1.3865e-02   3.4713e+00   3.8677e+00   2.0000e-01
   68  00:00:02   1.0916e+03     23  SN      7.7184e+02   8.0109e-04   1.2254e-02   3.4347e+00   3.7687e+00   2.0000e-01
   68  00:00:02   1.0916e+03     24  FP      7.7184e+02   8.0109e-04   1.2254e-02   3.4347e+00   3.7687e+00   2.0000e-01
   70  00:00:02   1.0916e+03     25          7.7185e+02   7.9543e-04   1.2006e-02   3.4293e+00   3.7541e+00   2.0000e-01
   80  00:00:02   1.0917e+03     26          7.7191e+02   7.6482e-04   1.0842e-02   3.4051e+00   3.6874e+00   2.0000e-01
   90  00:00:02   1.0937e+03     27          7.7334e+02   5.9308e-04   6.5533e-03   3.3241e+00   3.4608e+00   2.0000e-01
  100  00:00:02   1.0987e+03     28          7.7687e+02   4.2281e-04   3.8519e-03   3.2685e+00   3.3271e+00   2.0000e-01
  110  00:00:02   1.1037e+03     29          7.8041e+02   3.3244e-04   2.7655e-03   3.2410e+00   3.2744e+00   2.0000e-01
  120  00:00:02   1.1087e+03     30          7.8394e+02   2.7442e-04   2.1610e-03   3.2235e+00   3.2452e+00   2.0000e-01
  130  00:00:03   1.1137e+03     31          7.8748e+02   2.3375e-04   1.7742e-03   3.2114e+00   3.2266e+00   2.0000e-01
  140  00:00:03   1.1187e+03     32          7.9101e+02   2.0362e-04   1.5050e-03   3.2023e+00   3.2137e+00   2.0000e-01
  150  00:00:03   1.1237e+03     33          7.9455e+02   1.8038e-04   1.3068e-03   3.1954e+00   3.2042e+00   2.0000e-01
  160  00:00:03   1.1287e+03     34          7.9809e+02   1.6191e-04   1.1547e-03   3.1899e+00   3.1969e+00   2.0000e-01
  170  00:00:03   1.1337e+03     35          8.0162e+02   1.4687e-04   1.0343e-03   3.1854e+00   3.1911e+00   2.0000e-01
  180  00:00:03   1.1387e+03     36          8.0516e+02   1.3439e-04   9.3669e-04   3.1817e+00   3.1865e+00   2.0000e-01
  190  00:00:03   1.1437e+03     37          8.0869e+02   1.2386e-04   8.5589e-04   3.1785e+00   3.1826e+00   2.0000e-01
  200  00:00:03   1.1487e+03     38          8.1223e+02   1.1486e-04   7.8793e-04   3.1758e+00   3.1793e+00   2.0000e-01
  210  00:00:03   1.1537e+03     39          8.1576e+02   1.0708e-04   7.2996e-04   3.1735e+00   3.1765e+00   2.0000e-01
  220  00:00:04   1.1587e+03     40          8.1930e+02   1.0029e-04   6.7993e-04   3.1715e+00   3.1742e+00   2.0000e-01
  230  00:00:04   1.1637e+03     41          8.2283e+02   9.4309e-05   6.3632e-04   3.1697e+00   3.1721e+00   2.0000e-01
  240  00:00:04   1.1687e+03     42          8.2637e+02   8.9000e-05   5.9797e-04   3.1681e+00   3.1702e+00   2.0000e-01
  250  00:00:04   1.1737e+03     43          8.2991e+02   8.4257e-05   5.6397e-04   3.1667e+00   3.1686e+00   2.0000e-01
  260  00:00:04   1.1787e+03     44          8.3344e+02   7.9994e-05   5.3364e-04   3.1654e+00   3.1671e+00   2.0000e-01
  270  00:00:04   1.1837e+03     45          8.3698e+02   7.6142e-05   5.0639e-04   3.1643e+00   3.1658e+00   2.0000e-01
  279  00:00:04   1.1879e+03     46  EP      8.3995e+02   7.3176e-05   4.8554e-04   3.1634e+00   3.1648e+00   2.0000e-01

 Run='SSMToolvonKarmanPlateIReps2.ep': Continue equilibria with varied omega at eps equal to 4.000000e-01.

    STEP   DAMPING               NORMS              COMPUTATION TIMES
  IT SIT     GAMMA     ||d||     ||f||     ||U||   F(x)  DF(x)  SOLVE
   0                          8.00e-13  1.08e+03    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:00   1.0799e+03      1  EP      7.6359e+02   4.4921e-03   1.1521e-02   5.2281e+00   6.0162e+00   4.0000e-01
   10  00:00:00   1.0759e+03      2          7.6073e+02   2.6890e-03   8.6112e-03   5.8365e+00   6.0770e+00   4.0000e-01
   20  00:00:00   1.0709e+03      3          7.5720e+02   1.5817e-03   6.2800e-03   6.0412e+00   6.1324e+00   4.0000e-01
   30  00:00:00   1.0659e+03      4          7.5367e+02   1.0811e-03   4.8469e-03   6.1205e+00   6.1669e+00   4.0000e-01
   40  00:00:00   1.0609e+03      5          7.5013e+02   8.1337e-04   3.9214e-03   6.1614e+00   6.1892e+00   4.0000e-01
   50  00:00:00   1.0559e+03      6          7.4660e+02   6.4976e-04   3.2849e-03   6.1861e+00   6.2045e+00   4.0000e-01
   60  00:00:00   1.0509e+03      7          7.4306e+02   5.4021e-04   2.8232e-03   6.2025e+00   6.2155e+00   4.0000e-01
   70  00:00:00   1.0459e+03      8          7.3953e+02   4.6198e-04   2.4740e-03   6.2142e+00   6.2239e+00   4.0000e-01
   80  00:00:00   1.0409e+03      9          7.3599e+02   4.0340e-04   2.2011e-03   6.2230e+00   6.2305e+00   4.0000e-01
   90  00:00:00   1.0359e+03     10          7.3245e+02   3.5794e-04   1.9822e-03   6.2298e+00   6.2357e+00   4.0000e-01
  100  00:00:00   1.0309e+03     11          7.2892e+02   3.2165e-04   1.8026e-03   6.2352e+00   6.2400e+00   4.0000e-01
  110  00:00:01   1.0260e+03     12  EP      7.2541e+02   2.9224e-04   1.6539e-03   6.2396e+00   6.2436e+00   4.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1         rho2          th1          th2          eps
    0  00:00:01   1.0799e+03     13  EP      7.6359e+02   4.4921e-03   1.1521e-02   5.2281e+00   6.0162e+00   4.0000e-01
   10  00:00:01   1.0837e+03     14          7.6625e+02   7.2276e-03   1.4894e-02   3.9816e+00   5.8904e+00   4.0000e-01
   20  00:00:01   1.0885e+03     15          7.6968e+02   4.7147e-03   2.0135e-02   3.2417e+00   5.7530e+00   4.0000e-01
   30  00:00:01   1.0935e+03     16          7.7321e+02   3.0473e-03   2.5432e-02   3.0184e+00   5.6114e+00   4.0000e-01
   40  00:00:01   1.0985e+03     17          7.7674e+02   2.1728e-03   3.0164e-02   2.8804e+00   5.4604e+00   4.0000e-01
   50  00:00:01   1.1035e+03     18          7.8027e+02   1.5503e-03   3.4432e-02   2.7948e+00   5.2955e+00   4.0000e-01
   60  00:00:01   1.1085e+03     19          7.8380e+02   1.0217e-03   3.8295e-02   2.8117e+00   5.0952e+00   4.0000e-01
   70  00:00:01   1.1127e+03     20          7.8676e+02   6.5205e-04   4.1125e-02   3.2504e+00   4.8079e+00   4.0000e-01
   80  00:00:01   1.1132e+03     21          7.8716e+02   6.8670e-04   4.1286e-02   3.5610e+00   4.6607e+00   4.0000e-01
   81  00:00:02   1.1132e+03     22  FP      7.8716e+02   6.8701e-04   4.1286e-02   3.5618e+00   4.6603e+00   4.0000e-01
   81  00:00:02   1.1132e+03     23  SN      7.8716e+02   6.8701e-04   4.1286e-02   3.5618e+00   4.6603e+00   4.0000e-01
   90  00:00:02   1.1132e+03     24          7.8712e+02   7.2344e-04   4.1153e-02   3.6336e+00   4.6139e+00   4.0000e-01
  100  00:00:02   1.1125e+03     25          7.8665e+02   8.5320e-04   4.0378e-02   3.7516e+00   4.4906e+00   4.0000e-01
  110  00:00:02   1.1077e+03     26          7.8326e+02   1.2399e-03   3.5526e-02   3.7439e+00   4.1680e+00   4.0000e-01
  120  00:00:02   1.1027e+03     27          7.7973e+02   1.4291e-03   2.9864e-02   3.6226e+00   3.9414e+00   4.0000e-01
  130  00:00:02   1.0977e+03     28          7.7621e+02   1.4487e-03   2.2192e-02   3.4700e+00   3.7024e+00   4.0000e-01
  140  00:00:02   1.0957e+03     29          7.7478e+02   1.2513e-03   1.4864e-02   3.3571e+00   3.5056e+00   4.0000e-01
  141  00:00:02   1.0957e+03     30  SN      7.7478e+02   1.2452e-03   1.4718e-02   3.3551e+00   3.5018e+00   4.0000e-01
  141  00:00:02   1.0957e+03     31  FP      7.7478e+02   1.2452e-03   1.4718e-02   3.3551e+00   3.5018e+00   4.0000e-01
  150  00:00:02   1.0959e+03     32          7.7488e+02   1.1671e-03   1.2990e-02   3.3332e+00   3.4579e+00   4.0000e-01
  160  00:00:03   1.0999e+03     33          7.7773e+02   8.2087e-04   7.4029e-03   3.2659e+00   3.3198e+00   4.0000e-01
  170  00:00:03   1.1049e+03     34          7.8126e+02   6.4120e-04   5.2784e-03   3.2377e+00   3.2683e+00   4.0000e-01
  180  00:00:03   1.1099e+03     35          7.8480e+02   5.3037e-04   4.1432e-03   3.2209e+00   3.2409e+00   4.0000e-01
  190  00:00:03   1.1149e+03     36          7.8833e+02   4.5320e-04   3.4188e-03   3.2093e+00   3.2235e+00   4.0000e-01
  200  00:00:03   1.1199e+03     37          7.9187e+02   3.9597e-04   2.9128e-03   3.2007e+00   3.2114e+00   4.0000e-01
  210  00:00:03   1.1249e+03     38          7.9540e+02   3.5170e-04   2.5383e-03   3.1941e+00   3.2024e+00   4.0000e-01
  220  00:00:03   1.1299e+03     39          7.9894e+02   3.1639e-04   2.2495e-03   3.1888e+00   3.1955e+00   4.0000e-01
  230  00:00:03   1.1349e+03     40          8.0248e+02   2.8755e-04   2.0200e-03   3.1845e+00   3.1900e+00   4.0000e-01
  240  00:00:03   1.1399e+03     41          8.0601e+02   2.6354e-04   1.8331e-03   3.1809e+00   3.1855e+00   4.0000e-01
  250  00:00:03   1.1449e+03     42          8.0955e+02   2.4324e-04   1.6779e-03   3.1779e+00   3.1818e+00   4.0000e-01
  260  00:00:04   1.1499e+03     43          8.1308e+02   2.2585e-04   1.5469e-03   3.1753e+00   3.1786e+00   4.0000e-01
  270  00:00:04   1.1549e+03     44          8.1662e+02   2.1079e-04   1.4350e-03   3.1730e+00   3.1760e+00   4.0000e-01
  280  00:00:04   1.1599e+03     45          8.2015e+02   1.9761e-04   1.3381e-03   3.1711e+00   3.1736e+00   4.0000e-01
  290  00:00:04   1.1649e+03     46          8.2369e+02   1.8598e-04   1.2535e-03   3.1693e+00   3.1716e+00   4.0000e-01
  300  00:00:04   1.1699e+03     47  EP      8.2722e+02   1.7564e-04   1.1790e-03   3.1678e+00   3.1698e+00   4.0000e-01
Calculate FRC in physical domain at epsilon 2.000000e-01
the forcing frequency 7.2541e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2546e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2582e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2617e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2652e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2688e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2723e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2758e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2794e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2829e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2864e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2900e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2935e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2970e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3006e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3041e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3077e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3112e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3147e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3183e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3218e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3253e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3289e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3324e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3359e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3395e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3430e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3465e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3501e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3536e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3572e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3607e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3642e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3678e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3713e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3748e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3784e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3819e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3854e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3890e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3925e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3960e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3996e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4031e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4066e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4102e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4137e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4173e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4208e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4243e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4279e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4314e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4349e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4385e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4420e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4455e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4491e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4526e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4561e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4597e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4632e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4668e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4703e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4738e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4774e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4809e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4844e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4880e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4915e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4950e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.4986e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5021e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5056e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5092e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5127e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5162e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5198e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5233e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5269e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5304e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5339e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5375e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5410e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5445e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5481e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5516e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5551e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5587e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5622e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5657e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5693e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5728e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5763e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5799e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5834e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5869e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5905e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5940e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.5975e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6010e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6045e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6080e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6115e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6150e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6185e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6219e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6253e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6286e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6318e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6341e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6353e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6359e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6359e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6365e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6377e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6399e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6430e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6461e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6492e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6524e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6557e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6591e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6625e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6660e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6694e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6729e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6764e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6799e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6834e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6869e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6904e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6939e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.6974e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7009e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7044e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7079e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7114e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7148e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7183e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7217e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7251e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7282e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7288e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7293e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7295e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7297e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7298e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7298e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7299e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7299e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7299e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7299e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7299e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7299e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7299e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7299e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7299e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7299e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7298e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7298e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7297e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7296e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7295e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7294e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7293e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7291e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7288e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7284e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7279e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7272e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7260e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7241e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7210e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7196e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7191e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7188e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7186e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7185e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7185e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7185e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7185e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7184e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7184e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7184e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7184e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7184e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7185e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7185e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7185e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7185e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7186e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7186e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7187e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7187e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7188e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7189e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7191e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7193e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7195e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7198e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7203e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7210e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7220e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7236e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7264e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7299e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7334e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7369e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7405e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7440e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7475e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7511e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7546e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7581e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7617e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7652e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7687e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7723e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7758e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7793e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7829e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7864e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7899e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7935e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.7970e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8005e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8041e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8076e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8112e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8147e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8182e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8218e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8253e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8288e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8324e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8359e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8394e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8430e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8465e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8500e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8536e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8571e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8606e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8642e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8677e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8713e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8748e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8783e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8819e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8854e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8889e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8925e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8960e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.8995e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9031e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9066e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9101e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9137e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9172e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9208e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9243e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9278e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9314e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9349e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9384e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9420e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9455e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9490e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9526e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9561e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9596e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9632e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9667e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9702e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9738e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9773e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9809e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9844e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9879e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9915e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9950e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.9985e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0021e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0056e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0091e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0127e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0162e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0197e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0233e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0268e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0304e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0339e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0374e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0410e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0445e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0480e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0516e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0551e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0586e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0622e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0657e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0692e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0728e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0763e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0798e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0834e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0869e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0905e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0940e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.0975e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1011e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1046e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1081e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1117e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1152e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1187e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1223e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1258e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1293e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1329e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1364e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1400e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1435e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1470e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1506e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1541e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1576e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1612e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1647e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1682e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1718e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1753e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1788e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1824e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1859e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1895e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1930e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.1965e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2001e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2036e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2071e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2107e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2142e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2177e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2213e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2248e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2283e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2319e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2354e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2389e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2425e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2460e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2496e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2531e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2566e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2602e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2637e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2672e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2708e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2743e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2778e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2814e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2849e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2884e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2920e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2955e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.2991e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3026e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3061e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3097e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3132e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3167e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3203e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3238e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3273e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3309e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3344e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3379e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3415e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3450e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3485e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3521e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3556e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3592e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3627e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3662e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3698e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3733e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3768e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3804e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3839e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3874e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3910e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3945e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3980e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 8.3995e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
Calculate FRC in physical domain at epsilon 4.000000e-01
the forcing frequency 7.2541e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2574e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2609e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2644e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2680e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2715e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2751e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2786e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2821e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2857e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2892e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2927e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2963e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.2998e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3033e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3069e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3104e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3139e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3175e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3210e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3245e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3281e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3316e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3352e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3387e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3422e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3458e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3493e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3528e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3564e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3599e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcing frequency 7.3634e+02 is nearly resonant with the eigenvalue -1.6662e+00 + i7.6359e+02
the forcin...
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
```

![figure_32.png
](README_images/figure_32.png
)

# Shallow-curved shell structure with geometric nonlinearities

This a preview of the live-script `vonkarmanshell.mlx`. See [1] for the details of this model, and [2] the description of this example.

[1] Jain, S., \& Tiso, P. (2018). Simulation-free hyper-reduction for geometrically nonlinear structural dynamics: a quadratic manifold lifting approach. *Journal of Computational and Nonlinear Dynamics*, *13*(7), 071003. [https://doi.org/10.1115/1.4040021](https://doi.org/10.1115/1.4040021)

[2] Cenedese, M., Marconi, J., Haller, G., \& Jain, S. (2023). Data-assisted non-intrusive model reduction for forced nonlinear finite elements models. Preprint: [arXiv: 2311.17865](https://arxiv.org/abs/2311.17865) 

The finite element code taken from the following package:

Jain, S., Marconi, J., Tiso P. (2020). YetAnotherFEcode (Version v1.1). Zenodo. [http://doi.org/10.5281/zenodo.4011282](http://doi.org/10.5281/zenodo.4011282)

![images/image_0.png
](images/image_0.png
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

![images/figure_0.png
](images/figure_0.png
)

```text:Output
Assembling M,C,K matrices
Applying boundary conditions
Solving undamped eigenvalue problem
```

![images/figure_1.png
](images/figure_1.png
)

![images/figure_2.png
](images/figure_2.png
)

![images/figure_3.png
](images/figure_3.png
)

![images/figure_4.png
](images/figure_4.png
)

```text:Output
Using Rayleigh damping
a = 2x1    
    0.0000
    0.4022

Assembling external force vector
Getting nonlinearity coefficients
Loaded tensors from storage
Total time spent on model assembly = 00:00:13
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

```matlab:Code
masterModes = [1 2]; % Modal displacements and modal velocities
Ve = V(:,masterModes); % Mode shape
We = W(masterModes,:); % Projection to mode shape
Ae = full(We*A*Ve) % Reduced, linearized dynamics
```

```text:Output
Ae = 2x2    
1.0e+04 *

         0    0.0001
   -2.1743   -0.0001

```

```matlab:Code
SSMDim = length(masterModes);

displacementVector = zeros(1,n); displacementVector(outdof) = 1;
loadVector = f_0;  %  could also be set as modal ones
```

# Load static analysis and plot

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

![images/figure_5.png
](images/figure_5.png
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

Define observables and timescales

```matlab:Code
observable = @(x) x; % Observe the full phase space
slowTimeScale = 2*pi/abs(lambda(1));
fastTimeScale = 2*pi/abs(lambda(round(SSMDim/2)));
```

# **Time Integration using Generalized alpha**

The computation of integration time is estimated from the linear decay that gets from the nonlinear amplitude to linear regime.

```matlab:Code
numberPeriodsSlow = floor(log(desiredAmplitudeDecay)/...
    (2*pi*(-real(lambda(1))/abs(lambda(1)))))
```

```text:Output
numberPeriodsSlow = 201
```

```matlab:Code
endTime = numberPeriodsSlow*slowTimeScale;
```

Set the sampling time to capture approximately 100 points per period on the faster time scale

```matlab:Code
numberPeriodsFast = floor(endTime/fastTimeScale);
numberPointsPerPeriod = 100;
nSamp = numberPeriodsFast*numberPointsPerPeriod+1;
dt = endTime/(nSamp-1);
```

Integrate (here we load the data set to avoid excessive waiting time for integration).

```matlab:Code
% xData = integrateTrajectoriesGalphaDirect(Model, DS, endTime, IC, nSamp, observable);
% % xData = integrateTrajectoriesNewmarkDirect(Model, DS, endTime, IC, nSamp, observable);
% DataInfo = struct('nElements', Model.Mesh.nElements, 'loadShape', loadShape);
% save('dataVKDecay2DGalpha.mat', 'DataInfo', 'xData', 'dt', 'endTime', 'nSamp')
load("../../data/vonkarmanshell/dataVKDecay2DGalphaStatic.mat")
```

# Visualize data

```matlab:Code
yData = coordinatesEmbedding(xData, SSMDim, 'ForceEmbedding', 1);
```

```text:Output
The embedding coordinates consist of the measured states.
```

Data filtering: We need to make sure that the data that we use to identify the slow manifold lies close to it. We can do this by plotting a spectrogram of the beam tip displacement. In general, there may be many vibratory modes present at first, but the faster ones quickly die out.

```matlab:Code
showSpectrogram(yData(indTrain,:), outdof);
ylim([0,abs(lambda(1))/2/pi*5])
```

![images/figure_6.png
](images/figure_6.png
)

We plot the tip displacement over time for closer inspection. 

```matlab:Code
customFigure();
plot(xData{1,1}, xData{1,2}(outdof,:), xData{2,1}, xData{2,2}(outdof,:), ':');
xlabel('$t \, [$s$]$','Interpreter','latex'); ylabel('$u \, [$m$]$','Interpreter','latex'); 
legend({'Trajectory 1', 'Trajectory 2'})
title('Generated data')
```

![images/figure_7.png
](images/figure_7.png
)

# Truncate transient data from trajectories

Next, we use the information from the spectrogram and plot to remove the first part of the training data. After the first few oscillations have passed, there is one dominant mode left in the frequency spectrum. In this case, faster modes die out very quickly, so we can use almost all of the data. We must however remove the first transient to fulfill the assumption that trajectories lie close to the SSM. We keep only the time interval |sliceInt|.

```matlab:Code
sliceInt = [5*slowTimeScale, endTime];
yDataTrunc = sliceTrajectories(yData, sliceInt);
```

# Datadriven manifold fitting

The measured trajectories are assumed to lie close to a two-dimensional manifold that is tangent at the origin to the eigenspace corresponding to the slowest pair of eigenvalues. We now fit a polynomial of order $M$ to the data points to approximate the manifold. Here we use a graph style parametrization, meaning that the manifold is parametrized using coordinates $\eta =W_e y$ of the eigenspace represenation $V_e$. This excludes the possibility of folds on the manifold (however, if necessary we may circumvent this problem by increasing the dimensionality of the observed space). 

We seek the $2n\times m_M$ polynomial coefficient matrix $H$ of the manifold on the form

$y=V_e V_e^{\top } y+H\phi_{m,2:M} (V_e^{\top } y)$,

where the function $\phi_{m,2:M} (q)$ computes a vector of all $m_M$ monomials from orders 2 up to $M$ of an $m$-vector $q$. The coefficients $H$are obtained via regression.

```matlab:Code
SSMOrder = 7;

% Get projection or modal coordinates 
SSMDim = size(Ve,2);
etaDataTrunc = yDataTrunc;
nTraj = size(yDataTrunc,1);
for iTraj = 1:nTraj
    etaDataTrunc{iTraj,2} = We*yDataTrunc{iTraj,2};    
end

% Plot reduced coordinates
plotReducedCoordinates(etaDataTrunc);
legend({'Test set trajectory', 'Training set trajectory'})
```

![images/figure_8.png
](images/figure_8.png
)

```matlab:Code
if SSMDim>2
   view(3) 
end

% Compute nonlinear part of the parametrization
IMInfo = IMGeometry(yDataTrunc(indTrain,:), SSMDim,SSMOrder,...
         'reducedCoordinates',etaDataTrunc(indTrain,:),'Ve',Ve,'outdof',outdof); 
IMInfo.chart.map = @(x) We*x;                          

% Parametrization error on test trajectory
normedTrajDist = computeTrajectoryErrors(liftTrajectories(IMInfo,...
    etaDataTrunc), yDataTrunc);
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
plotSSMandTrajectories(IMInfo, idxPlot, yDataTrunc(indTest,:), etaDataTrunc(indTest,:))
view(-100,20); legend('off')
```

![images/figure_9.png
](images/figure_9.png
)

```matlab:Code
gPlot = @(x) W(5,:)*x; % 3D Plot: eta_1, eta_2 and gPlot(x) 
plotSSMandTrajectories(IMInfo, gPlot, yDataTrunc(indTest,:), etaDataTrunc(indTest,:))
view(-100,20); legend('off')
```

![images/figure_10.png
](images/figure_10.png
)

```matlab:Code
gPlot = @(x) W([3],:)*x; % 3D Plot with the three values of gPlot(x) 
plotSSMandTrajectories(IMInfo, gPlot, yDataTrunc(indTest,:), etaDataTrunc(indTest,:))
view(-100,20); legend('off')
xlabel('$u_1$','interpreter','latex')
ylabel('$\dot{u}_1$','interpreter','latex')
zlabel('$u_2$','interpreter','latex')
```

![images/figure_11.png
](images/figure_11.png
)

# Reduced dynamics on the manifold

We compute a model for the reduced dynamics with the truncated training data projected onto the manifold. The function `IMDynamicsFlow` fits a polynomial map = 

$$
\dot{\eta} =W_r \phi (\eta )
$$

where $\phi (\eta )$ again computes a vector of all monomials of $\eta$, and  $W_r$ is a matrix of polynomial coefficients. 

We are also specifying that we want the reduced dynamics on a normal form, and seek the Taylor expansion of a map $N$ fulfilling

$$
\dot{z} =N(z)=Dz+W_n \phi (z)
$$

with $D$ a diagonal matrix and $W_n$ containing coefficients for the near-resonant nonlinear terms, after a near-identity change of coordinates

$$
z=T^{-1} (\eta )=\eta +W_t \phi (\eta )
$$

with $W_t$ containing the coefficients for the nonlinear terms of the coordinate change. Here we fix the linear part.

```matlab:Code
ROMOrder = 7;

RDInfo = IMDynamicsMech(etaDataTrunc(indTrain,:), ...
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
    13          16         0.944215              1          0.551  
    14          17         0.940702              1          0.859  
    15          18         0.938716              1           2.62  
    16          19         0.933349              1          0.564  
    17          20         0.930356              1          0.643  
    18          21         0.929047              1          0.413  
    19          22         0.928417              1          0.224  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    20          23         0.927589              1          0.333  
    21          24         0.926516              1          0.566  
    22          25         0.925183              1           0.54  
    23          26         0.924128              1          0.274  
    24          27         0.923647              1         0.0886  
    25          28         0.923482              1         0.0908  
    26          29         0.923359              1         0.0915  
    27          30         0.923144              1          0.134  
    28          31         0.922795              1          0.178  
    29          32         0.922322              1          0.182  
    30          33         0.921845              1          0.206  
    31          34         0.921435              1          0.225  
    32          35         0.921035              1          0.206  
    33          36         0.920547              1           0.21  
    34          37         0.919969              1          0.216  
    35          38         0.919428              1          0.228  
    36          39         0.919057              1          0.151  
    37          40         0.918835              1          0.106  
    38          41         0.918665              1          0.102  
    39          42         0.918457              1          0.142  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    40          43         0.918128              1          0.221  
    41          44         0.917592              1          0.272  
    42          45         0.916876              1          0.249  
    43          46         0.916256              1          0.221  
    44          47         0.915962              1           0.12  
    45          48         0.915867              1         0.0538  
    46          49         0.915812              1         0.0534  
    47          50         0.915723              1          0.105  
    48          51         0.915547              1          0.184  
    49          52         0.915164              1          0.276  
    50          53         0.914397              1          0.354  
    51          54         0.913193              1          0.413  
    52          55         0.912043              1          0.305  
    53          56           0.9115              1         0.0998  
    54          57         0.911328              1         0.0913  
    55          58         0.911183              1          0.102  
    56          59         0.910817              1          0.222  
    57          60         0.909964              1           0.39  
    58          61          0.90802              1          0.603  
    59          62         0.904513              1          0.752  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    60          63         0.900407              1          0.634  
    61          64         0.898125              1          0.274  
    62          65         0.897627              1         0.0811  
    63          66         0.897549              1         0.0826  
    64          67         0.897468              1         0.0817  
    65          68         0.897241              1          0.108  
    66          69         0.896767              1           0.17  
    67          70         0.895913              1          0.211  
    68          71         0.894975              1           0.17  
    69          72         0.894472              1         0.0648  
    70          73         0.894346              1         0.0625  
    71          74         0.894298              1         0.0595  
    72          75         0.894202              1          0.057  
    73          76         0.893977              1         0.0884  
    74          77         0.893448              1           0.17  
    75          78         0.892279              1          0.353  
    76          79         0.889843              1          0.635  
    77          80         0.885188              1          0.994  
    78          81         0.877277              1           1.31  
    79          82         0.866194              1           1.36  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    80          83         0.855044              1          0.948  
    81          84          0.84831              1          0.326  
    82          85         0.846116              1          0.133  
    83          86         0.845646              1           0.11  
    84          87         0.845469              1          0.116  
    85          88         0.845214              1          0.121  
    86          89         0.844689              1          0.128  
    87          90         0.843619              1          0.174  
    88          91         0.841841              1            0.3  
    89          92         0.839679              1          0.502  
    90          93         0.837646              1          0.632  
    91          94          0.83545              1          0.675  
    92          95         0.831761              1          0.652  
    93          96         0.824525              1          0.538  
    94          97         0.812762              1          0.686  
    95          98          0.80108              1          0.536  
    96          99         0.795932              1          0.192  
    97         100         0.795093              1         0.0873  
    98         101         0.794995              1         0.0863  
    99         102         0.794906              1         0.0848  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   100         103          0.79469              1         0.0816  
   101         104         0.794318              1          0.164  
   102         105         0.793701              1          0.278  
   103         106         0.792888              1          0.342  
   104         107          0.79201              1          0.292  
   105         108         0.791245              1          0.202  
   106         109         0.790622              1           0.22  
   107         110         0.790017              1          0.224  
   108         111         0.789224              1          0.337  
   109         112         0.788048              1          0.385  
   110         113         0.786511              1          0.305  
   111         114         0.785034              1          0.318  
   112         115         0.783965              1          0.346  
   113         116         0.783117              1          0.315  
   114         117         0.782136              1          0.252  
   115         118         0.781019              1          0.208  
   116         119         0.780212              1          0.189  
   117         120           0.7799              1          0.112  
   118         121           0.7798              1         0.0578  
   119         122         0.779715              1         0.0579  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   120         123         0.779557              1          0.152  
   121         124         0.779258              1          0.281  
   122         125         0.778691              1          0.407  
   123         126         0.777825              1          0.444  
   124         127         0.776936              1          0.303  
   125         128          0.77641              1          0.176  
   126         129         0.776155              1          0.176  
   127         130         0.775883              1          0.178  
   128         131         0.775273              1          0.304  
   129         132         0.773825              1          0.462  
   130         133         0.770366              1          0.648  
   131         134         0.762789              1          0.844  
   132         135          0.74882              1           1.26  
   133         136         0.730663              1           1.41  
   134         137         0.717222              1           1.04  
   135         138         0.712049              1          0.477  
   136         139         0.710661              1          0.207  
   137         140         0.710152              1           0.11  
   138         141          0.70988              1          0.108  
   139         142         0.709663              1          0.108  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   140         143         0.709254              1          0.147  
   141         144         0.708285              1          0.246  
   142         145         0.705965              1          0.374  
   143         146         0.701027              1          0.505  
   144         147         0.692967              1          0.531  
   145         148         0.685014              1          0.342  
   146         149         0.681128              1           0.36  
   147         150         0.679811              1          0.356  
   148         151         0.678593              1           0.34  
   149         152         0.675498              1          0.307  
   150         153         0.668172              1          0.434  
   151         154         0.651307              1          0.636  
   152         155         0.620269              1          0.771  
   153         156          0.58292              1          0.633  
   154         157         0.561564              1          0.292  
   155         158         0.556851              1          0.177  
   156         159         0.556349              1          0.165  
   157         160         0.556099              1          0.161  
   158         161         0.555219              1          0.154  
   159         162          0.55331              1          0.242  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   160         163         0.548539              1          0.498  
   161         164         0.538721              1          0.834  
   162         165         0.522841              1           1.08  
   163         166         0.507747              1          0.897  
   164         167         0.500974              1          0.431  
   165         168         0.499369              1          0.418  
   166         169         0.498606              1          0.399  
   167         170         0.496948              1          0.382  
   168         171         0.493557              1          0.758  
   169         172         0.486848              1           1.19  
   170         173         0.477104              1           1.36  
   171         174         0.468702              1           0.96  
   172         175         0.465451              1          0.332  
   173         176         0.464937              1         0.0577  
   174         177          0.46488              1         0.0577  
   175         178         0.464848              1         0.0577  
   176         179         0.464786              1         0.0578  
   177         180         0.464652              1         0.0578  
   178         181         0.464317              1         0.0578  
   179         182         0.463534              1         0.0642  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   180         183         0.461902              1         0.0666  
   181         184         0.459401              1          0.105  
   182         185         0.457215              1          0.163  
   183         186         0.456359              1          0.176  
   184         187         0.456161              1          0.169  
   185         188         0.456016              1          0.161  
   186         189         0.455609              1          0.145  
   187         190         0.454621              1          0.121  
   188         191         0.452034              1         0.0808  
   189         192         0.445733              1         0.0919  
   190         193         0.431607              1          0.135  
   191         194         0.406484              1          0.161  
   192         195         0.378406              1          0.173  
   193         196         0.363954              1         0.0894  
   194         197         0.361255              1         0.0908  
   195         198         0.361068              1         0.0891  
   196         199         0.361021              1         0.0883  
   197         200         0.360836              1         0.0863  
   198         201         0.360429              1         0.0833  
   199         202         0.359314              1         0.0775  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   200         203         0.356618              1         0.0676  
   201         204         0.350418              1         0.0755  
   202         205         0.338892              1         0.0941  
   203         206         0.324721              1         0.0789  
   204         207         0.316351              1         0.0825  
   205         208         0.314524              1         0.0774  
   206         209         0.314354              1         0.0732  
   207         210         0.314295              1         0.0717  
   208         211         0.314087              1         0.0677  
   209         212         0.313623              1         0.0617  
   210         213         0.312423              1         0.0507  
   211         214         0.309823              1         0.0397  
   212         215         0.305245              1         0.0396  
   213         216         0.300323              1         0.0541  
   214         217         0.297911              1           0.06  
   215         218         0.297486              1         0.0513  
   216         219         0.297448              1         0.0472  
   217         220         0.297425              1         0.0452  
   218         221         0.297343              1         0.0402  
   219         222          0.29716              1         0.0373  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   220         223         0.296694              1         0.0369  
   221         224         0.295718              1          0.036  
   222         225         0.294127              1         0.0343  
   223         226         0.292637              1         0.0337  
   224         227         0.292033              1          0.034  
   225         228         0.291943              1         0.0321  
   226         229          0.29193              1         0.0322  
   227         230         0.291908              1         0.0323  
   228         231         0.291849              1         0.0324  
   229         232           0.2917              1         0.0324  
   230         233         0.291326              1         0.0321  
   231         234         0.290468              1         0.0307  
   232         235          0.28882              1         0.0449  
   233         236         0.286689              1         0.0736  
   234         237         0.285306              1         0.0795  
   235         238         0.284956              1         0.0687  
   236         239         0.284913              1         0.0623  
   237         240          0.28489              1         0.0593  
   238         241         0.284814              1         0.0525  
   239         242         0.284644              1         0.0422  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   240         243         0.284222              1         0.0375  
   241         244         0.283382              1         0.0333  
   242         245         0.282156              1         0.0237  
   243         246          0.28121              1         0.0259  
   244         247          0.28092              1         0.0135  
   245         248          0.28089              1         0.0138  
   246         249         0.280888              1         0.0137  
   247         250         0.280887              1         0.0137  
   248         251         0.280881              1         0.0135  
   249         252         0.280867              1         0.0132  
   250         253         0.280831              1         0.0123  
   251         254         0.280751              1         0.0106  
   252         255         0.280598              1         0.0148  
   253         256          0.28041              1          0.014  
   254         257         0.280296              1        0.00715  
   255         258          0.28027              1        0.00685  
   256         259         0.280268              1        0.00636  
   257         260         0.280268              1        0.00622  
   258         261         0.280267              1        0.00592  
   259         262         0.280265              1        0.00545  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   260         263         0.280259              1        0.00455  
   261         264         0.280247              1        0.00356  
   262         265         0.280223              1        0.00333  
   263         266         0.280193              1        0.00208  
   264         267         0.280175              1        0.00192  
   265         268         0.280171              1        0.00164  
   266         269         0.280171              1        0.00165  
   267         270         0.280171              1        0.00165  
   268         271         0.280171              1        0.00165  
   269         272          0.28017              1        0.00165  
   270         273          0.28017              1        0.00166  
   271         274         0.280168              1        0.00166  
   272         275         0.280165              1        0.00166  
   273         276         0.280156              1        0.00222  
   274         277         0.280135              1        0.00365  
   275         278         0.280096              1        0.00483  
   276         279         0.280049              1        0.00539  
   277         280         0.280021              1        0.00722  
   278         281         0.280015              1        0.00737  
   279         282         0.280014              1        0.00715  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   280         283         0.280014              1        0.00704  
   281         284         0.280013              1        0.00678  
   282         285         0.280011              1         0.0064  
   283         286         0.280005              1        0.00576  
   284         287         0.279989              1        0.00471  
   285         288          0.27995              1          0.003  
   286         289         0.279854              1         0.0045  
   287         290         0.279635              1         0.0066  
   288         291         0.279235              1        0.00803  
   289         292         0.278761              1        0.00789  
   290         293         0.278497              1        0.00412  
   291         294         0.278444              1        0.00179  
   292         295         0.278441              1        0.00174  
   293         296         0.278441              1        0.00173  
   294         297         0.278441              1        0.00173  
   295         298         0.278441              1        0.00172  
   296         299         0.278441              1        0.00171  
   297         300          0.27844              1        0.00169  
   298         301         0.278438              1        0.00164  
   299         302         0.278434              1        0.00156  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   300         303         0.278424              1        0.00139  
   301         304         0.278401              1        0.00105  
   302         305         0.278355              1       0.000878  
   303         306         0.278295              1        0.00103  
   304         307         0.278254              1        0.00138  
   305         308         0.278244              1        0.00141  
   306         309         0.278243              1        0.00136  
   307         310         0.278243              1        0.00135  
   308         311         0.278243              1        0.00134  
   309         312         0.278243              1        0.00132  
   310         313         0.278243              1        0.00129  
   311         314         0.278242              1        0.00121  
   312         315         0.278241              1        0.00106  
   313         316         0.278239              1       0.000721  
   314         317         0.278235              1       0.000573  
   315         318         0.278231              1       0.000624  
   316         319         0.278228              1        0.00101  
   317         320         0.278227              1       0.000955  
   318         321         0.278227              1       0.000867  
   319         322         0.278227              1       0.000845  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   320         323         0.278227              1       0.000805  
   321         324         0.278227              1       0.000745  
   322         325         0.278227              1       0.000645  
   323         326         0.278226              1       0.000559  
   324         327         0.278225              1       0.000559  
   325         328         0.278223              1       0.000558  
   326         329         0.278217              1       0.000855  
   327         330         0.278203              1         0.0018  
   328         331         0.278175              1        0.00286  
   329         332         0.278134              1        0.00326  
   330         333         0.278101              1        0.00317  
   331         334         0.278091              1        0.00351  
   332         335          0.27809              1        0.00343  
   333         336          0.27809              1        0.00338  
   334         337          0.27809              1        0.00335  
   335         338         0.278089              1        0.00328  
   336         339         0.278089              1        0.00317  
   337         340         0.278088              1        0.00299  
   338         341         0.278084              1         0.0027  
   339         342         0.278075              1         0.0022  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   340         343         0.278052              1        0.00138  
   341         344            0.278              1         0.0014  
   342         345         0.277906              1        0.00183  
   343         346         0.277798              1        0.00195  
   344         347          0.27774              1        0.00122  
   345         348         0.277729              1        0.00072  
   346         349         0.277728              1       0.000726  
   347         350         0.277728              1       0.000724  
   348         351         0.277728              1       0.000723  
   349         352         0.277728              1        0.00072  
   350         353         0.277728              1       0.000716  
   351         354         0.277728              1       0.000705  
   352         355         0.277728              1       0.000682  
   353         356         0.277727              1        0.00063  
   354         357         0.277726              1       0.000513  
   355         358         0.277723              1       0.000352  
   356         359         0.277718              1       0.000346  
   357         360         0.277714              1       0.000395  
   358         361         0.277712              1       0.000449  
   359         362         0.277712              1       0.000399  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   360         363         0.277712              1        0.00038  
   361         364         0.277712              1       0.000376  
   362         365         0.277712              1       0.000369  
   363         366         0.277712              1       0.000369  
   364         367         0.277712              1        0.00037  
   365         368         0.277712              1        0.00037  
   366         369         0.277712              1       0.000369  
   367         370         0.277711              1       0.000367  
   368         371         0.277709              1       0.000361  
   369         372         0.277705              1       0.000343  
   370         373         0.277695              1        0.00069  
   371         374         0.277677              1        0.00103  
   372         375         0.277653              1        0.00102  
   373         376         0.277638              1       0.000538  
   374         377         0.277635              1       0.000113  
   375         378         0.277634              1       0.000103  
   376         379         0.277634              1       0.000103  
   377         380         0.277634              1       0.000103  
   378         381         0.277634              1       0.000103  
   379         382         0.277634              1       0.000103  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   380         383         0.277634              1       0.000103  
   381         384         0.277634              1       0.000103  
   382         385         0.277634              1       0.000103  
   383         386         0.277634              1       0.000103  
   384         387         0.277634              1       0.000103  
   385         388         0.277634              1       0.000103  
   386         389         0.277634              1       0.000103  
   387         390         0.277633              1       0.000103  
   388         391         0.277633              1       0.000108  
   389         392         0.277633              1       0.000126  
   390         393         0.277633              1       0.000117  
   391         394         0.277633              1       0.000113  
   392         395         0.277633              1       0.000113  
   393         396         0.277633              1       0.000111  
   394         397         0.277633              1       0.000109  
   395         398         0.277633              1       0.000105  
   396         399         0.277633              1       0.000103  
   397         400         0.277633              1       0.000103  
   398         401         0.277633              1       0.000103  
   399         402         0.277633              1       0.000103  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   400         403         0.277632              1       0.000103  
   401         404         0.277632              1       0.000103  
   402         405          0.27763              1       0.000168  
   403         406         0.277626              1       0.000325  
   404         407         0.277618              1       0.000516  
   405         408         0.277603              1       0.000634  
   406         409         0.277589              1        0.00051  
   407         410         0.277584              1       0.000279  
   408         411         0.277583              1       0.000272  
   409         412         0.277583              1       0.000266  
   410         413         0.277583              1       0.000265  
   411         414         0.277583              1       0.000265  
   412         415         0.277583              1       0.000263  
   413         416         0.277583              1       0.000262  
   414         417         0.277583              1       0.000259  
   415         418         0.277583              1       0.000255  
   416         419         0.277583              1       0.000248  
   417         420         0.277582              1       0.000236  
   418         421         0.277581              1       0.000218  
   419         422         0.277579              1       0.000191  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   420         423         0.277575              1       0.000158  
   421         424         0.277569              1       0.000369  
   422         425         0.277564              1       0.000522  
   423         426         0.277562              1       0.000526  
   424         427         0.277562              1        0.00049  
   425         428         0.277562              1       0.000479  
   426         429         0.277562              1       0.000477  
   427         430         0.277562              1       0.000471  
   428         431         0.277562              1       0.000463  
   429         432         0.277562              1       0.000448  
   430         433         0.277562              1       0.000425  
   431         434         0.277562              1       0.000387  
   432         435         0.277561              1       0.000325  
   433         436         0.277559              1       0.000225  
   434         437         0.277554              1       0.000208  
   435         438         0.277541              1       0.000182  
   436         439         0.277514              1       0.000492  
   437         440         0.277472              1       0.000707  
   438         441         0.277434              1       0.000574  
   439         442         0.277419              1       0.000316  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   440         443         0.277418              1        0.00032  
   441         444         0.277417              1       0.000318  
   442         445         0.277417              1       0.000317  
   443         446         0.277417              1       0.000317  
   444         447         0.277417              1       0.000317  
   445         448         0.277417              1       0.000316  
   446         449         0.277417              1       0.000314  
   447         450         0.277417              1       0.000312  
   448         451         0.277417              1       0.000308  
   449         452         0.277416              1       0.000301  
   450         453         0.277414              1        0.00029  
   451         454         0.277409              1       0.000446  
   452         455         0.277396              1       0.000725  
   453         456         0.277363              1        0.00116  
   454         457         0.277279              1        0.00182  
   455         458         0.277086              1        0.00265  
   456         459          0.27672              1        0.00331  
   457         460         0.276252              1        0.00293  
   458         461         0.275958              1        0.00137  
   459         462         0.275889              1       0.000257  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   460         463         0.275884              1       7.41e-05  
   461         464         0.275884              1       7.36e-05  
   462         465         0.275884              1       7.35e-05  
   463         466         0.275884              1       7.35e-05  
   464         467         0.275884              1       7.35e-05  
   465         468         0.275884              1       7.35e-05  
   466         469         0.275884              1       0.000151  
   467         470         0.275884              1        0.00028  
   468         471         0.275884              1       0.000498  
   469         472         0.275884              1        0.00085  
   470         473         0.275884              1        0.00141  
   471         474         0.275884              1        0.00228  
   472         475         0.275883              1         0.0035  
   473         476         0.275883              1        0.00479  
   474         477         0.275882              1        0.00511  
   475         478         0.275882              1        0.00337  
   476         479         0.275882              1          0.001  
   477         480         0.275881              1       7.89e-05  
   478         481         0.275881              1       1.26e-05  
   479         482         0.275881              1       1.24e-05  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   480         483         0.275881              1       1.24e-05  
   481         484         0.275881              1       1.22e-05  
   482         485         0.275881              1       1.52e-05  
   483         486         0.275881              1       2.99e-05  
   484         487         0.275881              1       5.23e-05  
   485         488         0.275881              1       8.91e-05  
   486         489         0.275881              1       0.000147  
   487         490         0.275881              1       0.000238  
   488         491         0.275881              1       0.000371  
   489         492         0.275881              1       0.000536  
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

![images/figure_12.png
](images/figure_12.png
)

```matlab:Code

% We transform the truncated initial condition of our test trajectory according to 
% the obtained change of coordinates, and integrate our reduced order evolution rule 
% to predict the development of the trajectory. 
[yRec, etaRec, zRec] = advect(IMInfo, RDInfo, yDataTrunc);

% Evaluation of reduced dynamics
% The error NMTE is computed as the average distance of the predicted trajectory 
% to the measured one in the full state space.
normedTrajDist = computeTrajectoryErrors(yRec, yDataTrunc);
NMTE = mean(normedTrajDist(indTest))*100;
disp(['Normalized mean trajectory error = ' num2str(NMTE) '%'])
```

```text:Output
Normalized mean trajectory error = 6.3316%
```

```matlab:Code

% We plot the true test set trajectory in the reduced coordinates and compare it to the prediction. 
plotReducedCoordinates(etaDataTrunc(indTest,:), etaRec(indTest,:))
legend({'Test set (truncated)', 'Prediction'})
```

![images/figure_13.png
](images/figure_13.png
)

```matlab:Code
if size(Ae,1)==2
    % Plot SSM with trajectories in the normal form reduced coordinates
    plotSSMandTrajectories(IMInfo, outdof, yDataTrunc(indTest,:), ...
        zRec(indTest,:), 'NFT', RDInfo.transformation.map)
    view(-100,20); legend('off')
else
    view(3)
end
```

![images/figure_14.png
](images/figure_14.png
)

We also plot the measured and predicted tip displacement. The reduced model seems to do well on previously unseen data, provided that it is close to the 2D manifold.

```matlab:Code
plotTrajectories(yData(indTest,:), yRec(indTest,:), 'm','PlotCoordinate', outdof(1), 'DisplayName', {'Test set', 'Prediction'})
xlabel('$t \, [$s$]$','Interpreter','latex')
ylabel('$q_A \, [$m$]$','Interpreter','latex')
```

![images/figure_15.png
](images/figure_15.png
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
   1   1  1.00e+00  8.40e-04  5.48e-08  2.09e+02    0.0    0.0    0.0
   2   1  1.00e+00  1.12e-07  2.88e-15  2.09e+02    0.0    0.0    0.0

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1          th1          eps
    0  00:00:00   2.0867e+02      1  EP      1.4745e+02   1.1416e-01   5.3150e+00   1.0000e-01
   10  00:00:00   2.0702e+02      2          1.4625e+02   1.8703e-01   6.3141e+00   1.0000e-01
   14  00:00:00   2.0700e+02      3  FP      1.4623e+02   1.8456e-01   6.4671e+00   1.0000e-01
   14  00:00:00   2.0700e+02      4  SN      1.4623e+02   1.8456e-01   6.4671e+00   1.0000e-01
   20  00:00:00   2.0718e+02      5          1.4634e+02   1.6262e-01   6.8494e+00   1.0000e-01
   27  00:00:00   2.0752e+02      6  SN      1.4656e+02   1.0226e-01   7.3309e+00   1.0000e-01
   27  00:00:00   2.0752e+02      7  FP      1.4656e+02   1.0226e-01   7.3309e+00   1.0000e-01
   30  00:00:00   2.0751e+02      8          1.4655e+02   9.1414e-02   7.3955e+00   1.0000e-01
   40  00:00:00   2.0733e+02      9          1.4641e+02   6.6737e-02   7.5309e+00   1.0000e-01
   50  00:00:00   2.0458e+02     10          1.4445e+02   2.0906e-02   7.7574e+00   1.0000e-01
   60  00:00:00   1.9959e+02     11          1.4092e+02   9.5960e-03   7.8109e+00   1.0000e-01
   70  00:00:00   1.9460e+02     12          1.3738e+02   6.2298e-03   7.8267e+00   1.0000e-01
   75  00:00:01   1.9217e+02     13  EP      1.3566e+02   5.3207e-03   7.8310e+00   1.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1          th1          eps
    0  00:00:01   2.0867e+02     14  EP      1.4745e+02   1.1416e-01   5.3150e+00   1.0000e-01
   10  00:00:01   2.1263e+02     15          1.5028e+02   2.2000e-02   4.8181e+00   1.0000e-01
   20  00:00:01   2.1763e+02     16          1.5381e+02   9.8609e-03   4.7607e+00   1.0000e-01
   30  00:00:01   2.2262e+02     17          1.5735e+02   6.3429e-03   4.7442e+00   1.0000e-01
   32  00:00:01   2.2323e+02     18  EP      1.5778e+02   6.0789e-03   4.7430e+00   1.0000e-01
```

![images/figure_16.png
](images/figure_16.png
)

![images/figure_17.png
](images/figure_17.png
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
   80  00:00:00   1.9816e+02     13          1.3990e+02   1.6632e-02   7.8169e+00   2.0000e-01
   90  00:00:00   1.9317e+02     14          1.3637e+02   1.1324e-02   7.8294e+00   2.0000e-01
   93  00:00:00   1.9217e+02     15  EP      1.3566e+02   1.0644e-02   7.8310e+00   2.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1          th1          eps
    0  00:00:00   2.0866e+02     16  EP      1.4745e+02   1.4968e-01   5.1079e+00   2.0000e-01
   10  00:00:01   2.1276e+02     17          1.5037e+02   4.1978e-02   4.8139e+00   2.0000e-01
   20  00:00:01   2.1776e+02     18          1.5390e+02   1.9408e-02   4.7601e+00   2.0000e-01
   30  00:00:01   2.2275e+02     19          1.5744e+02   1.2562e-02   4.7439e+00   2.0000e-01
   31  00:00:01   2.2323e+02     20  EP      1.5778e+02   1.2153e-02   4.7430e+00   2.0000e-01
```

![images/figure_18.png
](images/figure_18.png
)

![images/figure_19.png
](images/figure_19.png
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
   1.0e+02 *

  -0.0029 + 1.4745i
  -0.0029 - 1.4745i
  -0.0063 + 3.1598i
  -0.0063 - 3.1598i
  -0.0094 + 4.1319i
  -0.0094 - 4.1319i
  -0.0148 + 5.4515i
  -0.0148 - 5.4515i
  -0.0173 + 5.9601i
  -0.0173 - 5.9601i

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
   1.0e+02 *

  -0.0063 + 3.1598i
  -0.0063 + 3.1598i
  -0.0063 - 3.1598i
  -0.0063 - 3.1598i
  -0.0094 + 4.1319i
  -0.0094 + 4.1319i
  -0.0094 - 4.1319i
  -0.0094 - 4.1319i
  -0.0148 + 5.4515i
  -0.0148 - 5.4515i
  -0.0173 + 5.9601i
  -0.0173 - 5.9601i

sigma_out = 5
(near) inner resonance detected for the following combination of master eigenvalues
     2     1
     3     2
     1     2
     2     3

These are in resonance with the follwing eigenvalues of the master subspace
   1.0e+02 *

  -0.0029 + 1.4745i
  -0.0029 + 1.4745i
  -0.0029 - 1.4745i
  -0.0029 - 1.4745i

sigma_in = 5
Due to (near) outer resonance, the exisitence of the manifold is questionable and the underlying computation may suffer.
Attempting manifold computation
Manifold computation time at order 2 = 00:00:00
Estimated memory usage at order  2 = 1.30E+01 MB
Manifold computation time at order 3 = 00:00:00
Estimated memory usage at order  3 = 1.41E+01 MB
Manifold computation time at order 4 = 00:00:01
Estimated memory usage at order  4 = 1.52E+01 MB
Manifold computation time at order 5 = 00:00:01
Estimated memory usage at order  5 = 1.70E+01 MB
Manifold computation time at order 6 = 00:00:02
Estimated memory usage at order  6 = 1.91E+01 MB
Manifold computation time at order 7 = 00:00:03
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
    0  00:00:00   2.0860e+02     14  EP      1.4745e+02   3.9582e-02   3.7344e+00   1.0000e-01
   10  00:00:00   2.1260e+02     15          1.5030e+02   7.5474e-03   3.2445e+00   1.0000e-01
   20  00:00:00   2.1760e+02     16          1.5383e+02   3.3968e-03   3.1878e+00   1.0000e-01
   30  00:00:00   2.2260e+02     17          1.5737e+02   2.1876e-03   3.1713e+00   1.0000e-01
   32  00:00:00   2.2318e+02     18  EP      1.5778e+02   2.1012e-03   3.1702e+00   1.0000e-01

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
   90  00:00:00   1.9488e+02     14          1.3766e+02   4.4293e-03   6.2531e+00   2.0000e-01
   96  00:00:00   1.9205e+02     15  EP      1.3566e+02   3.6791e-03   6.2582e+00   2.0000e-01

 STEP      TIME        ||U||  LABEL  TYPE            om         rho1          th1          eps
    0  00:00:00   2.0859e+02     16  EP      1.4745e+02   5.1962e-02   3.5277e+00   2.0000e-01
   10  00:00:00   2.1271e+02     17          1.5037e+02   1.4494e-02   3.2408e+00   2.0000e-01
   20  00:00:01   2.1770e+02     18          1.5391e+02   6.7051e-03   3.1872e+00   2.0000e-01
   30  00:00:01   2.2270e+02     19          1.5744e+02   4.3408e-03   3.1711e+00   2.0000e-01
   31  00:00:01   2.2318e+02     20  EP      1.5778e+02   4.2007e-03   3.1702e+00   2.0000e-01
Calculate FRC in physical domain at epsilon 1.000000e-01
the forcing frequency 1.3566e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3597e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3632e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3667e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3703e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3738e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3773e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3809e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3844e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3879e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3915e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3950e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3986e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4021e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4056e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4092e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4127e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4162e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4198e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4233e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4268e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4304e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4339e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4374e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4410e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4445e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4480e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4516e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4551e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4583e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4604e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4617e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4626e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4632e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4637e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4641e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4644e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4646e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4648e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4649e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4651e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4652e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4653e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4653e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4654e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4655e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4655e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4655e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4656e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4656e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4656e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4656e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4656e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4655e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4653e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4649e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4641e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4634e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4629e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4625e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4623e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4621e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4620e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4619e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4619e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4619e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4619e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4619e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4619e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4620e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4621e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4622e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4625e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4629e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4635e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4648e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4674e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4704e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4727e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4739e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4745e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4745e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4752e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4764e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4788e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4822e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4854e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4889e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4924e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4959e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4994e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5030e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5065e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5100e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5136e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5171e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5207e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5242e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5277e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5313e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5348e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5383e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5419e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5454e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5489e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5525e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5560e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5595e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5631e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5666e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5701e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5737e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5772e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5778e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
Calculate FRC in physical domain at epsilon 2.000000e-01
the forcing frequency 1.3566e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3589e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3624e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3659e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3695e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3730e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3766e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3801e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3836e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3872e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3907e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3942e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.3978e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4013e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4048e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4084e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4119e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4154e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4190e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4225e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4261e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4296e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4331e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4367e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4402e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4437e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4473e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4508e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4539e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4559e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4571e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4579e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4583e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4586e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4589e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4590e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4592e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4593e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4593e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4594e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4594e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4595e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4595e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4596e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4596e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4596e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4596e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4596e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4596e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4596e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4596e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4595e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4595e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4594e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4591e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4587e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4575e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4559e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4526e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4493e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4460e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4437e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4422e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4412e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4405e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4400e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4396e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4393e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4391e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4389e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4388e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4387e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4386e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4386e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4385e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4385e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4385e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4385e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4385e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4385e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4385e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4385e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4385e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4386e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4387e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4388e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4389e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4392e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4399e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4429e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4461e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4495e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4529e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4563e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4597e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4631e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4666e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4700e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4725e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4739e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4745e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4745e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4752e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4766e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4791e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4826e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4861e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4896e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4931e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.4967e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5002e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5037e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5073e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5108e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5143e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5179e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5214e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5249e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5285e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5320e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5355e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5391e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5426e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5461e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5497e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5532e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5568e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5603e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5638e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5674e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5709e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5744e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
the forcing frequency 1.5778e+02 is nearly resonant with the eigenvalue -2.9491e-01 + i1.4745e+02
```

![images/figure_20.png
](images/figure_20.png
)

![images/figure_21.png
](images/figure_21.png
)

![images/figure_22.png
](images/figure_22.png
)

![images/figure_23.png
](images/figure_23.png
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

![images/figure_24.png
](images/figure_24.png
)

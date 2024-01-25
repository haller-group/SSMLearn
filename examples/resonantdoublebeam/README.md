This is a preview of the livescript `resonantdoublebeam.mlx`.

# Finding a 4D SSM from resonant vibration data

We analyze the dynamics along a slow 4D resonant SSM. This example uses velocity measurements from laser scanner vibrometry [1].

```matlab:Code
clearvars
close all
format shortg
clc
```

[1] M. Cenedese, J. Axås, H. Yang, M. Eriten \& G. Haller, Data-driven nonlinear model reduction to spectral  submanifolds in mechanical systems, *arXiv:2110.01929 *(2022). [https://arxiv.org/abs/2110.01929](https://arxiv.org/abs/2110.01929) 

# Example setup

The data consists of 12 decaying velocity signals measured on the tip of the inner beam, excited via hammer impacts. These occur on three positions (4 repetitions per impact) along the inner beam, spaced from the joint to the free end. An initial data inspection reveal two fundamental resonant frequencies, hereby justifying the identification of the slow 4D SSM of the system. Additional information can be found in [1].

![image_0.png
](README_images/image_0.png
)

```matlab:Code
load('data.mat')
% Plot data
customFigure; colororder(jet(size(xData,1)))
for indTraj = 1:size(xData,1)
    plot(xData{indTraj,1},xData{indTraj,2},'Linewidth',1)
end
xlabel('time [s]')
ylabel('velocity [m/s]')
```

![figure_0.png
](README_images/figure_0.png
)

```matlab:Code
xlim([xData{indTraj,1}(1) xData{indTraj,1}(end)])

% Train and test
indTest = [2 6 10];
indTrain = setdiff(1:size(xData,1),indTest);
```

Other than time inspection, we also look at time-frequency plots obtained via spectrogram. We can clearly see two frequency components, one which is more excited on some trajectories, less on others.

```matlab:Code
customFigure('subPlot',[2 1]); subplot(211); 
indTraj = 10;
numberWindows = round(length(xData{indTraj,1})/50); 
samplingTime = xData{indTraj,1}(2)-xData{indTraj,1}(1);
[Ss,frequencySpecg,timeSpecg]=spectrogram(xData{indTraj,2}(1,:),numberWindows,round(numberWindows*0.5),[],1/samplingTime);
[timeSpecgGrid,frequencySpecgGrid] = meshgrid(linspace(timeSpecg(1),timeSpecg(end),801),linspace(frequencySpecg(1),frequencySpecg(end),2401));
amplitudeSpecg = interp2(timeSpecg,frequencySpecg,abs(Ss),timeSpecgGrid,frequencySpecgGrid);
amplitudeSpecg(amplitudeSpecg<9) = NaN;
surf(timeSpecgGrid,frequencySpecgGrid,(amplitudeSpecg)); colormap(jet)
ylim([0 400])
xlim([timeSpecg(1) timeSpecg(end)])
xlabel('time [s]')
ylabel('frequency [Hz]')
set(gca,'fontname','helvetica')
set(gca,'fontsize',18)
set(gca,'zscale','log')
colorbar('northoutside')
view(2)
shading interp
climits = caxis;

subplot(212); 
indTraj = 12;
[Ss,frequencySpecg,timeSpecg]=spectrogram(xData{indTraj,2}(1,:),numberWindows,round(numberWindows*0.5),[],1/samplingTime);
[timeSpecgGrid,frequencySpecgGrid] = meshgrid(linspace(timeSpecg(1),timeSpecg(end),801),linspace(frequencySpecg(1),frequencySpecg(end),2401));
amplitudeSpecg = interp2(timeSpecg,frequencySpecg,abs(Ss),timeSpecgGrid,frequencySpecgGrid);
amplitudeSpecg(amplitudeSpecg<9) = NaN;
surf(timeSpecgGrid,frequencySpecgGrid,(amplitudeSpecg)); colormap(jet)
ylim([0 400])
xlim([timeSpecg(1) timeSpecg(end)])
xlabel('time [s]')
ylabel('frequency [Hz]')
set(gca,'fontname','helvetica')
set(gca,'fontsize',18)
set(gca,'zscale','log')
caxis(climits)
view(2)
shading interp
```

![figure_1.png
](README_images/figure_1.png
)

```matlab:Code

```

# Delay embedding

Now we arrange the scalar measurements in an observable space of dimension at least $2m+1$, with $m$ the dimension of the manifold. This guarantees that the manifold from the full state space can be embedded in the observable space. We form a multi-dimensional observable by stacking $d$ subsequent scalar measurements $x$ in a vector $y$, and we will later use the trajectories in this augmented space for the manifold fitting.

$$
\left\lbrack \begin{array}{c}
y_1 \\
\vdots \\
y_{i+1} \\
\vdots \\
y_{d+1} 
\end{array}\right\rbrack =\left\lbrack \begin{array}{c}
x(t^k )\\
\vdots \\
x(t^k +i\Delta t)\\
\vdots \\
x(t^k +d\Delta t)
\end{array}\right\rbrack
$$

The total number of embedded coordinates is $2m+1+{\textrm{overEmbed}}$. We impose additional delays so that an embedding vector covers 2 pseudo-periods of a slow oscillation.

```matlab:Code
samplingTime = xData{1,1}(2)-xData{1,1}(1);
slowTime = 1/120; % 1 divided by the (approximate) slow frequency
embeddingTime = 2 * slowTime;
SSMDim = 4;
overEmbed = round(embeddingTime/samplingTime);
[yData, opts_embd] = coordinatesEmbedding(xData, SSMDim, 'OverEmbedding', overEmbed);
```

```text:Output
The 94 embedding coordinates consist of the measured state and its 93 time-delayed measurements.
```

# Manifold fitting

The measured trajectories are assumed to lie close to a two-dimensional manifold that is tangent at the origin to the eigenspace corresponding to the slowest pair of eigenvalues. We now want to fit a polynomial of order $M$ to the data points to approximate the manifold. Here we use a graph style parametrization, meaning that the manifold is parametrized using coordinates in the computed eigenspace $V_e$. This excludes the possibility of folds on the manifold (however, if necessary we may circumvent this problem by increasing the dimensionality of the observed space). 

![image_1.png
](README_images/image_1.png
)

We seek the $2n\times m_M$ polynomial coefficient matrix $H$ of the manifold on the form

$y=V_e V_e^{\top } y+H\phi_{m,2:M} (V_e^{\top } y)$,

where the function $\phi_{m,2:M} (q)$ computes a vector of all $m_M$ monomials from orders 2 up to $M$ of an $m$-vector $q$. This optimization problem is formulated as the minimization of a cost function computed from the measured data $y_k$,

$$
C_r (V_e ,H)=\sum_{k=1}^N ||y_k -V_e V_e^{\top } y_k -H\phi_{m,2:M} (V_e^{\top } y_k )||
$$

to find $H$ and the $2n\times 2$ eigenspace matrix $V_e$, under the constraints that $V_e^{\top } H=0$ and $V_e^{\top } V_e =I$. The minimization is performed with the Matlab function `fmincon` in `IMGeometry`. In this case, we set $M=1$, so that our manifold is the plane spanned by the 4 most relevant directions of the Principal Component Anlysis of data.

```matlab:Code
SSMOrder = 1;
IMInfo = IMGeometry(yData(indTrain,:), SSMDim, SSMOrder);
```

# Plot and validation

Now that we have computed the tangent space of the manifold, we pass to the reduced coordinates $y$ by projecting all trajectories onto it:

$\eta =V_e^{\top } y$.

```matlab:Code
etaData = projectTrajectories(IMInfo, yData);
```

We plot the trajectories in the reduced coordinates.

```matlab:Code
plotReducedCoordinates(etaData);
view(-50,5); legend('off'); colororder(jet(size(xData,1)))
```

![figure_2.png
](README_images/figure_2.png
)

We also compute the static reconstruction error.

```matlab:Code
normedTrajDist = computeTrajectoryErrors(liftTrajectories(IMInfo,etaData), yData);
staticNMTE = mean(normedTrajDist(indTest))*100 % in percentage
```

```text:Output
staticNMTE = 
      0.16515

```

# Reduced order model

We compute a model for the reduced dynamics with the training data projected onto the manifold. The function `IMDynamicsFlow` fits a polynomial map

$$
\dot{\eta} =W_r \phi (\eta )
$$

where $\phi (\eta )$ again computes a vector of all monomials of $\eta$, and $W_r$ is a matrix of polynomial coefficients. 

We are also specifying that we want the reduced dynamics on a normal form, and seek the Taylor expansion of a map $N$ fulfilling

$$
\dot{z} =N(z)=Dz+W_n \phi (z)
$$

with $D$ a diagonal matrix and $W_n$ containing coefficients for the near-resonant nonlinear terms, after a near-identity change of coordinates

$$
z=T^{-1} (\eta )=\eta +W_t \phi (\eta )
$$

with $W_t$ containing the coefficients for the nonlinear terms of the coordinate change. We feed the algorithm with the specification of the ratio of the resonant frequencies to facilitate the recognition of the normal form.

```matlab:Code
ROMOrder = 3;
RDInfo = IMDynamicsFlow(etaData(indTrain,:), ...
    'R_PolyOrd', ROMOrder, 'style', 'normalform','frequencies_norm',[1 2],'IC_nf',1,'MaxIter',3e3);
```

```text:Output
Estimation of the reduced dynamics...  Done. 
Estimation of the reduced dynamics in normal form...
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
     0           1          0.27056                          6.19
     1           4         0.134393     0.00268634           8.63  
     2           6         0.121096       0.181728           2.83  
     3           7         0.108353              1           7.32  
     4           9        0.0985289       0.350751           1.73  
     5          10        0.0909289              1           1.01  
     6          11        0.0849516              1          0.955  
     7          12        0.0742036              1          0.777  
     8          13        0.0687933              1          0.898  
     9          14        0.0607143              1           2.18  
    10          16        0.0579362       0.223347           2.94  
    11          17        0.0503186              1          0.769  
    12          19        0.0486252       0.483655            2.3  
    13          20        0.0472985              1           1.79  
    14          21        0.0458836              1          0.215  
    15          22        0.0455654              1          0.202  
    16          23        0.0443543              1          0.656  
    17          24        0.0436788              1          0.538  
    18          25        0.0431721              1          0.155  
    19          26        0.0429441              1          0.128  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    20          27        0.0427482              1          0.217  
    21          28        0.0425197              1          0.218  
    22          29        0.0422245              1          0.169  
    23          30        0.0419212              1          0.134  
    24          31        0.0416478              1          0.149  
    25          32        0.0413859              1          0.161  
    26          33        0.0410709              1          0.173  
    27          34        0.0406415              1          0.218  
    28          35        0.0400899              1           0.18  
    29          36        0.0395712              1          0.207  
    30          37        0.0391068              1          0.179  
    31          38        0.0387328              1           0.17  
    32          39        0.0384905              1          0.131  
    33          40        0.0383618              1         0.0735  
    34          41        0.0382791              1         0.0602  
    35          42        0.0381939              1         0.0925  
    36          43        0.0380936              1          0.115  
    37          44         0.037999              1         0.0832  
    38          45        0.0379343              1         0.0659  
    39          46        0.0378924              1         0.0519  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    40          47        0.0378504              1         0.0622  
    41          48        0.0377905              1         0.0757  
    42          49        0.0377097              1          0.105  
    43          50         0.037617              1          0.126  
    44          51        0.0375177              1          0.152  
    45          52        0.0373976              1          0.145  
    46          53         0.037244              1          0.132  
    47          54        0.0370936              1          0.145  
    48          55         0.037007              1         0.0767  
    49          56        0.0369728              1         0.0367  
    50          57         0.036948              1          0.049  
    51          58        0.0369033              1          0.108  
    52          59        0.0368296              1          0.151  
    53          60        0.0367415              1          0.133  
    54          61        0.0366797              1         0.0549  
    55          62        0.0366467              1         0.0379  
    56          63        0.0366175              1         0.0701  
    57          64        0.0365705              1          0.112  
    58          65         0.036503              1          0.119  
    59          66        0.0364307              1         0.0712  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    60          67        0.0363738              1         0.0655  
    61          68        0.0363299              1         0.0543  
    62          69        0.0362881              1         0.0738  
    63          70        0.0362464              1         0.0528  
    64          71         0.036207              1         0.0507  
    65          72         0.036166              1         0.0539  
    66          73        0.0361169              1         0.0838  
    67          74        0.0360632              1         0.0803  
    68          75        0.0360177              1         0.0583  
    69          76        0.0359821              1         0.0485  
    70          77        0.0359431              1         0.0695  
    71          78        0.0358848              1          0.105  
    72          79        0.0358055              1           0.12  
    73          80        0.0357281              1          0.116  
    74          81        0.0356751              1         0.0601  
    75          82        0.0356374              1         0.0479  
    76          83        0.0355929              1         0.0805  
    77          84         0.035527              1          0.147  
    78          85        0.0354459              1          0.167  
    79          86        0.0353802              1           0.11  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    80          87        0.0353487              1         0.0334  
    81          88        0.0353358              1         0.0242  
    82          89        0.0353254              1         0.0449  
    83          90        0.0353089              1         0.0582  
    84          91        0.0352821              1         0.0683  
    85          92        0.0352452              1         0.0768  
    86          93        0.0352076              1         0.0612  
    87          94        0.0351786              1         0.0463  
    88          95        0.0351577              1          0.038  
    89          96        0.0351405              1         0.0452  
    90          97        0.0351256              1         0.0314  
    91          98        0.0351115              1         0.0386  
    92          99        0.0350941              1         0.0387  
    93         100        0.0350703              1          0.063  
    94         101        0.0350444              1         0.0646  
    95         102        0.0350264              1          0.037  
    96         103        0.0350163              1         0.0375  
    97         104        0.0350065              1         0.0348  
    98         105        0.0349897              1         0.0662  
    99         106        0.0349658              1         0.0791  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   100         107        0.0349456              1         0.0533  
   101         108        0.0349374              1         0.0159  
   102         109        0.0349348              1          0.016  
   103         110        0.0349321              1         0.0214  
   104         111        0.0349258              1         0.0412  
   105         112        0.0349136              1         0.0595  
   106         113        0.0348952              1         0.0621  
   107         114        0.0348788              1         0.0382  
   108         115         0.034871              1         0.0278  
   109         116        0.0348682              1         0.0147  
   110         117        0.0348664              1         0.0133  
   111         118        0.0348648              1        0.00927  
   112         119        0.0348636              1         0.0154  
   113         120        0.0348623              1         0.0164  
   114         121        0.0348606              1         0.0143  
   115         122        0.0348587              1         0.0146  
   116         123        0.0348574              1        0.00751  
   117         124        0.0348569              1        0.00647  
   118         125        0.0348564              1        0.00613  
   119         126        0.0348555              1         0.0121  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   120         127        0.0348538              1          0.018  
   121         128        0.0348515              1         0.0181  
   122         129        0.0348498              1        0.00977  
   123         130        0.0348491              1        0.00588  
   124         131        0.0348486              1        0.00536  
   125         132        0.0348478              1          0.011  
   126         133        0.0348461              1         0.0182  
   127         134        0.0348431              1         0.0228  
   128         135        0.0348397              1         0.0183  
   129         136        0.0348378              1        0.00686  
   130         137        0.0348372              1        0.00516  
   131         138        0.0348367              1        0.00579  
   132         139        0.0348356              1         0.0112  
   133         140        0.0348337              1         0.0156  
   134         141        0.0348312              1         0.0149  
   135         142        0.0348291              1         0.0141  
   136         143        0.0348281              1         0.0114  
   137         144        0.0348275              1        0.00702  
   138         145        0.0348269              1        0.00613  
   139         146        0.0348262              1        0.00496  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   140         147        0.0348258              1        0.00826  
   141         148        0.0348254              1        0.00875  
   142         149        0.0348248              1        0.00741  
   143         150         0.034824              1        0.00764  
   144         151        0.0348231              1         0.0061  
   145         152        0.0348225              1        0.00427  
   146         153        0.0348222              1        0.00302  
   147         154        0.0348219              1        0.00331  
   148         155        0.0348214              1        0.00561  
   149         156        0.0348205              1         0.0103  
   150         157        0.0348193              1         0.0129  
   151         158        0.0348179              1         0.0105  
   152         159        0.0348167              1         0.0121  
   153         160        0.0348155              1         0.0118  
   154         161        0.0348143              1         0.0109  
   155         162        0.0348132              1         0.0123  
   156         163        0.0348123              1        0.00772  
   157         164        0.0348119              1        0.00473  
   158         165        0.0348116              1        0.00286  
   159         166        0.0348114              1        0.00407  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   160         167        0.0348113              1        0.00323  
   161         168        0.0348111              1        0.00301  
   162         169        0.0348109              1         0.0045  
   163         170        0.0348106              1        0.00558  
   164         171        0.0348102              1        0.00734  
   165         172        0.0348097              1        0.00617  
   166         173        0.0348094              1        0.00305  
   167         174        0.0348093              1        0.00237  
   168         175        0.0348092              1        0.00209  
   169         176        0.0348091              1        0.00192  
   170         177         0.034809              1        0.00321  
   171         178        0.0348088              1        0.00357  
   172         179        0.0348086              1        0.00259  
   173         180        0.0348085              1        0.00215  
   174         181        0.0348085              1        0.00119  
   175         182        0.0348084              1        0.00161  
   176         183        0.0348083              1        0.00202  
   177         184        0.0348082              1         0.0038  
   178         185        0.0348081              1        0.00495  
   179         186        0.0348079              1        0.00499  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   180         187        0.0348077              1        0.00346  
   181         188        0.0348075              1        0.00314  
   182         189        0.0348074              1        0.00301  
   183         190        0.0348073              1        0.00254  
   184         191        0.0348073              1        0.00177  
   185         192        0.0348072              1        0.00171  
   186         193        0.0348071              1        0.00114  
   187         194        0.0348071              1        0.00117  
   188         195        0.0348071              1         0.0013  
   189         196        0.0348071              1        0.00129  
   190         197         0.034807              1        0.00205  
   191         198        0.0348069              1         0.0024  
   192         199        0.0348069              1        0.00164  
   193         200        0.0348068              1       0.000484  
   194         201        0.0348068              1       0.000257  
   195         202        0.0348068              1       0.000254  
   196         203        0.0348068              1       0.000248  
   197         204        0.0348068              1       0.000522  
   198         205        0.0348068              1       0.000934  
   199         206        0.0348068              1        0.00114  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   200         207        0.0348068              1       0.000861  
   201         208        0.0348068              1       0.000791  
   202         209        0.0348068              1       0.000754  
   203         210        0.0348068              1       0.000827  
   204         211        0.0348068              1         0.0011  
   205         212        0.0348068              1       0.000879  
   206         213        0.0348068              1       0.000484  
   207         214        0.0348068              1         0.0003  
   208         215        0.0348068              1       0.000385  
   209         216        0.0348068              1       0.000314  
   210         217        0.0348068              1       0.000147  
   211         218        0.0348068              1       0.000204  
   212         219        0.0348068              1        0.00051  
   213         220        0.0348068              1       0.000798  
   214         221        0.0348068              1        0.00087  
   215         222        0.0348068              1       0.000574  
   216         223        0.0348068              1       0.000393  
   217         224        0.0348068              1       0.000377  
   218         225        0.0348068              1       0.000499  
   219         226        0.0348067              1       0.000684  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   220         227        0.0348067              1       0.000668  
   221         228        0.0348067              1       0.000368  
   222         229        0.0348067              1       0.000332  
   223         230        0.0348067              1       0.000312  
   224         231        0.0348067              1       0.000371  
   225         232        0.0348067              1       0.000479  
   226         233        0.0348067              1       0.000451  
   227         234        0.0348067              1       0.000586  
   228         235        0.0348067              1        0.00055  
   229         236        0.0348067              1       0.000352  
   230         237        0.0348067              1       0.000186  
   231         238        0.0348067              1       0.000184  
   232         239        0.0348067              1       0.000238  
   233         240        0.0348067              1       0.000332  
   234         241        0.0348067              1       0.000409  
   235         242        0.0348067              1       0.000463  
   236         243        0.0348067              1       0.000473  
   237         244        0.0348067              1       0.000455  
   238         245        0.0348067              1       0.000412  
   239         246        0.0348067              1       0.000426  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   240         247        0.0348067              1       0.000392  
   241         248        0.0348067              1       0.000324  
   242         249        0.0348067              1       0.000347  
   243         250        0.0348067              1       0.000384  
   244         251        0.0348067              1       0.000233  
   245         252        0.0348067              1       8.76e-05  
   246         253        0.0348067              1       8.78e-05  
   247         254        0.0348067              1        8.7e-05  

Local minimum possible.

fminunc stopped because the size of the current step is less than
the value of the step size tolerance.

<stopping criteria details>
Plotting figure with the polar normal form equations ... Done. 
```

![figure_3.png
](README_images/figure_3.png
)

```matlab:Code
frequenciesHz = imag(RDInfo.eigenvaluesLinPartFlow(1:2))/2/pi
```

```text:Output
frequenciesHz = 2x1    
       122.38
        243.4

```

We transform the initial condition of our test trajectories according to the obtained change of coordinates, and integrate our reduced order evolution rule to predict the development of the trajectory. 

```matlab:Code
[yRec, etaRec, zRec] = advect(IMInfo, RDInfo, yData);
```

# Evaluation

The error NMTE is computed as the average distance of the predicted trajectory to the measured one in the full observable.

```matlab:Code
fullTrajDist = computeTrajectoryErrors(yRec, yData);
NMTE = mean(fullTrajDist(indTest))*100 % in percentage
```

```text:Output
RRMSE = 
       3.6498

```

We plot the measured and predicted signal for the test trajectories. 

```matlab:Code
plotTrajectories(yData(indTest,:), yRec(indTest,:), 'm', 1)
ylabel('$\hat{X} \, [\%]$','Interpreter','latex'); title('');
```

![figure_4.png
](README_images/figure_4.png
)

# Decays along the SSM and instantaneous damping

We plot the trajectories decay on the $(\rho_1 ,\rho_2 )$ plane, splitting hence the slow and fast directions.

```matlab:Code
customFigure; colororder(jet(size(xData,1)))
for indTraj = 1:size(xData,1)
    plot(abs(zRec{indTraj,2}(1,:)),abs(zRec{indTraj,2}(2,:)),'Linewidth',1.5)
end
xlabel('\rho_1')
ylabel('\rho_2')
set(gca,'xscale','log')
set(gca,'yscale','log')
```

![figure_5.png
](README_images/figure_5.png
)

We now focus on the decaying trajectories excited from the third impact position. First, we display the relative modal energy content over time:

```matlab:Code
customFigure; 
for indTraj = 9:size(xData,1)
    instEnergy = abs(zRec{indTraj,2}(1,:))+abs(zRec{indTraj,2}(2,:));
    yyaxis left
    plot(zRec{indTraj,1},abs(zRec{indTraj,2}(1,:))./instEnergy*100,'Linewidth',1.5)
    yyaxis right
    plot(zRec{indTraj,1},abs(zRec{indTraj,2}(2,:))./instEnergy*100,'Linewidth',1.5)
end
xlim([0 2])
xlabel('time')
yyaxis left
ylabel('\rho_1/(\rho_1+\rho_2) [%]')
yyaxis right
ylabel('\rho_2/(\rho_1+\rho_2) [%]')
```

![figure_6.png
](README_images/figure_6.png
)

and then we show the instantaneous damping evolution. Instantaneous damping is already a function handle present in the normal form reduced dynamics, so one just needs to be evaluated it and plot the results.

```matlab:Code
customFigure; 
for indTraj = 9:size(xData,1)
    instDamping = RDInfo.conjugateDynamics.damping(abs(zRec{indTraj,2}(1:2,:)),atan2(imag(zRec{indTraj,2}(1:2,:)),real(zRec{indTraj,2}(1:2,:))));    
    yyaxis left
    plot(zRec{indTraj,1},instDamping(1,:),'Linewidth',1.5)
    yyaxis right
    plot(zRec{indTraj,1},instDamping(2,:),'Linewidth',1.5)
end
xlim([0 2])
xlabel('time')
yyaxis left
ylabel('c_1')
yyaxis right
ylabel('c_2')
```

![figure_7.png
](README_images/figure_7.png
)

This is a preview of the livescript `brb.mlx`.

# Brake-Reuß Beam - DIC data

We consider digital image correlation measurements from a shaker ringdown test on a beam with a bolted lap joint [1,2]. Oscillations are initialized close to its slowest 2D spectral submanifold (SSM), which will be identified and used for model reduction.

The dataset contains the displacement field (206 points in 720 mm of beam length) of the bottom beam side measured via DIC and accelerations measured via accelerometers (ACC) in 2 locations. All data refers to a single experimental trial, acquired with different measurement devices (as described in [1]). The dataset also includes backbone curves results obtained using the Peak Finding and Fitting (PFF) method based on accelerometer signals. 

Here we construct a reduced-order model using DIC data and we validate it on the acceleration signals.

  

![image_0.png
](README_images/image_0.png
)

[1] W. Chen, D. Jana, A. Singh, M. Jin, M. Cenedese, G. Kosova, M. W. R. Brake, C. W. Schwingshackl, S. Nagarajaiah, K. Moore, and J. P. Noël. Measurement and identification of the nonlinear dynamics of a jointed structure using full-field data; Part I - Measurement of nonlinear dynamics. *Mechanical Systems and Signal Processing*, 166:108401, 2022. [https://doi.org/10.1016/j.ymssp.2021.108401](https://doi.org/10.1016/j.ymssp.2021.108401)

[2] M. Jin, G. Kosova, M. Cenedese, W. Chen, D. Jana, A. Singh, M. W. R. Brake, C. W. Schwingshackl, S. Nagarajaiah, K. Moore, and J. P. Noël. Measurement and identification of the nonlinear dynamics of a jointed structure using full-field data; Part II - Nonlinear system identification. *Mechanical Systems and Signal Processing*, 166:108402, 2022. [https://doi.org/10.1016/j.ymssp.2021.108402](https://doi.org/10.1016/j.ymssp.2021.108402)

```matlab:Code
clearvars
close all
```

# Experimental Data

We will divide the data into a training set, for model fitting, and a test set, for validation. In this case the experimental data consists of a single trajectory, with different measruements. We train on DIC data and validate on acceleration data.

```matlab:Code
load data.mat

data_BRB
```

```text:Output
data_BRB = 
            TimeDIC: [15260x1 double]
    DisplacementDIC: [206x15260 double]
            TimeACC: [1x272128 double]
    AccelerationACC: [2x272128 double]
              Xmesh: [1x206 double]
              Units: [1x1 struct]
        LocationACC: [2x1 double]
      PFFResultsACC: [1x1 struct]

```

```matlab:Code

% Store data
Xmesh = data_BRB.Xmesh;
% Displacements
uData = cell(1,2);
uData{1,1} = data_BRB.TimeDIC;
uData{1,2} = data_BRB.DisplacementDIC;
% Accelerations
aData{1,1} = data_BRB.TimeACC;
aData{1,2} = data_BRB.AccelerationACC;

% Plot data: acceleration & displacements
customFigure;
yyaxis left
plot(aData{1,1}, aData{1,2}(2,:))
ylabel('acceleration [m/s^2]'); ylim([-1 1]*max(abs(aData{1,2}(2,:))))
xlabel('$t \, [$s$]$','interpreter','latex');  
ylabel('$a \, [$m/s$^2]$','interpreter','latex'); 
yyaxis right
plot(uData{1,1}, uData{1,2}(183,:));  ylim([-1 1]*max(abs(uData{1,2}(183,:))))
ylabel('$u \, [$mm$]$','interpreter','latex'); 
```

![figure_0.png
](README_images/figure_0.png
)

```matlab:Code
customFigure;
map = pinkgreen(101);
[TT,XX] = meshgrid(uData{1,1},Xmesh);
surf(TT,XX,uData{1,2})
xlabel('$t \, [$s$]$','interpreter','latex'); 
ylabel('$x \, [$mm$]$','interpreter','latex'); 
zlabel('$u \, [$mm$]$','interpreter','latex'); 
shading  interp
colormap(map)
colorbar
xlim([0 1/80*3]+uData{1,1}(1))
ylim([Xmesh(1) Xmesh(end)])
view(3)
```

![figure_1.png
](README_images/figure_1.png
)

# Delay embedding

To diversify the data, we augment displacements with four delayed measurements, which will be also useful for recovering velocity and acceleration fields later.

```matlab:Code
overEmbed = 5;
SSMDim = 2;

yData = coordinatesEmbedding(uData, SSMDim, 'OverEmbedding', overEmbed);
```

```text:Output
The 1236 embedding coordinates consist of the 206 measured states and their 5 time-delayed measurements.
```

# Datadriven manifold fitting

The measured trajectories are assumed to lie close to a two-dimensional manifold that is tangent at the origin to the eigenspace corresponding to the slowest pair of eigenvalues. We now want to fit a polynomial of order $M$ to the data points to approximate the manifold. Here we use a graph style parametrization, meaning that the manifold is parametrized using coordinates in the computed eigenspace $V_e$. This excludes the possibility of folds on the manifold (however, if necessary we may circumvent this problem by increasing the dimensionality of the observed space). 

 ![image_1.png
](README_images/image_1.png
)

We seek the $2N\times m_M$ polynomial coefficient matrix $H$ of the manifold on the form

$y=V_e V_e^{\top } y+H\phi_{m,2:M} (V_e^{\top } y)$,

where the function $\phi_{m,2:M} (q)$ computes a vector of all $m_M$ monomials from orders 2 up to $M$ of an $m$-vector $q$. This optimization problem amounts to minimizing a cost function computed from the measured data $y_k$,

$$
C_r (V_e ,H)=\sum_{k=1}^N ||y_k -V_e V_e^{\top } y_k -H\phi_{m,2:M} (V_e^{\top } y_k )||
$$

to find $H$ and the $2N\times 2$ eigenspace matrix $V_e$, under the constraints that $V_e^{\top } H=0$ and $V_e^{\top } V_e =I$. The minimization is performed with the Matlab function `fmincon` in `IMparametrization`. In this case, we set $M=1$, so that our manifold is the plane spanned by the 2 most relevant directions of the Principal Component Anlysis of data.

```matlab:Code
SSMOrder = 1;
[IMInfo, SSMChart, SSMFunction] = IMGeometry(yData, SSMDim, SSMOrder);
```

# Plot and validation

Now that we have computed the eigenspace of the manifold, we pass to the reduced coordinates $y$ by projecting all trajectories onto the eigenspace. 

```matlab:Code
etaData = projectTrajectories(IMInfo, yData);
```

We plot the test and training set trajectories projected onto the plane $V$.

```matlab:Code
plotReducedCoordinates(etaData);
```

![figure_2.png
](README_images/figure_2.png
)

Furthermore, we draw the first component of the manifold shape along with the trajectory from the training set. 

```matlab:Code
outdof = 1;
plotSSMandTrajectories(IMInfo, outdof, yData, etaData)
view(-230,20); zlabel('$u_1 \, [$m$]$','Interpreter','latex')
```

![figure_3.png
](README_images/figure_3.png
)

# Reduced order model

We compute a model for the reduced dynamics with the truncated training data projected onto the manifold. The function `IMDynamicsFlow` fits a polynomial map

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

with $W_t$ containing the coefficients for the nonlinear terms of the coordinate change.

The optimization is performed in `IMDynamicsFlow` with the Matlab function `fminunc`. We apply a regularization on the polynomial coefficients to prevent overfitting. The regularization parameter is optimized from a logarithmic range of values with ridge regression. 

```matlab:Code
ROMOrder = 11;
Nfolds = 5;
Nregvals = 30; 
RDInfo = IMDynamicsFlow(etaData, 'R_PolyOrd', ROMOrder, 'n_folds', Nfolds, 'l_vals', logspace(-6,0,Nregvals), 'style', 'normalform');
```

```text:Output
Estimation of the reduced dynamics...  Done. 
Estimation of the reduced dynamics in normal form...
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
     0           1         0.286682                          16.3
     1           4        0.0959592     0.00122636           12.7  
     2           6        0.0829207       0.147846          0.977  
     3           7        0.0825095              1          0.979  
     4           8        0.0820495              1          0.492  
     5           9        0.0816302              1          0.489  
     6          10        0.0801643              1            1.3  
     7          12        0.0791611       0.352853           2.55  
     8          13        0.0745167              1           3.62  
     9          14        0.0670175              1           3.84  
    10          15        0.0447264              1           3.18  
    11          16         0.037179              1           1.81  
    12          17        0.0323693              1           1.03  
    13          18        0.0290042              1           1.31  
    14          19        0.0254562              1           1.58  
    15          20        0.0233247              1           0.89  
    16          21        0.0219112              1          0.546  
    17          22        0.0204096              1          0.851  
    18          23        0.0183494              1           1.42  
    19          24        0.0164742              1           1.23  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    20          25        0.0154975              1          0.424  
    21          26         0.015142              1          0.229  
    22          27        0.0149714              1          0.274  
    23          28         0.014773              1          0.271  
    24          29        0.0146472              1           0.14  
    25          30        0.0145767              1         0.0987  
    26          31        0.0145247              1          0.114  
    27          32        0.0144606              1          0.142  
    28          33        0.0144005              1          0.137  
    29          34        0.0143636              1         0.0694  
    30          35        0.0143431              1          0.075  
    31          36        0.0143191              1         0.0819  
    32          37        0.0142717              1          0.177  
    33          38        0.0141878              1          0.261  
    34          39        0.0140827              1          0.255  
    35          40         0.014011              1          0.131  
    36          41         0.013984              1         0.0736  
    37          42          0.01397              1         0.0664  
    38          43        0.0139474              1          0.115  
    39          44         0.013909              1          0.162  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    40          45        0.0138652              1          0.145  
    41          46        0.0138406              1         0.0639  
    42          47        0.0138341              1         0.0193  
    43          48        0.0138317              1         0.0192  
    44          49        0.0138275              1         0.0413  
    45          50        0.0138184              1         0.0698  
    46          51        0.0138019              1         0.0906  
    47          52        0.0137826              1         0.0778  
    48          53        0.0137714              1         0.0336  
    49          54        0.0137683              1         0.0193  
    50          55        0.0137672              1         0.0185  
    51          56        0.0137651              1         0.0219  
    52          57        0.0137603              1          0.038  
    53          58        0.0137509              1         0.0529  
    54          59        0.0137377              1          0.056  
    55          60        0.0137278              1         0.0336  
    56          61        0.0137243              1         0.0137  
    57          62        0.0137232              1         0.0135  
    58          63        0.0137218              1         0.0147  
    59          64        0.0137183              1         0.0309  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    60          65        0.0137103              1         0.0524  
    61          66        0.0136942              1         0.0736  
    62          67        0.0136716              1         0.0745  
    63          68        0.0136545              1         0.0425  
    64          69        0.0136489              1         0.0144  
    65          70        0.0136475              1         0.0143  
    66          71        0.0136462              1         0.0141  
    67          72        0.0136426              1         0.0268  
    68          73        0.0136341              1         0.0461  
    69          74        0.0136134              1         0.0725  
    70          75        0.0135688              1            0.1  
    71          76        0.0134854              1          0.136  
    72          77        0.0133662              1           0.19  
    73          78        0.0132465              1          0.196  
    74          79        0.0131601              1          0.142  
    75          80        0.0131065              1         0.0687  
    76          81        0.0130757              1         0.0501  
    77          82        0.0130621              1         0.0312  
    78          83         0.013056              1         0.0285  
    79          84        0.0130495              1          0.029  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
    80          85        0.0130376              1         0.0384  
    81          86         0.013019              1         0.0442  
    82          87        0.0129977              1         0.0565  
    83          88        0.0129795              1         0.0816  
    84          89        0.0129616              1         0.0953  
    85          90         0.012933              1          0.104  
    86          91          0.01288              1          0.109  
    87          92        0.0127955              1          0.106  
    88          93        0.0126966              1          0.169  
    89          94        0.0126111              1          0.206  
    90          95        0.0125394              1          0.187  
    91          96        0.0124683              1           0.12  
    92          97        0.0124078              1          0.058  
    93          98        0.0123767              1         0.0607  
    94          99        0.0123641              1         0.0635  
    95         100        0.0123514              1         0.0653  
    96         101        0.0123208              1         0.0677  
    97         102        0.0122471              1          0.111  
    98         103         0.012073              1          0.176  
    99         104        0.0117254              1          0.237  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   100         105        0.0112255              1          0.235  
   101         106        0.0108242              1           0.14  
   102         107        0.0106653              1          0.152  
   103         108        0.0106063              1          0.144  
   104         109        0.0105319              1          0.128  
   105         110        0.0103824              1           0.11  
   106         111         0.010171              1           0.12  
   107         112         0.010003              1         0.0738  
   108         113       0.00994677              1          0.047  
   109         114       0.00993589              1         0.0421  
   110         115       0.00992861              1         0.0359  
   111         116       0.00990936              1         0.0302  
   112         117       0.00987372              1         0.0442  
   113         118        0.0098197              1         0.0471  
   114         119       0.00977289              1         0.0695  
   115         120       0.00974778              1         0.0725  
   116         121       0.00973101              1         0.0656  
   117         122       0.00970103              1         0.0705  
   118         123       0.00963086              1         0.0806  
   119         124       0.00947627              1          0.094  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   120         125       0.00920548              1          0.116  
   121         126        0.0089104              1         0.0932  
   122         127       0.00876251              1         0.0372  
   123         128       0.00873398              1         0.0303  
   124         129       0.00872944              1          0.029  
   125         130       0.00872374              1         0.0277  
   126         131       0.00870745              1         0.0251  
   127         132       0.00866908              1         0.0288  
   128         133       0.00858042              1         0.0473  
   129         134       0.00841663              1          0.062  
   130         135       0.00821194              1          0.083  
   131         136       0.00807826              1          0.112  
   132         137       0.00803149              1          0.115  
   133         138       0.00800707              1           0.11  
   134         139       0.00795979              1         0.0971  
   135         140       0.00785835              1          0.072  
   136         141         0.007677              1         0.0539  
   137         142       0.00747366              1         0.0454  
   138         143       0.00736589              1         0.0303  
   139         144       0.00734267              1          0.029  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   140         145       0.00733762              1         0.0284  
   141         146        0.0073298              1         0.0277  
   142         147       0.00730866              1         0.0265  
   143         148       0.00725583              1         0.0268  
   144         149       0.00712471              1          0.047  
   145         150       0.00683215              1         0.0799  
   146         151       0.00630668              1          0.105  
   147         152       0.00570917              1         0.0902  
   148         153       0.00539069              1         0.0415  
   149         154       0.00532614              1         0.0423  
   150         155       0.00531743              1         0.0408  
   151         156       0.00530932              1         0.0389  
   152         157       0.00528515              1         0.0333  
   153         158       0.00523456              1         0.0216  
   154         159       0.00513827              1         0.0272  
   155         160       0.00502323              1         0.0283  
   156         161       0.00495328              1         0.0434  
   157         162       0.00493393              1         0.0426  
   158         163       0.00492725              1         0.0393  
   159         164       0.00491513              1         0.0343  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   160         165       0.00488425              1         0.0259  
   161         166        0.0048072              1          0.027  
   162         167       0.00462643              1         0.0363  
   163         168       0.00426692              1         0.0483  
   164         169       0.00376584              1         0.0473  
   165         170       0.00340087              1         0.0343  
   166         171       0.00329633              1         0.0201  
   167         172       0.00328489              1         0.0198  
   168         173       0.00328235              1         0.0197  
   169         174       0.00327519              1         0.0195  
   170         175       0.00325837              1         0.0192  
   171         176       0.00321457              1         0.0185  
   172         177       0.00311264              1         0.0227  
   173         178       0.00290537              1         0.0328  
   174         179       0.00260702              1         0.0348  
   175         180       0.00237587              1           0.04  
   176         181       0.00230314              1         0.0399  
   177         182       0.00229322              1         0.0369  
   178         183       0.00228972              1         0.0354  
   179         184       0.00227924              1         0.0323  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   180         185       0.00225526              1         0.0275  
   181         186       0.00219515              1         0.0192  
   182         187       0.00206918              1         0.0152  
   183         188       0.00186584              1          0.018  
   184         189       0.00167783              1         0.0253  
   185         190       0.00160308              1         0.0264  
   186         191        0.0015919              1         0.0239  
   187         192       0.00158998              1         0.0226  
   188         193       0.00158659              1          0.021  
   189         194        0.0015779              1         0.0178  
   190         195        0.0015578              1         0.0138  
   191         196       0.00151771              1          0.013  
   192         197       0.00146209              1         0.0113  
   193         198       0.00142147              1         0.0172  
   194         199       0.00140947              1         0.0156  
   195         200        0.0014077              1         0.0133  
   196         201       0.00140677              1         0.0121  
   197         202       0.00140382              1         0.0103  
   198         203       0.00139701              1         0.0107  
   199         204       0.00137949              1         0.0112  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   200         205       0.00134046              1         0.0123  
   201         206       0.00126838              1         0.0255  
   202         207       0.00118288              1         0.0332  
   203         208       0.00113415              1         0.0284  
   204         209       0.00112358              1         0.0207  
   205         210       0.00112224              1         0.0179  
   206         211       0.00112124              1          0.017  
   207         212       0.00111805              1         0.0172  
   208         213       0.00111121              1         0.0164  
   209         214       0.00109676              1         0.0134  
   210         215       0.00107607              1        0.00861  
   211         216       0.00105972              1         0.0104  
   212         217       0.00105453              1        0.00675  
   213         218       0.00105384              1        0.00704  
   214         219       0.00105363              1        0.00707  
   215         220       0.00105302              1        0.00711  
   216         221       0.00105159              1        0.00713  
   217         222       0.00104779              1        0.00708  
   218         223       0.00103848              1        0.00991  
   219         224       0.00101707              1         0.0172  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   220         225      0.000976943              1         0.0235  
   221         226      0.000927071              1         0.0217  
   222         227      0.000896968              1         0.0106  
   223         228      0.000890266              1        0.00457  
   224         229      0.000889764              1        0.00436  
   225         230      0.000889709              1        0.00429  
   226         231      0.000889575              1        0.00418  
   227         232      0.000889243              1        0.00399  
   228         233      0.000888366              1        0.00367  
   229         234      0.000886149              1         0.0031  
   230         235      0.000880723              1        0.00393  
   231         236      0.000868833              1        0.00611  
   232         237      0.000848443              1        0.00758  
   233         238      0.000827151              1        0.00711  
   234         239      0.000817177              1        0.00816  
   235         240      0.000815456              1        0.00789  
   236         241      0.000815249              1        0.00765  
   237         242      0.000815037              1        0.00743  
   238         243      0.000814381              1         0.0069  
   239         244      0.000812875              1        0.00594  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   240         245      0.000809352              1        0.00423  
   241         246      0.000802899              1        0.00584  
   242         247      0.000795051              1        0.00555  
   243         248      0.000790481              1        0.00358  
   244         249      0.000789469              1        0.00304  
   245         250      0.000789355              1        0.00309  
   246         251      0.000789287              1         0.0031  
   247         252      0.000789051              1        0.00311  
   248         253      0.000788503              1        0.00312  
   249         254       0.00078701              1        0.00443  
   250         255       0.00078322              1        0.00804  
   251         256      0.000773581              1         0.0139  
   252         257      0.000750609              1         0.0224  
   253         258      0.000702667              1         0.0317  
   254         259      0.000629136              1         0.0344  
   255         260       0.00056626              1         0.0227  
   256         261      0.000544268              1        0.00707  
   257         262      0.000541734              1        0.00166  
   258         263      0.000541635              1        0.00164  
   259         264      0.000541618              1        0.00163  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   260         265       0.00054154              1        0.00161  
   261         266      0.000541378              1        0.00159  
   262         267      0.000540945              1         0.0033  
   263         268      0.000539999              1        0.00588  
   264         269      0.000538242              1        0.00843  
   265         270      0.000536156              1        0.00823  
   266         271      0.000534944              1        0.00441  
   267         272      0.000534663              1        0.00355  
   268         273      0.000534611              1        0.00339  
   269         274       0.00053455              1        0.00323  
   270         275      0.000534377              1        0.00327  
   271         276      0.000533951              1        0.00669  
   272         277      0.000532851              1         0.0125  
   273         278       0.00053019              1         0.0212  
   274         279       0.00052434              1          0.032  
   275         280       0.00051431              1         0.0388  
   276         281      0.000503783              1         0.0307  
   277         282      0.000498818              1         0.0122  
   278         283      0.000497955              1        0.00562  
   279         284      0.000497854              1        0.00556  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   280         285      0.000497757              1        0.00548  
   281         286      0.000497452              1        0.00527  
   282         287      0.000496741              1         0.0092  
   283         288      0.000494973              1          0.016  
   284         289      0.000491236              1         0.0241  
   285         290      0.000485124              1         0.0283  
   286         291      0.000479321              1         0.0209  
   287         292       0.00047695              1        0.00748  
   288         293      0.000476616              1         0.0016  
   289         294      0.000476592              1         0.0016  
   290         295      0.000476579              1        0.00159  
   291         296      0.000476529              1        0.00158  
   292         297      0.000476419              1         0.0032  
   293         298      0.000476126              1        0.00611  
   294         299      0.000475439              1         0.0103  
   295         300       0.00047397              1         0.0154  
   296         301      0.000471628              1         0.0176  
   297         302      0.000469454              1         0.0126  
   298         303      0.000468588              1        0.00432  
   299         304      0.000468453              1        0.00304  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   300         305      0.000468425              1        0.00293  
   301         306      0.000468371              1        0.00276  
   302         307      0.000468232              1         0.0035  
   303         308      0.000467879              1        0.00653  
   304         309      0.000467017              1         0.0111  
   305         310      0.000465104              1         0.0167  
   306         311       0.00046174              1         0.0204  
   307         312      0.000458041              1         0.0165  
   308         313      0.000456163              1        0.00668  
   309         314      0.000455796              1        0.00319  
   310         315      0.000455739              1        0.00296  
   311         316      0.000455671              1        0.00306  
   312         317      0.000455474              1        0.00341  
   313         318      0.000455004              1        0.00611  
   314         319      0.000453863              1           0.01  
   315         320      0.000451532              1         0.0141  
   316         321      0.000448029              1         0.0149  
   317         322      0.000445161              1        0.00939  
   318         323      0.000444211              1        0.00263  
   319         324      0.000444107              1       0.000894  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   320         325      0.000444102              1       0.000882  
   321         326      0.000444099              1       0.000877  
   322         327      0.000444088              1       0.000861  
   323         328      0.000444063              1       0.000839  
   324         329      0.000443995              1        0.00126  
   325         330      0.000443827              1        0.00221  
   326         331      0.000443417              1        0.00364  
   327         332      0.000442548              1        0.00526  
   328         333      0.000441149              1        0.00592  
   329         334      0.000439851              1        0.00411  
   330         335      0.000439339              1         0.0015  
   331         336      0.000439268              1        0.00146  
   332         337      0.000439261              1        0.00142  
   333         338      0.000439254              1        0.00139  
   334         339      0.000439232              1        0.00133  
   335         340      0.000439178              1        0.00123  
   336         341      0.000439042              1        0.00184  
   337         342      0.000438734              1        0.00282  
   338         343      0.000438152              1        0.00364  
   339         344       0.00043743              1        0.00329  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   340         345      0.000436994              1        0.00159  
   341         346      0.000436896              1       0.000432  
   342         347      0.000436889              1       0.000377  
   343         348      0.000436888              1       0.000376  
   344         349      0.000436887              1       0.000373  
   345         350      0.000436882              1       0.000368  
   346         351      0.000436871              1       0.000466  
   347         352      0.000436843              1       0.000824  
   348         353      0.000436776              1        0.00137  
   349         354       0.00043663              1        0.00202  
   350         355      0.000436392              1        0.00233  
   351         356      0.000436161              1        0.00167  
   352         357      0.000436064              1        0.00104  
   353         358      0.000436049              1       0.000817  
   354         359      0.000436047              1       0.000754  
   355         360      0.000436043              1       0.000689  
   356         361      0.000436033              1       0.000617  
   357         362      0.000436008              1       0.000971  
   358         363      0.000435942              1        0.00166  
   359         364      0.000435777              1        0.00275  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   360         365       0.00043539              1        0.00426  
   361         366      0.000434617              1        0.00574  
   362         367      0.000433531              1        0.00571  
   363         368       0.00043273              1        0.00328  
   364         369        0.0004325              1       0.000823  
   365         370      0.000432478              1       0.000515  
   366         371      0.000432477              1       0.000513  
   367         372      0.000432476              1        0.00051  
   368         373      0.000432473              1       0.000502  
   369         374      0.000432465              1       0.000486  
   370         375      0.000432443              1       0.000749  
   371         376      0.000432391              1        0.00133  
   372         377      0.000432276              1         0.0021  
   373         378      0.000432076              1        0.00265  
   374         379      0.000431863              1        0.00221  
   375         380       0.00043176              1        0.00095  
   376         381      0.000431742              1       0.000684  
   377         382      0.000431741              1       0.000651  
   378         383       0.00043174              1       0.000634  
   379         384      0.000431736              1       0.000595  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   380         385      0.000431728              1       0.000615  
   381         386      0.000431705              1        0.00121  
   382         387      0.000431649              1        0.00218  
   383         388      0.000431509              1        0.00371  
   384         389       0.00043119              1        0.00573  
   385         390      0.000430589              1        0.00743  
   386         391      0.000429844              1        0.00674  
   387         392      0.000429394              1         0.0033  
   388         393      0.000429294              1       0.000678  
   389         394      0.000429287              1       0.000429  
   390         395      0.000429286              1       0.000426  
   391         396      0.000429285              1       0.000422  
   392         397      0.000429283              1       0.000414  
   393         398      0.000429276              1       0.000541  
   394         399      0.000429259              1        0.00104  
   395         400      0.000429218              1         0.0018  
   396         401      0.000429129              1        0.00278  
   397         402      0.000428976              1        0.00341  
   398         403      0.000428816              1        0.00274  
   399         404      0.000428741              1        0.00111  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   400         405      0.000428729              1       0.000255  
   401         406      0.000428729              1       0.000236  
   402         407      0.000428729              1       0.000231  
   403         408      0.000428728              1       0.000219  
   404         409      0.000428726              1       0.000302  
   405         410      0.000428722              1       0.000622  
   406         411      0.000428711              1        0.00115  
   407         412      0.000428684              1          0.002  
   408         413      0.000428618              1        0.00325  
   409         414      0.000428478              1        0.00468  
   410         415      0.000428253              1        0.00526  
   411         416      0.000428043              1        0.00369  
   412         417       0.00042796              1        0.00122  
   413         418      0.000427949              1       0.000514  
   414         419      0.000427948              1        0.00051  
   415         420      0.000427948              1       0.000505  
   416         421      0.000427945              1       0.000495  
   417         422       0.00042794              1       0.000539  
   418         423      0.000427925              1       0.000999  
   419         424      0.000427887              1        0.00172  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   420         425      0.000427795              1         0.0028  
   421         426      0.000427594              1        0.00408  
   422         427      0.000427253              1        0.00476  
   423         428      0.000426903              1         0.0036  
   424         429      0.000426745              1        0.00133  
   425         430       0.00042672              1       0.000199  
   426         431      0.000426719              1       0.000197  
   427         432      0.000426719              1       0.000196  
   428         433      0.000426718              1       0.000193  
   429         434      0.000426717              1        0.00019  
   430         435      0.000426715              1       0.000185  
   431         436      0.000426709              1       0.000269  
   432         437      0.000426694              1        0.00049  
   433         438      0.000426653              1       0.000859  
   434         439      0.000426549              1        0.00146  
   435         440       0.00042629              1        0.00237  
   436         441      0.000425693              1        0.00357  
   437         442      0.000424553              1        0.00456  
   438         443      0.000423092              1        0.00414  
   439         444      0.000422163              1        0.00204  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   440         445      0.000421944              1       0.000403  
   441         446      0.000421928              1       0.000216  
   442         447      0.000421928              1       0.000215  
   443         448      0.000421927              1       0.000215  
   444         449      0.000421927              1       0.000213  
   445         450      0.000421925              1       0.000211  
   446         451       0.00042192              1       0.000208  
   447         452      0.000421907              1       0.000201  
   448         453      0.000421873              1       0.000324  
   449         454      0.000421788              1       0.000566  
   450         455      0.000421579              1       0.000931  
   451         456      0.000421121              1        0.00137  
   452         457      0.000420338              1        0.00163  
   453         458      0.000419525              1        0.00125  
   454         459      0.000419151              1       0.000478  
   455         460      0.000419092              1       0.000207  
   456         461      0.000419089              1       0.000201  
   457         462      0.000419089              1         0.0002  
   458         463      0.000419089              1       0.000198  
   459         464      0.000419088              1       0.000196  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   460         465      0.000419085              1       0.000191  
   461         466      0.000419078              1       0.000182  
   462         467       0.00041906              1       0.000163  
   463         468      0.000419018              1       0.000271  
   464         469      0.000418925              1       0.000407  
   465         470      0.000418768              1       0.000485  
   466         471      0.000418608              1       0.000373  
   467         472      0.000418536              1       0.000212  
   468         473      0.000418525              1        0.00019  
   469         474      0.000418525              1        0.00018  
   470         475      0.000418525              1       0.000177  
   471         476      0.000418524              1        0.00017  
   472         477      0.000418523              1       0.000161  
   473         478      0.000418519              1       0.000159  
   474         479      0.000418508              1       0.000158  
   475         480      0.000418482              1       0.000195  
   476         481      0.000418414              1       0.000332  
   477         482      0.000418244              1       0.000542  
   478         483      0.000417846              1       0.000826  
   479         484      0.000417062              1        0.00108  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   480         485      0.000415989              1        0.00104  
   481         486      0.000415234              1       0.000559  
   482         487       0.00041503              1       0.000238  
   483         488      0.000415012              1       0.000231  
   484         489      0.000415012              1       0.000228  
   485         490      0.000415012              1       0.000227  
   486         491      0.000415011              1       0.000225  
   487         492      0.000415009              1       0.000221  
   488         493      0.000415003              1       0.000216  
   489         494      0.000414988              1       0.000206  
   490         495       0.00041495              1       0.000192  
   491         496      0.000414854              1       0.000169  
   492         497      0.000414626              1       0.000219  
   493         498      0.000414158              1       0.000305  
   494         499      0.000413462              1       0.000342  
   495         500        0.0004129              1       0.000497  
   496         501      0.000412717              1       0.000523  
   497         502      0.000412698              1         0.0005  
   498         503      0.000412696              1       0.000491  
   499         504      0.000412696              1       0.000486  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   500         505      0.000412692              1       0.000475  
   501         506      0.000412685              1       0.000459  
   502         507      0.000412664              1       0.000429  
   503         508      0.000412611              1        0.00038  
   504         509      0.000412475              1       0.000291  
   505         510      0.000412139              1       0.000225  
   506         511      0.000411372              1       0.000292  
   507         512      0.000409947              1       0.000495  
   508         513      0.000408211              1       0.000797  
   509         514      0.000407195              1       0.000796  
   510         515      0.000406977              1       0.000646  
   511         516      0.000406962              1       0.000583  
   512         517      0.000406961              1       0.000573  
   513         518      0.000406958              1       0.000562  
   514         519      0.000406952              1       0.000542  
   515         520      0.000406934              1        0.00051  
   516         521      0.000406888              1       0.000469  
   517         522       0.00040677              1       0.000445  
   518         523      0.000406475              1       0.000398  
   519         524      0.000405796              1       0.000311  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   520         525      0.000404494              1       0.000266  
   521         526      0.000402815              1       0.000418  
   522         527      0.000401738              1       0.000296  
   523         528       0.00040148              1       0.000147  
   524         529      0.000401461              1       0.000147  
   525         530      0.000401461              1       0.000147  
   526         531      0.000401461              1       0.000147  
   527         532       0.00040146              1       0.000147  
   528         533      0.000401458              1       0.000146  
   529         534      0.000401453              1       0.000145  
   530         535      0.000401439              1       0.000143  
   531         536      0.000401404              1       0.000141  
   532         537      0.000401314              1       0.000136  
   533         538       0.00040109              1       0.000126  
   534         539      0.000400571              1       0.000185  
   535         540      0.000399572              1        0.00046  
   536         541      0.000398267              1       0.000806  
   537         542      0.000397413              1       0.000996  
   538         543      0.000397203              1       0.000972  
   539         544      0.000397186              1       0.000927  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   540         545      0.000397185              1       0.000917  
   541         546      0.000397183              1       0.000907  
   542         547      0.000397179              1       0.000889  
   543         548      0.000397166              1       0.000861  
   544         549      0.000397133              1       0.000812  
   545         550      0.000397049              1        0.00073  
   546         551      0.000396833              1       0.000591  
   547         552      0.000396312              1       0.000361  
   548         553      0.000395189              1       0.000105  
   549         554       0.00039334              1       0.000345  
   550         555      0.000391551              1       0.000458  
   551         556      0.000390807              1       0.000257  
   552         557      0.000390702              1        0.00012  
   553         558      0.000390698              1       0.000118  
   554         559      0.000390698              1       0.000118  
   555         560      0.000390697              1       0.000118  
   556         561      0.000390697              1       0.000117  
   557         562      0.000390696              1       0.000116  
   558         563      0.000390692              1       0.000114  
   559         564      0.000390682              1       0.000111  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   560         565      0.000390657              1       0.000106  
   561         566      0.000390592              1       9.75e-05  
   562         567      0.000390429              1       9.56e-05  
   563         568      0.000390043              1       0.000134  
   564         569      0.000389261              1       0.000214  
   565         570      0.000388129              1       0.000396  
   566         571      0.000387255              1       0.000572  
   567         572      0.000386989              1       0.000603  
   568         573      0.000386962              1        0.00058  
   569         574      0.000386961              1       0.000572  
   570         575       0.00038696              1       0.000568  
   571         576      0.000386957              1       0.000559  
   572         577       0.00038695              1       0.000545  
   573         578       0.00038693              1       0.000522  
   574         579      0.000386881              1       0.000483  
   575         580      0.000386752              1       0.000418  
   576         581       0.00038643              1       0.000309  
   577         582      0.000385668              1       0.000134  
   578         583      0.000384122              1       0.000104  
   579         584       0.00038188              1         0.0003  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   580         585      0.000380141              1       0.000279  
   581         586      0.000379608              1       0.000117  
   582         587      0.000379555              1       0.000102  
   583         588      0.000379554              1       0.000101  
   584         589      0.000379554              1       0.000101  
   585         590      0.000379553              1         0.0001  
   586         591      0.000379553              1       9.92e-05  
   587         592      0.000379551              1        9.8e-05  
   588         593      0.000379548              1       9.59e-05  
   589         594      0.000379538              1       9.23e-05  
   590         595      0.000379512              1       8.65e-05  
   591         596      0.000379444              1       7.89e-05  
   592         597      0.000379271              1       7.37e-05  
   593         598      0.000378842              1       6.11e-05  
   594         599      0.000377868              1       6.28e-05  
   595         600      0.000376064              1        9.4e-05  
   596         601      0.000373887              1       0.000166  
   597         602       0.00037263              1       0.000204  
   598         603      0.000372367              1       0.000201  
   599         604       0.00037235              1       0.000194  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   600         605      0.000372349              1       0.000193  
   601         606      0.000372349              1       0.000193  
   602         607      0.000372348              1       0.000191  
   603         608      0.000372347              1        0.00019  
   604         609      0.000372342              1       0.000186  
   605         610       0.00037233              1        0.00018  
   606         611      0.000372299              1       0.000168  
   607         612      0.000372219              1       0.000144  
   608         613       0.00037202              1       9.11e-05  
   609         614      0.000371561              1       7.58e-05  
   610         615      0.000370683              1       0.000209  
   611         616      0.000369553              1       0.000446  
   612         617      0.000368832              1       0.000569  
   613         618       0.00036866              1       0.000548  
   614         619      0.000368647              1       0.000518  
   615         620      0.000368647              1       0.000512  
   616         621      0.000368646              1       0.000509  
   617         622      0.000368644              1         0.0005  
   618         623       0.00036864              1       0.000489  
   619         624      0.000368628              1       0.000469  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   620         625      0.000368597              1       0.000437  
   621         626      0.000368516              1       0.000384  
   622         627      0.000368305              1       0.000298  
   623         628      0.000367764              1       0.000158  
   624         629      0.000366416              1       0.000106  
   625         630      0.000363306              1       0.000392  
   626         631      0.000357345              1       0.000753  
   627         632      0.000349645              1       0.000885  
   628         633        0.0003447              1       0.000601  
   629         634      0.000343517              1       0.000359  
   630         635      0.000343429              1       0.000373  
   631         636      0.000343427              1       0.000372  
   632         637      0.000343426              1       0.000371  
   633         638      0.000343424              1       0.000369  
   634         639      0.000343421              1       0.000366  
   635         640      0.000343409              1       0.000361  
   636         641      0.000343381              1       0.000351  
   637         642      0.000343306              1       0.000331  
   638         643      0.000343116              1        0.00029  
   639         644      0.000342658              1       0.000207  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   640         645      0.000341674              1       0.000161  
   641         646      0.000340069              1       0.000188  
   642         647       0.00033854              1       0.000389  
   643         648      0.000337919              1        0.00042  
   644         649      0.000337833              1       0.000375  
   645         650       0.00033783              1       0.000358  
   646         651      0.000337829              1       0.000356  
   647         652      0.000337828              1       0.000351  
   648         653      0.000337826              1       0.000344  
   649         654      0.000337819              1       0.000333  
   650         655      0.000337803              1       0.000314  
   651         656      0.000337758              1       0.000284  
   652         657      0.000337643              1       0.000235  
   653         658      0.000337344              1       0.000216  
   654         659      0.000336575              1       0.000204  
   655         660      0.000334662              1       0.000178  
   656         661      0.000330262              1       0.000456  
   657         662       0.00032187              1       0.000738  
   658         663      0.000311149              1       0.000778  
   659         664      0.000304381              1       0.000456  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   660         665      0.000302797              1       0.000328  
   661         666      0.000302682              1       0.000333  
   662         667      0.000302679              1       0.000331  
   663         668      0.000302679              1        0.00033  
   664         669      0.000302677              1       0.000329  
   665         670      0.000302674              1       0.000327  
   666         671      0.000302665              1       0.000323  
   667         672      0.000302642              1       0.000315  
   668         673      0.000302581              1       0.000301  
   669         674      0.000302426              1       0.000271  
   670         675      0.000302039              1        0.00021  
   671         676      0.000301146              1       0.000166  
   672         677      0.000299428              1       0.000159  
   673         678        0.0002972              1       0.000414  
   674         679      0.000295757              1       0.000557  
   675         680      0.000295409              1       0.000531  
   676         681      0.000295382              1       0.000495  
   677         682      0.000295381              1       0.000488  
   678         683      0.000295381              1       0.000485  
   679         684      0.000295378              1       0.000478  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   680         685      0.000295373              1       0.000468  
   681         686       0.00029536              1       0.000451  
   682         687      0.000295325              1       0.000424  
   683         688      0.000295232              1       0.000379  
   684         689      0.000294991              1       0.000306  
   685         690      0.000294367              1       0.000187  
   686         691       0.00029278              1       0.000152  
   687         692      0.000288912              1       0.000296  
   688         693      0.000280471              1       0.000678  
   689         694      0.000266173              1       0.000983  
   690         695      0.000251599              1       0.000879  
   691         696      0.000245068              1       0.000622  
   692         697      0.000244065              1       0.000644  
   693         698      0.000244019              1       0.000635  
   694         699      0.000244018              1       0.000633  
   695         700      0.000244017              1       0.000632  
   696         701      0.000244014              1       0.000629  
   697         702      0.000244006              1       0.000625  
   698         703      0.000243985              1       0.000617  
   699         704      0.000243932              1       0.000604  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   700         705      0.000243792              1       0.000579  
   701         706      0.000243436              1       0.000534  
   702         707       0.00024256              1       0.000448  
   703         708      0.000240609              1       0.000292  
   704         709      0.000237154              1       0.000127  
   705         710      0.000233324              1       0.000142  
   706         711      0.000231383              1       0.000172  
   707         712      0.000231037              1       0.000156  
   708         713      0.000231018              1       0.000157  
   709         714      0.000231018              1       0.000157  
   710         715      0.000231018              1       0.000157  
   711         716      0.000231018              1       0.000157  
   712         717      0.000231017              1       0.000156  
   713         718      0.000231014              1       0.000156  
   714         719      0.000231007              1       0.000154  
   715         720      0.000230988              1       0.000152  
   716         721       0.00023094              1       0.000149  
   717         722      0.000230814              1       0.000143  
   718         723      0.000230495              1       0.000134  
   719         724      0.000229713              1       0.000117  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   720         725      0.000227987              1       0.000184  
   721         726       0.00022499              1       0.000271  
   722         727      0.000221793              1       0.000293  
   723         728      0.000220261              1       0.000371  
   724         729      0.000220006              1       0.000372  
   725         730      0.000219993              1       0.000364  
   726         731      0.000219993              1       0.000362  
   727         732      0.000219993              1       0.000361  
   728         733      0.000219991              1       0.000359  
   729         734      0.000219988              1       0.000356  
   730         735      0.000219979              1        0.00035  
   731         736      0.000219957              1        0.00034  
   732         737      0.000219899              1       0.000324  
   733         738      0.000219749              1       0.000297  
   734         739      0.000219368              1        0.00025  
   735         740      0.000218452              1        0.00017  
   736         741      0.000216507              1       7.39e-05  
   737         742      0.000213411              1          8e-05  
   738         743      0.000210585              1       0.000125  
   739         744      0.000209504              1       6.43e-05  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   740         745      0.000209366              1       3.05e-05  
   741         746      0.000209361              1       2.99e-05  
   742         747      0.000209361              1       2.99e-05  
   743         748      0.000209361              1       2.99e-05  
   744         749      0.000209361              1       2.99e-05  
   745         750      0.000209361              1       2.99e-05  
   746         751       0.00020936              1       2.99e-05  
   747         752      0.000209359              1       2.99e-05  
   748         753      0.000209356              1       2.98e-05  
   749         754      0.000209349              1       2.98e-05  
   750         755       0.00020933              1       2.97e-05  
   751         756      0.000209282              1       2.97e-05  
   752         757      0.000209172              1       2.97e-05  
   753         758      0.000208961              1       3.01e-05  
   754         759      0.000208685              1       6.31e-05  
   755         760      0.000208507              1       8.64e-05  
   756         761      0.000208463              1       8.62e-05  
   757         762       0.00020846              1        8.2e-05  
   758         763       0.00020846              1       8.11e-05  
   759         764       0.00020846              1       8.09e-05  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   760         765       0.00020846              1       8.03e-05  
   761         766      0.000208459              1       7.94e-05  
   762         767      0.000208459              1        7.8e-05  
   763         768      0.000208457              1       7.57e-05  
   764         769      0.000208452              1       7.19e-05  
   765         770      0.000208439              1       6.57e-05  
   766         771      0.000208406              1       5.58e-05  
   767         772       0.00020832              1       3.96e-05  
   768         773      0.000208098              1       3.04e-05  
   769         774      0.000207545              1       3.01e-05  
   770         775      0.000206262              1       8.18e-05  
   771         776      0.000203773              1       0.000137  
   772         777      0.000200478              1       0.000157  
   773         778      0.000198277              1       0.000111  
   774         779      0.000197724              1       5.55e-05  
   775         780       0.00019768              1       4.05e-05  
   776         781      0.000197679              1          4e-05  
   777         782      0.000197679              1          4e-05  
   778         783      0.000197679              1       3.99e-05  
   779         784      0.000197679              1       3.99e-05  
                                                        First-order 
 Iteration  Func-count       f(x)        Step-size       optimality
   780         785      0.000197679              1       3.98e-05  
   781         786      0.000197679              1       3.97e-05  
   782         787      0.000197678              1       3.95e-05  
   783         788      0.000197675              1       3.92e-05  
   784         789      0.000197669              1       3.94e-05  
   785         790      0.000197652              1       3.96e-05  
   786         791      0.000197611              1        3.8e-05  
   787         792      0.000197512              1       3.39e-05  
   788         793      0.000197309              1       3.09e-05  
   789         794      0.000197008              1        5.6e-05  
   790         795      0.000196764              1       7.18e-05  
   791         796      0.000196684              1       7.93e-05  
   792         797      0.000196676              1       7.74e-05  
   793         798      0.000196676              1       7.59e-05  
   794         799      0.000196676              1       7.57e-05  
   795         800      0.000196676              1       7.52e-05  
   796         801      0.000196676              1       7.44e-05  
   797         802      0.000196675              1       7.32e-05  
   798         803      0.000196674              1       7.12e-05  
   799         804      0.000196671              1       6.79e-05  
                       ...
```

![figure_4.png
](README_images/figure_4.png
)

We transform the truncated initial condition of our test trajectory according to the obtained change of coordinates, and integrate our reduced order evolution rule to predict the development of the trajectory. 

```matlab:Code
[yRec, etaRec, zRec] = advect(IMInfo, RDInfo, yData);
```

# Evaluation of reduced dynamics

The error NMTE is computed as the average distance of the predicted trajectories to the measured ones in the observable space.

```matlab:Code
normedTrajDist = computeTrajectoryErrors(yRec, yData);
NMTE = mean(normedTrajDist)*100
```

```text:Output
NMTE = 
      0.90285

```

We plot the measured test set trajectory in the reduced coordinates and compare it to the prediction. 

```matlab:Code
plotReducedCoordinates(etaData, etaRec)
legend({'Trajectory', 'Prediction'})
```

![figure_5.png
](README_images/figure_5.png
)

We also plot the measured and predicted displacement of the first DIC point.

```matlab:Code
plotTrajectories(yData, yRec, 'm','PlotCoordinate', outdof)
legend({'Trajectory', 'Prediction'})
xlabel('$t \, [$s$]$','interpreter','latex');
ylabel('$u_1 \, [$m$]$','Interpreter','latex')
```

![figure_6.png
](README_images/figure_6.png
)

# Energy as amplitude metric

In our delay-embedding space, we can define numerical differentiation on a state vector. Therefore, we can compute velocities as central differences and compute the kinetic energy of the beam, which is defined as our amplitude metric.

```matlab:Code
nDisplacements = size(uData{1,2},1);
samplingTime = uData{1,1}(2)-uData{1,1}(1);
amplitudeFunction = @(y) mean( ( (y([1:nDisplacements]+2*nDisplacements,:)-y([1:nDisplacements],:))/2/samplingTime).^2)*1.796/2*1e-3; % energy in [mJ]
% Plot SSM with the kinetic energy
plotSSMandTrajectories(IMInfo, amplitudeFunction, yData, zRec, 'NFT', RDInfo.transformation.map)
view(-100,20); legend('off')
zlabel('$K \, [$mJ$]$','Interpreter','latex'); 
```

![figure_7.png
](README_images/figure_7.png
)

# Backbone curves

With the knowledge of the coefficients of the normal form, we extract backbone curves for the instantaneous damping and frequency. The instantaneous damping is strongly nonlinear, but the frequency remains virtually unchanged in our model.

```matlab:Code
zData = transformTrajectories(RDInfo.inverseTransformation.map, etaData);
rhoMax = abs(zData{1,2}(1,1));
BBCInfo = backboneCurves(IMInfo, RDInfo, amplitudeFunction, rhoMax,'Hz');
subplot(121); 
xlabel('$\zeta \, [$\%$]$','Interpreter','latex');
ylabel('$K \, [$mJ$]$','Interpreter','latex'); 
xlim([0.1 0.6]); ylim([0.01 10])
set(gca,'yscale','log')
subplot(122); 
xlabel('$\omega \, [$Hz$]$','Interpreter','latex');
ylabel('$K \, [$mJ$]$','Interpreter','latex'); 
xlim([80 80.4]); ylim([0.01 10])
set(gca,'yscale','log')
```

![figure_8.png
](README_images/figure_8.png
)

# Validation with Accelerations

We now validate our model using the accelerometer data, and the identication results of backbone curves via Peak Finding and Fitting, reported on [2] above.

```matlab:Code
% Select Accelerometer Location
idx_acc = 1;
aData{1,1} = data_BRB.TimeACC;
aData{1,2} = data_BRB.AccelerationACC;
% Linear Interpolation among nodes
tAcc = aData{1,1}; aAcc = aData{1,2}(idx_acc,:); 
loc_acc = data_BRB.LocationACC(idx_acc);
idxs_acc_DIC =sum(Xmesh<loc_acc)+[0 1];
% Define accelerations
accelerationFunctions = @(y) 1e-3*(y([1:nDisplacements]+2*nDisplacements,:)-2*y([1:nDisplacements]+nDisplacements,:)+y([1:nDisplacements],:))/samplingTime.^2; % accelerations in [m/s^2]
accelerationFunction = @(a) (a(idxs_acc_DIC(2),:)-a(idxs_acc_DIC(1),:))*(loc_acc - Xmesh(idxs_acc_DIC(1)))/(Xmesh(idxs_acc_DIC(2)) - Xmesh(idxs_acc_DIC(1))) + a(idxs_acc_DIC(1),:);
tRec = yRec{1,1};
aRec = accelerationFunction( accelerationFunctions (yRec{1,2}));
% Plot results
customFigure('subPlot',[1 2]); subplot(121);
plot(tAcc,aAcc,tRec,aRec)
xlim([tRec(1) tRec(end)])
xlabel('$t \, [$s$]$','interpreter','latex'); 
ylabel('$a \, [$m/s$^2]$','interpreter','latex'); 
xlim([tRec(1) tRec(end)])
subplot(122);
plot(tAcc,aAcc,tRec,aRec)
xlabel('$t \, [$s$]$','interpreter','latex'); 
legend('Measured Acceleration','Reconstructed Acceleration')
ylabel('$a \, [$m/s$^2]$','interpreter','latex'); 
xlim([0.03 0.07])
```

![figure_9.png
](README_images/figure_9.png
)

```matlab:Code
accelerationError = sqrt(mean( (aRec - interp1(tAcc,aAcc,tRec)).^2 ))/max(abs(aRec))*100
```

```text:Output
accelerationError = 
       3.8439

```

```matlab:Code
% Instantaneous damping and frequency in time
instDamping = RDInfo.conjugateDynamics.damping(abs(zRec{1,2}(1,:)));    
instFrequency = RDInfo.conjugateDynamics.frequency(abs(zRec{1,2}(1,:)));  
instDampingRatio = -instDamping./instFrequency*100;   
customFigure('subPlot',[2 1]);
subplot(211);
plot(data_BRB.PFFResultsACC.Time(idx_acc,:),data_BRB.PFFResultsACC.Damping(idx_acc,:)*100,'Linewidth',2,'DisplayName','Peak Finding & Fitting')
plot(tRec,instDampingRatio,'Linewidth',2,'DisplayName',['SSMlearn O(' num2str(ROMOrder) ')'])
xlabel('$t \, [$s$]$','interpreter','latex'); 
ylabel('$\xi \, [$\%$]$','interpreter','latex'); 
legend
xlim([tRec(1) tRec(end)])
subplot(212);
plot(data_BRB.PFFResultsACC.Time(idx_acc,:),data_BRB.PFFResultsACC.Frequency(idx_acc,:),'Linewidth',2,'DisplayName','PFF')
plot(tRec,instFrequency/2/pi,'Linewidth',2,'DisplayName','SSMlearn')
xlabel('$t \, [$s$]$','interpreter','latex'); 
ylabel('$\omega \, [$Hz$]$','interpreter','latex');
xlim([tRec(1) tRec(end)])
```

![figure_10.png
](README_images/figure_10.png
)

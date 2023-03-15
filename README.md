<img src="docs/images/SSMLearnLogo.png" width="350" align="right">

<h1 style="font-family:Helvetica;" align="left">
    SSMLearn
</h1>

### Data-driven Reduced Order Models for Nonlinear Dynamical Systems

This package identifies reduced-order models on spectral submanifolds (SSMs) from data. The required input consists of trajectory data of generic system observables close to an SSM, the SSM dimension, and the polynomial orders of approximation for the parametrization and reduced dynamics.
In addition, an alternative simplified method, fastSSM, is included. See the fastSSM subfolder for documentation.

The computational steps for achieving a reduced-order model are:

1. Embedding of the measurements in a suitable observable space;
2. Computation of the invariant manifold parametrization and its reduced order coordinates;
3. Identification of the reduced dynamics and its normal form.

Once the normal form dynamics has been determined, the code can run analytics and predictions on the reduced-order model, such as backbone curves and forced responses, as shown in our examples.
We have included a demonstration of SSM identification on the following examples.

- Oscillator chain: *n* degrees of freedom with trajectories on or off specific SSMs;
- Von Kármán straight beam in 2D: geometrically nonlinear finite element model from SSMTool, with reduced order models constructed using different observables;
- Brake-Reuss beam: benchmark system for the dynamics of jointed structures, data from experimental measurements (DIC and accelerometers);
- Resonant double beam: structure with a 1:2 internal resonance, data from laser vibrometry;
- Vortex Shedding behind a cylinder: data from CFD simulations, projected on a low-dimensional linear subspace of the phase space;
- Plane Couette flow: reduced order modeling of transitions between exact coherent states;
- Liquid sloshing of a water tank: data from experimental measurements.

This package uses the following external open-source packages for some of the examples and post-processing capabilities:

1. Continuation core (coco) https://sourceforge.net/projects/cocotools/
2. SSMTool 2.3: Computation of invariant manifolds & their reduced dynamics in high-dimensional mechanics problems https://github.com/haller-group/SSMTool-2.3

## Installation
1. Once located in the main folder, install the package:  
    ```sh
    install
    ```
2. If external packages are not yet installed, download SSMTool from the link above, which also include coco, and install it. 
3. (optional) Figure specifications can be edited in the function customFigure.m located in the src folder.
4. You are ready.

## References
Please consider to cite this article when using this code:

- M. Cenedese, J. Axås, B. Bäuerlein, K. Avila and G. Haller. Data-driven modeling and prediction of non-linearizable dynamics via spectral submanifolds. [*Nature Communications*](https://doi.org/10.1038/s41467-022-28518-y), **13** (2022) 872. [[PDF]](https://www.nature.com/articles/s41467-022-28518-y.pdf) [[Supplementary Information]](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-28518-y/MediaObjects/41467_2022_28518_MOESM1_ESM.pdf)

Additional works appear in the references:

- M. Cenedese, J. Axås, H. Yang, M. Eriten and G. Haller. Data-driven nonlinear model reduction to spectral submanifolds in mechanical systems, [*Philosophical Transactions of the Royal Society A: Mathematical, Physical and Engineering Sciences*](https://doi.org/10.1038/s41467-022-28518-y) **380** (2022) 20210194. [[PDF]](http://www.georgehaller.com/reprints/Cenedeseetal_DataDrivenNonlinearModelReduction.pdf) 

- G. Haller, S. Jain and M. Cenedese. Dynamics-based machine learning for nonlinearizable phenomena. Data-driven reduced models on spectral submanifolds, [*SIAM News*](https://sinews.siam.org/Details-Page/dynamics-based-machine-learning-for-nonlinearizable-phenomena) **55/5** (2022) 20210194. [[PDF]](http://www.georgehaller.com/reprints/HallerJainCenedese_dynamics_based_machine_learning.pdf) 

- B. Kaszás, M. Cenedese & G. Haller, Dynamics-based machine learning of transitions in Couette flow, [*Physical Review Fluids*](https://doi.org/10.1103/PhysRevFluids.7.L082402) **7** (2022) L082402. [[PDF]](http://www.georgehaller.com/reprints/dynamicsbasedmachinelearning.pdf) [[Supplemental Material]](https://journals.aps.org/prfluids/supplemental/10.1103/PhysRevFluids.7.L082402/supplemental_couette.pdf)

- J. Axås, M. Cenedese & G. Haller, Fast data-driven model reduction for nonlinear dynamical systems, [*Nonlinear Dynamics*](https://doi.org/10.1007/s11071-022-08014-0) (2022). [[PDF]](https://link.springer.com/content/pdf/10.1007/s11071-022-08014-0.pdf)

Please report any issues/bugs to Mattia Cenedese (mattiac@ethz.ch) or Joar Axås (jgoeransson@ethz.ch).

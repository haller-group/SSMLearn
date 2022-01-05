<img src="docs/images/SSMLearnLogoBlack.pdf" width="290" align="right">

<h1 style="font-family:Helvetica;" align="left">
    SSMLearn
</h1>

### Data-driven Reduced Order Models for Nonlinear Dynamical Systems

This package perform data-driven identification of reduced order model based on spectral submanifolds (SSMs). The required input consists of trajectory data for generic system observables, the SSM dimension and the polynomial orders for approximation of the parametrization and reduced dynamics. The code implicitly assumes that the data lays close to an autonomous SSM of the system.

The computational steps for achieving a reduced-order model are:

1. Embedding of the measurements in a suitable observable space
2. Computation of the invariant manifold parametrization and its reduced order coordinates
3. Identification of the reduced dynamics and its normal form

Once the normal form dynamics has been determined, the code can run analytics and predictions on the reduced-order model, such as backbone curves and forced responses, as shown in our examples. There are no constraints on the types of masurements, on the kind of nonlinearities or on the problem dimensions.

We have included a demonstration of SSM computation over the following examples.

- Stuart-Landau: trivial example for a 1D dynamical system
- Oscillator chain: n degrees of freedom with trajectories on or off specific SSMs
- Von Kármán straight beam in 2D: model from SSMTool, with reduced order models constructed using different observables
- Vortex Shedding behind a cylinder: data from CFD simulations, projected on a low-dimensional linear subspace of the phase space

This package uses the following external open-source packages for some of the examples and post-processing capabilities:

1. Continuation core (coco) https://sourceforge.net/projects/cocotools/
2. SSMTool 2.1: Computation of invariant manifolds & their reduced dynamics in high-dimensional mechanics problems https://github.com/haller-group/SSMTool-2.1

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

- M. Cenedese, J. Axås, B. Bäuerlein, K. Avila and G. Haller. Data-driven modeling and prediction of non-linearizable dynamics via spectral submanifolds. *Submitted* (2021).

Additional works appear in the references:

- M. Cenedese, J. Axås, H. Yang, M. Eriten and G. Haller. [Data-driven nonlinear model reduction to spectral submanifolds in mechanical systems](https://arxiv.org/pdf/2110.01929.pdf),  *arXiv:2110.01929* (2021). 

Please report any issues/bugs to Mattia Cenedese (mattiac@ethz.ch) or Joar Axås (jgoeransson@ethz.ch)

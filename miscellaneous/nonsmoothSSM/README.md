# Nonsmooth SSMs

Repository for sharing equation- and data-driven reduced order modeling of non smooth mechanical systems based on primary Spectral Submanifolds (SSMs). We discuss simple lumped systems and a finite element discretisation of a nonlinear beam, featuring different types of friction and contact nonlinearities. 

This repository is temporary, as the documentation is subject to potential modifications. 

This repository uses the following external open-source packages for some of the examples and post-processing capabilities:

1. `SSMLearn`: learning of invariant manifolds & their reduced dynamics from data of high-dimensional dynamical systems https://github.com/haller-group/SSMLearn
2. `SSMTool`: computation of invariant manifolds & their reduced dynamics in high-dimensional dynamical systems https://github.com/haller-group/SSMTool-2.4

## Installation 
1. Once located in the main folder, install the package:
	```sh
	install
	```
2. If external packages are not installed yet, download and install them.
3. (optional) Figure specifications can be edited in the function customFigure.m located in the src folder of `SSMLearn`.
4. You are ready.

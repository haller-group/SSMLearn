fastSSM and fastSSMplus locate spectral submanifolds (SSMs) from data and 
compute their reduced dynamics in the normal form. [1]

To add fastSSM(plus) to the matlab path, run the following matlab commands:

    cd SSMLearn/
    install

fastSSM runs without any dependencies. To run fastSSMplus, SSMTool is needed.
The latest version of SSMTool is available at https://github.com/haller-group
At the time of writing, the following commands work:

    git clone https://github.com/haller-group/SSMTool-2.1.git
    cd SSMTool/
    install

You can then run the following three examples: 
    
    examples/
        sloshing/sloshing.m
        vonkarmanbeam/vonkarman.m
        resonantbeam/resonantdoublebeam.m

fastSSM and fastSSMplus are open-source and free to use.
Please consider citing this article when using the code:

[1] J. Axås, M. Cenedese, and G. Haller. Fast data-driven model reduction for nonlinear dynamical systems. Nonlinear Dynamics, 2022.

Nov 1st, 2022
Joar Axås
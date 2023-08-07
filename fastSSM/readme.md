fastSSM and fastSSMplus locate spectral submanifolds (SSMs) from data and 
compute their reduced dynamics in the normal form. [1]

To add fastSSM(plus) to the matlab path, run the following matlab commands:

    cd SSMLearn/
    install

fastSSM runs without any dependencies. To run fastSSMplus, SSMTool 2.2 is needed.
The latest version of SSMTool is available at github.com/jain-shobhit/SSMTool
At the time of writing, the following commands work:

    git clone https://github.com/jain-shobhit/SSMTool.git
    cd SSMTool/
    install

You can then run the following examples: 
    
    examples/
        sloshing/sloshing.m
        vonkarmanbeam/vonkarman.m
        resonantbeam/resonantdoublebeam.m
        
 These examples cover how to employ Vandermonde matrices for delay-embedded data: [2]
    
    examples/
        oscillator/oscillator.m
        sloshing/sloshingDelay.m
        vonkarmanbeam/vonkarmanDelay.m

fastSSM and fastSSMplus are open-source and free to use.
Please consider citing this work when using the code:

[1] J. Axås, M. Cenedese, and G. Haller. Fast data-driven model reduction for nonlinear dynamical systems. Nonlinear Dynamics, 2022.

[2] J. Axås and G. Haller.Model reduction for nonlinearizable dynamics via delay-embedded spectral submanifolds. Nonlinear Dynamics, 2023.

Aug 7th, 2023
Joar Axås
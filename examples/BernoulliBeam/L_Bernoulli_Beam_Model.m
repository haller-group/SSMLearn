function [M,C,K]=L_Bernoulli_Beam_Model(n)
%This code gives the Mass, Damping and Stiffness Matrix for a Forced linear
%Bernoulli Beam with nonlinear cubic spring attached to the last node. It
%also gives a forcing amplitude vector f. The beam is forced in terms of an
%displacement imposed on the last node. The equations of motion are
%M*u''+C*u'+K*u+f_nl=f. Structural damping with coefficients alpha and beta
%is used for the damping matrix C.

L=1;                                 %Length of beam [mm]
h=20e-3;                                   %Height of beam [mm]
b=50e-3;                                   %Width of beam [mm]
Lel=L/n;                                %Auxiliary variable needed for integration
E=70e9;                             %Young's modulus [kPa]   
rho=2700;                       %Density [kg/mm^3]    
m0=b*h*rho;                             %Mass inertia
I=h^3*b/12;                             %Area moment of inertia
co = 0.125;
alpha=0.003*co;                         %Structural Damping Parameters
alpha=0;
beta=0.006*co;
beta=0.4286e-4;


M_el=m0*Lel/420*[156, -22*Lel, 54, 13*Lel;-22*Lel, 4*Lel^2, -13*Lel, -3*Lel^2; 54, -13*Lel, 156, 22*Lel; 13*Lel, -3*Lel^2, 22*Lel, 4*Lel^2];
K_el=2*E*I/Lel^3*[6, -3*Lel, -6, -3*Lel;-3*Lel, 2*Lel^2, 3*Lel, Lel^2; -6, 3*Lel, 6, 3*Lel; -3*Lel, Lel^2, 3*Lel, 2*Lel^2];


if n==1
    M=M_el;
    K=K_el;
else
    n_el=2*(1+n);
    M=zeros(n_el-2);
    K=M;
    M(1:2,1:2)=M_el(3:4,3:4);
    K(1:2,1:2)=K_el(3:4,3:4);
    for i=1:n-1
        K(2*(i-1)+1:2*(i-1)+4,2*(i-1)+1:2*(i-1)+4)=K(2*(i-1)+1:2*(i-1)+4,2*(i-1)+1:2*(i-1)+4)+K_el;
        M(2*(i-1)+1:2*(i-1)+4,2*(i-1)+1:2*(i-1)+4)=M(2*(i-1)+1:2*(i-1)+4,2*(i-1)+1:2*(i-1)+4)+M_el;
    end
end
        

C=alpha*M+beta*K;

end
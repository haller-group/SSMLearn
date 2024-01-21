function [alpha, omega, nu] = cal_parameters(N,l)

syms x
% N = 2; % number of modes
% l = 2;

% construct modal functions
phis   = [];
dphis  = [];
ddphis = [];
alphal = zeros(N,1);
alphal(1:4) = [3.927 7.069 10.210 13.352]';
alphal(5:N) = (4*(5:N)'+1)*pi/4;
alpha  = alphal/l;
R = sin(alphal)./sinh(alphal);
E = (0.5*l*(1-R.^2)+(R.^2.*sinh(2*alphal)-sin(2*alphal))/(4*alpha)).^(-0.5);
for n=1:N
   phin   = E(n)*(sin(alpha(n)*x)-R(n)*sinh(alpha(n)*x));
   dphin  = diff(phin,x);
   ddphin = diff(dphin,x);
   phis   = [phis; phin];
   dphis  = [dphis; dphin];
   ddphis = [ddphis; ddphin];
   normphi = int(phin^2, x, 0, l);
   fprintf('L2 norm of phi%d is %d\n', n, normphi);
end

% compute nonlinear coeffficients
alpha1 = zeros(N);
alpha2 = zeros(N);
alpha  = zeros(N,N,N,N);

for i=1:N
    for j=1:N
        alpha1(i,j) = int(phis(i)*ddphis(j),x,0,l);
        alpha2(i,j) = int(dphis(i)*dphis(j),x,0,l);
    end
end
        

for n=1:N
    for m=1:N
        for p=1:N
            for q=1:N
                alpha(n,m,p,q) = alpha1(n,q)*alpha2(m,p);
            end
        end
    end
end

% compute natural frequencies
omega = zeros(N,1);
for i=1:N
    int1 = int(phis(i)^2,x,0,l);
    int2 = int(diff(phis(i),x,4)*phis(i),x,0,l);
    omega(i) = sqrt(int2/int1);
end

nu = 1/(2*l);

end
    
                
                
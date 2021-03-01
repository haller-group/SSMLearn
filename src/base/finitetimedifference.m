function [dXdt,X,t] = finitetimedifference(X_in,t_in,halfw)
% Central finite time difference with uniform monodimensional grid spacing
% of customized accuracy, equal to 2*halfw. X is a matrix of n_variables 
% x n_instances. Thefore, the finite difference is implemented columnwise.

if halfw>4
    disp('The cofficients for this accuracy are not present in the current implementation. The finite difference is computed with accuracy O(Dt^8)')
    halfw = 4;
end
% Coefficients for the numerical derivative
coeff_mat = [1/2   2/3   3/4    4/5; ...
               0 -1/12 -3/20   -1/5; ...
               0     0  1/60  4/105; ...
               0     0     0 -1/280];
% Computation
base_int = halfw+1:size(X_in,2)-halfw;
X = X_in(:,base_int); t = t_in(base_int); dt = t(2)-t(1);
dX = zeros(size(X));
for ii = 1:halfw
    dX = dX + coeff_mat(ii,halfw) * ...
                 (X_in(:,base_int+ii) - X_in(:,base_int-ii));
end
dXdt = dX/dt;
end

function [omega,lambda,t_eig] = sigprocess_DMD(X,t,window_length,window_overlap,r)
% Computes the Dynamic Mode Decomposition for each
% time window with overlap
%
% INPUTS:
% X: data matrix, time indexes along columns
% window_length: samples per window
% window_overlap: real number in [0,1)
% r: target rank of SVD
% dt: time step
% OUTPUTS :
% omega: the continuous -time DMD eigenvalues
% lambda: the discrete -time DMD eigenvalues
% t_eig: times at which these eigenvalues are known

dt = t(2)-t(1);
lambda = []; omega = []; t_eig = [];
ind_end = window_length + 1;
stepping = round(window_length*(1-window_overlap));
while ind_end < size(X,2)
    ind_i  = ind_end-window_length:ind_end;
    X_i = X(:,ind_i);  Y_i = X(:,ind_i+1);
    [~,omega_i,lambda_i] = DMD(X_i,Y_i,r,dt);
    [~,pos] = sort(abs(omega_i),'ascend');
    lambda = [lambda lambda_i(pos)];
    omega = [omega omega_i(pos)];
    t_eig = [t_eig (t(ind_i(1))+t(ind_i(end)))/2];
    ind_end = ind_end + stepping ;
end 
end
               

function [V,s2] = PCA(X,k)
% Principal Compoenent Analysis of data via SVD. This is also known as the
% snapshot Principal Orthogonal Decomposition
%
% REQUIRED INPUT
%    X    - matrix of dimension n x N, where n is the number of features
%           and N that of the number of data-points or a cell array of
%           dimension (N_traj,2) where the first column contains time
%           instances (1 x mi each) and the second column the trajectories
%           (n x mi each)
%    k    - number of modes
%
% OUTPUT
%   V     - matrix of dimension n x k with the modes
%  s2     - square singular values

if iscell(X)==1
    X_cell = X; X = []; t = [];
    for ii = 1:size(X_cell,1)
        X = [X [X_cell{ii,2}]]; t = [t [X_cell{ii,1}]];
    end 
end
[~,S,V]=svds(transpose(X),k);
s2 = diag(S);

end
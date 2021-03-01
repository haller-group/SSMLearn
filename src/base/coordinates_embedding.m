function [x_embd,opts_embd] = coordinates_embedding(x,k,varargin)
% Returns the n-dim. time series x into a time series of properly embedded
% coordinate system y of dimension p. Optional inputs to be specified as
% 'field_name','field value'
%
% INPUT
% x - cell array of dimension (N_traj,2) where the first column contains
%     time instances (1 x mi each) and the second column the trajectories
%     (n x mi each)
% k - dimension of the invariant manifold to learn
%
% OPTIONAL INPUT
% OverEmbedding  - augment the minimal embedding dimension with a number of
%                  time delayed measurements, default 0
% ForceEmbedding - force the embedding in the states of x, default false
% TimeStepping   - time stepping in the time series, default 1
%
% If varargin is set to an integer value, it it set as OverEmbedding
%
% OUTPUT
% y - cell array of dimension (N_traj,2) where the first column contains
%     time instances (1 x mi each) and the second column the trajectories
%     (p x mi each)
% opts_embd - options containing the embedding information
%
% Developed by Mattia Cenedese. Updated March 2021.

opts_embd = struct('IMdimensions',k,'OverEmbedding',0,...
    'ForceEmbedding',0,'TimeStepping',1);
% Default case
if nargin == 3; opts_embd.OverEmbedding = varargin{:}; end
if rem(length(varargin),2) > 0 && length(varargin) > 1
    error('Error on input arguments. Missing or extra arguments.')
end
% Custom options
if nargin > 3
    for ii = 1:length(varargin)/2
        opts_embd = setfield(opts_embd,varargin{2*ii-1},...
            varargin{2*ii});
    end
end

% Determine number of embedding coordinate system
N_traj = size(x,1);
N = size([x{1,2}],1);
if length([x{1,1}])==N; N = size([x{1,2}],2); end
l = opts_embd.TimeStepping;
n_N = ( ceil( (2*k+1)/N ) + opts_embd.OverEmbedding);
% Construct embedding coordinate system
x_embd = cell(N_traj,2);
ind_traj = cell(N_traj,1); idx_end = 0;
if n_N > 1 && opts_embd.ForceEmbedding ~= 1
    p = n_N * N;
    % Augment embedding dimension with time-delays
    if N == 1
        disp(['The ' num2str(p) ' embedding coordinates consist of the '...
            'measured state and its ' num2str(n_N-1) ...
            ' time-delayed measurements.'])
    else
        disp(['The ' num2str(p) ' embedding coordinates consist of the '...
            num2str(N) ' measured states and their ' ...
            num2str(n_N-1) ' time-delayed measurements.'])
    end
    for jj = 1:N_traj
        
        t_j = [x{jj,1}]; x_j = [x{jj,2}];
        % Check dimensions
        if size(t_j,1)>size(t_j,2); t_j = transpose(t_j); end
        if length(t_j)~=size(x_j,2); x_j = transpose(x_j); end
        t_j = t_j(1:l:end); x_j = x_j(:,1:l:end);
        Y_j = x_j(:,1:end-(n_N-1));
        for ii = 1:(n_N-1)
            Y_j = [Y_j; x_j(:,1+ii:end-(n_N-1-ii))];
        end
        x_embd{jj,1} = t_j(1+ii:end-(n_N-1-ii)); x_embd{jj,2} = Y_j; 
    end
    
else
    p = N;
    disp(['The embedding coordinates consist of the measured states.'])
    for jj = 1:N_traj
        t_j = [x{jj,1}]; x_j = [x{jj,2}];
        % Check dimensions
        if size(t_j,1)>size(t_j,2); t_j = transpose(t_j); end
        if length(t_j)~=size(x_j,2); x_j = transpose(x_j); end
        t_j = t_j(1:l:end); x_j = x_j(:,1:l:end);
        x_embd{jj,1} = t_j; x_embd{jj,2} = x_j;
    end
end
% Generate a string to describe the states
if n_N > 1 && opts_embd.ForceEmbedding ~= 1
    strStateSpace = 'x(t), x(t+Dt)';
    if n_N > 2
        for ii = 2:(n_N-1)
            strStateSpace = append(strStateSpace,[', x(t+' ...
                                                      num2str(ii) 'Dt)']);
        end
    end
else
    strStateSpace = 'x(t)';
end
opts_embd.Observables = p;
opts_embd.ObservableStateSpace = strStateSpace;
end

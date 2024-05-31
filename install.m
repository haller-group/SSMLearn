function install
% initial notes
    disp('SSMLearn: Data-driven Reduced Order Models for Nonlinear Dynamical Systems')
    disp('License: GNUv3.0')
    
% Install the main software 
addpath(genpath('src'))
addpath('fastSSM')
disp('SSMLearn successfully installed!')
disp(' ')

% savepath % Uncomment if you would like to add SSMLearn in the Matlab default path. Run the command restoredefaultpath to restore the default path 

% check for SSMTool installation
disp('Checking SSMTool is installed...')
if exist('misc/frc_ab', 'file')
    disp('SSMTool v>=2.4 already installed!')
else
    warning(['SSMTool 2.4 not installed, cf. readme file. ' ...
        'Clone from github.com/jain-shobhit/SSMTool'])
end

end



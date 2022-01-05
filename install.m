function install
% initial notes
    disp('SSMLearn: Data-driven Reduced Order Models for Nonlinear Dynamical Systems')
    disp('Maintained by Mattia Cenedese (mattiac@ethz.ch) and Joar Axås (jgoeransson@ethz.ch)')
    disp('License: GNUv3.0')
    
% Install the main software 
addpath(genpath('src'))

% disp('Installing external packages ...')
% prompt = 'Do you want to install SSMTool and coco? Y/N [Y]: ';
% str = input(prompt,'s');
% if isempty(str)
%     str = 'Y';
% end
% if strcmp(str,'Y') == 1
%     % Install all external packages 
%     run([pwd '/ext/SSMTool/install.m']);
% end
% disp('Done.')

savepath
end



function [prob, args1, args2] = monitor_states(prob, ispolar, m)
% MONITOR_STATES This function add state of reduced dynamics as
% continuation parameters and define string array for the state. The prob
% here is a continuation problem. (args1,args2)=(rhoargs,thargs) or 
% (Reargs,Imargs) depending on the value of ispolar

args1 = cell(m,1);
args2 = cell(m,1);
if ispolar
    for k=1:m
        args1{k} = strcat('rho',num2str(k));
        args2{k}  = strcat('th',num2str(k));
    end
    prob = coco_add_pars(prob, 'radius', 1:2:2*m-1, args1(:)');
    prob = coco_add_pars(prob, 'angle', 2:2:2*m, args2(:)');
else
    for k=1:m
        args1{k} = strcat('Rez',num2str(k));
        args2{k} = strcat('Imz',num2str(k));
    end
    prob = coco_add_pars(prob, 'realParts', 1:2:2*m-1, args1(:)');
    prob = coco_add_pars(prob, 'imagParts', 2:2:2*m, args2(:)');
end 

end

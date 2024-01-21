function [mass,damp,stiff,fnl,fext] = buildModel(c,f,epsilon,n)

l = 2;
disp('Getting nonlinearity coefficients')
fileName = ['tensors_',num2str(n),'_',num2str(l),'.mat'];
try 
    load(fileName,'alpha','omega','nu');   
    disp('Loaded coefficients from storage');
catch
    disp('Calculating coefficients');
    [alpha, omega, nu] = cal_parameters(n,l);
    disp('Saving coefficients');
    save(fileName,'alpha','omega','nu','-v7.3')
end

mass = eye(n);
damp = 2*c*eye(n)*epsilon;
stiff = diag(omega.^2);

f3 = -epsilon*nu*sptensor(alpha);
fnl = {[],f3};
fext = f*epsilon;

end
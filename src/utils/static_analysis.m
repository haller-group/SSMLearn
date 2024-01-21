function [uNonlinear,uLinear,uNonlinearOut,uLinearOut, relativeDispDifference] = static_analysis(Model, loadVector,nsteps,outdof)
%STATIC_ANALYSIS computes nonlinear static equilibrium by incremental
%loading in predefined steps. 
loadCoefficients = (1:nsteps)/nsteps;
uInit = zeros(Model.Mesh.nDOFs,1);
[u_lin, U] = static_equilibrium(Model, uInit, Model.unconstrain_vector(loadVector), 'nsteps',nsteps,'method', 'newton');
grid on;
uNonlinear = Model.constrain_vector(U);
uLinear = Model.constrain_vector(u_lin)*loadCoefficients;

%% Plotting
uNonlinearOut = uNonlinear(outdof,:);
uLinearOut = uLinear(outdof,:);
% plot results
figure; hold on; grid on;
plot(loadCoefficients,uLinearOut,'r--','DisplayName','Linear')
plot(loadCoefficients,uNonlinearOut,'-k','DisplayName','Nonlinear')
legend('location','best')
ylabel('$$u_{\mathrm{out}}$$ [m]'); 
xlabel('Normalized load'); 
title('Static loading analysis');
axis tight

figure; hold on; grid on;
relativeDispDifference = 100*abs(uNonlinearOut-uLinearOut)./abs(uNonlinearOut);
plot(loadCoefficients,relativeDispDifference,'Linewidth',2)
axis tight
ylabel('$$\frac{|u_{\mathrm{out}} - u_{\mathrm{out,linear}}|}{|u_{\mathrm{out}}|}$$ [\%]'); 
xlabel('Normalized load'); 
end


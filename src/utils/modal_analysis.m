function [phi,relativeDiffForceNorm, relativeDiffForceOut] = modal_analysis(Model,scalingFactor,nsteps,outdof,varargin)

if nargin > 5
    FE = varargin{1};
    if nargin == 6
        iMode = varargin{2};
    else 
        iMode = 1;
    end
else
    FE = true;
    iMode =1;
end

scalingFactors = scalingFactor*(1:nsteps)/nsteps;
if FE
    M = Model.DATA.M;
    K = Model.DATA.K;
    
    M = Model.constrain_matrix(M);
    K = Model.constrain_matrix(K);
else
    M = Model.M;
    K = Model.K;
end
%% Eigenvalue problem
disp('Solving undamped eigenvalue problem')
n_VMs = 5; % first n_VMs modes with lowest frequency calculated
[V0,D] = eigs(K,M,n_VMs,'SM');
phi = V0(:,iMode);
disp('Mode shape')
disp(phi)
disp('Eigenfrequency')
disp(sqrt(D(iMode, iMode)))

n = size(M,1);
%% Internal force computation for different levels of modal displacement
intForceNormLinear = zeros(size(scalingFactors));
intForceNormNonlinear =  zeros(size(scalingFactors));
relativeDiffForceNorm = zeros(size(scalingFactors));
intForceNonlinearOut = zeros(size(scalingFactors));
intForceLinearOut = zeros(size(scalingFactors));
relativeDiffForceOut = zeros(size(scalingFactors));

for j = 1:length(scalingFactors)
    phiScaled = scalingFactors(j)*phi;
    intForceLinear = K*phiScaled;
    if FE
        u = Model.unconstrain_vector(phiScaled);
        [~,Fint] = Model.tangent_stiffness_and_force(u);
        intForceNonlinear = Model.constrain_vector(Fint);
    else
         rhs = Model.F(0,[phiScaled; zeros(n,1)]);
         intForceNonlinear = -M*rhs(n+1:2*n,:);
    end
    
    intForceNonlinearOut(j) = intForceNonlinear(outdof);
    intForceLinearOut(j) = intForceLinear(outdof);
    intForceNormLinear(j) = norm(intForceLinear);
    intForceNormNonlinear(j) = norm(intForceNonlinear);
    relativeDiffForceNorm(j) = 100*norm(intForceLinear-intForceNonlinear)/norm(intForceNonlinear);
    relativeDiffForceOut(j) = 100*abs(intForceLinearOut(j)-intForceNonlinearOut(j))/abs(intForceNonlinearOut(j));
end

%% Plotting
customFigure('subPlot', [1 2]);
subplot(121)
plot(scalingFactors,intForceNormLinear,'r--','DisplayName','Linear ($$ \| s\mathbf{K\phi}\|$$)','Linewidth',2)
plot(scalingFactors,intForceNormNonlinear,'-k','DisplayName','Nonlinear ($$ \| s\mathbf{K\phi} + \mathbf{f}(s\mathbf{\phi})\|$$)','Linewidth',2)
%legend('location','best','interpreter','latex')
ylabel('$$\|\mathbf{f}_{\mathrm{int}}\|$$','interpreter','latex');
xlabel('Modal scaling');
axis tight
subplot(122)
plot(scalingFactors*phi(outdof),intForceNormLinear,'r--','DisplayName','Linear ($$ \| s\mathbf{K\phi}\|$$)','Linewidth',2)
plot(scalingFactors*phi(outdof),intForceNormNonlinear,'-k','DisplayName','Nonlinear ($$ \| s\mathbf{K\phi} + \mathbf{f}(s\mathbf{\phi})\|$$)','Linewidth',2)
legend('location','best','interpreter','latex','interpreter','latex');
ylabel('$$\|\mathbf{f}_{\mathrm{int}}\|$$','interpreter','latex');
xlabel('$$q_{\mathrm{out}}$$','interpreter','latex');
axis tight

customFigure('subPlot', [1 2]);
subplot(121)
plot(scalingFactors,relativeDiffForceNorm,'Linewidth',2)
axis tight
ylabel('$$\frac{ \| \mathbf{f}(s\mathbf{\phi})\|}{ \| s\mathbf{K\phi} + \mathbf{f}(s\mathbf{\phi})\|}$$ [\%]','interpreter','latex');
xlabel('Modal scaling');
subplot(122)
plot(scalingFactors*phi(outdof),relativeDiffForceNorm,'Linewidth',2)
axis tight
ylabel('$$\frac{ \| \mathbf{f}(s\mathbf{\phi})\|}{ \| s\mathbf{K\phi} + \mathbf{f}(s\mathbf{\phi})\|}$$ [\%]','interpreter','latex');
xlabel('$q_{\mathrm{out}}$','interpreter','latex','interpreter','latex');

customFigure('subPlot', [1 2]);
subplot(121)
plot(scalingFactors,intForceLinearOut,'r--','DisplayName','Linear','Linewidth',2)
plot(scalingFactors,intForceNonlinearOut,'-k','DisplayName','Nonlinear','Linewidth',2)
%legend('location','best','interpreter','latex')
ylabel('$f_{\mathrm{int, out}}$','interpreter','latex');
xlabel('Modal scaling');
axis tight
subplot(122)
plot(scalingFactors,intForceLinearOut,'r--','DisplayName','Linear','Linewidth',2)
plot(scalingFactors,intForceNonlinearOut,'-k','DisplayName','Nonlinear','Linewidth',2)
legend('location','best')
ylabel('$f_{\mathrm{int, out}}$','interpreter','latex');
xlabel('$q_{\mathrm{out}}$','interpreter','latex');
axis tight

customFigure('subPlot', [1 2]);
subplot(121)
plot(scalingFactors,relativeDiffForceOut,'Linewidth',2)
axis tight
ylabel('Relative nodal force difference [%]');
xlabel('Modal scaling');
subplot(122)
plot(scalingFactors,relativeDiffForceOut,'Linewidth',2)
axis tight
ylabel('Relative nodal force difference [%]');
xlabel('$q_{\mathrm{out}}$','interpreter','latex');

disp(['Displacement at output DOF: ' num2str(scalingFactors(end)*phi(outdof))])
end


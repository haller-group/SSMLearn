function [Mmap, iMmap, Tmap, iTmap, Nflow, yRec, BBCInfo] = fastSSM(yData, mfdorder)
% [Mmap, iMmap, Tmap, iTmap, Nflow, yRec, BBCInfo] = fastSSM(yData, mfdorder)
% Computes a 2D SSM and O(3) normal form from data. By Joar Axas (jgoeransson@ethz.ch)
%
% INPUT:                                          OUTPUT:
% yData     {Nx2}  cell array with one row per    Mmap    function  SSM parametrization y=M(xi)^{1:m}
%                 input trajectory, where the     iMmap   function  SSM chart
%                 first column contains time      Tmap    function  normal form transformation xi=T(z)
%                 values (1xn each) and the       iTmap   function  inverse normal form transformation
%                 second column the trajectories  Nflow   function  normal form dynamics zdot=N(z)
%                 (pxn each)                      yRec    (pxn)     reconstruction of yData{1,2}
% mfdorder  int    SSM polynomial order           BBCInfo struct    computed backbone curves
%
% EXAMPLE (data from linear system):
% t = 0:0.1:10
% yData = {t, [sin(t);cos(t)]}
% [Mmap, iMmap, Tmap, iTmap, Nflow, yRec, BBCInfo] = fastSSM(yData, 1)

%% Input, settings
romorder = 3; nforder = 3; mfddim = 2;
t = horzcat(yData{:,1}); Y = horzcat(yData{:,2});
iStart = 1;
for iTraj = 1:size(yData,1); iStart(iTraj+1) = iStart(iTraj)+size(yData{iTraj,1},2); end

%% Fit manifold
[u,s,v] = svds(Y, mfddim);
V = (s\u'./max(abs(v'),[],2))'; % tangent space
iMmap = @(y) V'*y;
Xi = iMmap(Y);
M = Y/phi(Xi, 1:mfdorder);
Mmap = @(xi) M*phi(xi, 1:mfdorder);

%% Compute time derivative and reduced dynamics
Xidot = []; Xinew = [];
for iTraj = 1:size(yData,1) % finite time difference
    [Xdot,Xnew] = ftd(Xi(:,iStart(iTraj):iStart(iTraj+1)-1),t);
    Xidot = [Xidot, Xdot]; Xinew = [Xinew, Xnew];
end
R = Xidot/phi(Xinew, 1:romorder); % reduced dynamics
[W,Lambda] = eig(R(1:2,1:2));
G = W\Xidot/phi(W\Xinew,1:romorder); % diagonalize R

%% Compute transformation to normal form

ll = Lambda(1,1); % lambda_l
lc = conj(ll);    % lambda_l+1
lj = [ll; lc];    % lambda_j
T1 = eye(2);
T2 = [G(:,3)./(2*ll-lj), G(:,4)./(ll+lc-lj), G(:,5)./(2*lc-lj)];
T3 = horzcat(([2,1].*G(:,3:4)*T2(:,1) + G(:,6))./(3*ll-lj), ...
      [0;1].*([2,1].*G(:,3:4)*T2(:,2) + [1,2].*G(:,4:5)*T2(:,1) + G(:,7))./(2*ll+lc-lj), ...
      [1;0].*([2,1].*G(:,3:4)*T2(:,3) + [1,2].*G(:,4:5)*T2(:,2) + G(:,8))./(ll+2*lc-lj), ...
             ([1,2].*G(:,4:5)*T2(:,3) + G(:,9))./(3*lc-lj));
T = [T1, T2, T3];
Tmap = @(z) real(W*T*phi([z;conj(z)], 1:romorder));

gamma = [2,1].*G(1,3:4)*T2(:,2) + [1,2].*G(1,4:5)*T2(:,1) + G(1,7);
Nflow = @(t,z) ll*z + gamma*z.^2.*conj(z);

%% Approximate inverse to the normal form transformation
iT1 = [1,0];
iT2 = -T2(1,:);
iT3 = -T3(1,:) + horzcat(2*T(1,3).^2 + T(1,4)*conj(T(1,5)), ...
    3*T(1,3)*T(1,4) + T(1,4)*conj(T(1,4)) + 2*T(1,5)*conj(T(1,5)), ...
    2*T(1,3)*T(1,5) + T(1,4)*conj(T(1,3)) + T(1,4)*T(1,4) + 2*T(1,5)*conj(T(1,4)), ...
      T(1,4)*T(1,5) + 2*T(1,5)*conj(T(1,3)));
iT = [iT1, iT2, iT3];
iTmap = @(xi) iT*phi(W\xi, 1:romorder);

%% Evaluate model
fprintf('\\dot{\\rho} = %6.3f\\rho %+6.3f\\rho^3\n', real(ll), real(gamma))
fprintf('\\dot{\\theta} = %6.3f %+6.3f\\rho^2\n', imag(ll), imag(gamma))

[~, zRec] = ode45(Nflow, t(iStart(1):iStart(2)-1), iTmap(Xi(:,1)), odeset('RelTol', 1e-6));
yRec = Mmap(Tmap(zRec.'));

%% Backbone curves
yObservable = @(y) max(abs(y(1,:)));
rho = linspace(0,1);
theta = linspace(0,2*pi);
for ir = 1:length(rho)
    zBBC = rho(ir).*exp(1i*theta);
    uBBC(ir) = yObservable(Mmap(Tmap(zBBC)));
end
damp = @(rho) real([ll, gamma])*phi(rho,0:2:romorder);
freq = @(rho) imag([ll, gamma])*phi(rho,0:2:romorder);
figure
subplot(1,2,1)
plot(damp(rho), uBBC, 'LineWidth', 2)
set(gca, 'FontSize', 18)
set(gca, 'FontName', 'times')
grid on
xlabel('Damping [1/s]', 'interpreter', 'latex'); ylabel('Amplitude', 'interpreter', 'latex')
subplot(1,2,2)
plot(freq(rho), uBBC, 'LineWidth', 2)
set(gca, 'FontSize', 18)
set(gca, 'FontName', 'times')
grid on
xlabel('Frequency [rad/s]', 'interpreter', 'latex')
BBCInfo = struct('damp', damp, 'freq', freq, 'damping', damp(rho), 'frequency', freq(rho), 'amplitude', uBBC, 'amplitudeNormalForm', rho);
end

%% Multivariate polynomial
function exps = exponents(d, k)
    B = repmat({0:max(k)}, 1, d);
    A = combvec(B{:}).';
    exps = A(ismember(sum(A, 2), k),:);
    [~,ind] = sort(sum(exps, 2));
    exps = exps(ind,:);
end

function u = phi(xi, r) % return monomials
    x = reshape(xi, 1, size(xi, 1), []);
    exps = exponents(size(xi, 1), r);
    u = reshape(prod(x.^exps, 2), size(exps, 1), []);
end

function [dXdt, Xtrunc] = ftd(X, t) % finite time difference
    ind = 5:size(X,2)-4; Xtrunc = X(:,ind);
    dX = 4/5*(X(:,ind+1)-X(:,ind-1)) - 1/5*(X(:,ind+2)-X(:,ind-2)) + ...
        4/105*(X(:,ind+3)-X(:,ind-3)) - 1/280*(X(:,ind+4)-X(:,ind-4));
    dXdt = dX./(t(ind+1)-t(ind));
end
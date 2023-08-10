function [Mmap, iMmap, Tmap, iTmap, Nflow, yRec, BBCInfo] = fastSSMplus(yData, mfddim, mfdorder, romorder, nforder, varargin)
% [Mmap, iMmap, Tmap, iTmap, Nflow, yRec, BBCInfo] = fastSSMplus(yData, mfddim, mfdorder, romorder, nforder)
% [Mmap, iMmap, Tmap, iTmap, Nflow, yRec, BBCInfo] = fastSSMplus(yData, mfddim, mfdorder, romorder, nforder, invstyle)
% [Mmap, iMmap, Tmap, iTmap, Nflow, yRec, BBCInfo] = fastSSMplus(yData, mfddim, mfdorder, romorder, nforder, invstyle, tanspace)
% Computes an SSM and normal form from data with SSMTool. By Joar Axas (jgoeransson@ethz.ch)
%
% INPUT:                                 OUTPUT:
% yData    {Nx2}  cell array with time   Mmap    function  SSM parametrization y=M(xi)^{1:m}
%                values (1xn) and        iMmap   function  SSM chart
%                trajectories (pxn)      Tmap    function  normal form transformation xi=T(z)
% mfddim   int    SSM dimension d        iTmap   function  inverse normal form transformation
% mfdorder int    SSM polynomial order   Nflow   function  normal form dynamics zdot=N(z)
% romorder int    reduced dynamics order yRec    (pxn)     fastSSMplus reconstruction of yData{1,2}
% nforder  int    normal form order      BBCInfo struct    computed backbone curves (only 2D)
% invstyle 'poly' optional inverse 
%       or 'mini'  transformation style (default polynomial)
% tanspace (pxd)  optional tangent space
%
% Input should have the form {[t_1(1),t_1(2),t_1(3),...], [y_1(1),y_1(2),y_1(3),...];
%                             [t_2(1),t_2(2),t_2(3),...], [y_2(1),y_2(2),y_2(3),...];
%                              ... }
% where t_i are scalar time values y_i are column p-vectors with snapshots at the respective times.
%
% EXAMPLE (data from linear system):
% t = 0:0.1:10;
% yData = {t, [sin(t);cos(t)]}
% [Mmap, iMmap, Tmap, iTmap, Nflow, yRec, BBCInfo] = fastSSMplus(yData, 2, 1, 1, 1)

%% Input, settings
if ~exist('misc/frc_ab', 'file'); error('SSMTool 2.4 not installed, cf. readme file. Clone from github.com/jain-shobhit/SSMTool'); end
t = horzcat(yData{:,1}); Y = horzcat(yData{:,2});
iStart = 1;
for iTraj = 1:size(yData,1); iStart(iTraj+1) = iStart(iTraj)+size(yData{iTraj,1},2); end
invertstyle = 'poly';
if nargin>=6; invertstyle = lower(varargin{1}); assert(strcmpi(invertstyle, 'poly')||strcmpi(invertstyle, 'mini'), "invstyle argument must be 'poly' or 'mini'"); end

%% Fit manifold
if nargin < 7 % use SVD for tangent space
    [u,s,v] = svds(Y, mfddim);
    V = max(abs(v'),[],2).'.*u*s;
else % use user-supplied matrix as tangent space
    V = varargin{2}.*max((varargin{2}\Y)');
end
iMmap = @(y) V\y;
Xi = iMmap(Y);
M = Y/phi(Xi, 1:mfdorder);
Mmap = @(xi) M*phi(xi, 1:mfdorder);

%% Compute time derivative and reduced dynamics
Xidot = []; Xinew = [];
for iTraj = 1:size(yData,1)
    [Xdot,Xnew] = ftd(Xi(:,iStart(iTraj):iStart(iTraj+1)-1),t);
    Xidot = [Xidot, Xdot]; Xinew = [Xinew, Xnew];
end
R = Xidot/phi(Xinew, 1:romorder);

%% Rearrange coefficients into tensors for SSMTool
Rcol = 1;
for order = 1:romorder
    exps = exponents(mfddim, order);
    fnl{order}.ind = exps;
    fnl{order}.coeffs = R(:,Rcol:Rcol+size(exps, 1)-1);
    Rcol = Rcol + size(exps, 1);
    fnl{order} = multi_index_to_tensor(fnl{order}.coeffs,fnl{order}.ind);  
end

%% Compute normal form with SSMTool
DS = DynamicalSystem();
set(DS.Options, 'Emax', mfddim, 'Nmax', 100, 'notation', 'multiindex')
set(DS, 'A', R(1:mfddim,1:mfddim), 'B', eye(mfddim), 'fnl', fnl(2:end));
S = SSM(DS);
set(S.Options, 'reltol', 0.1, 'notation', 'multiindex');
S.choose_E(1:mfddim)
[Tcoeff, Ncoeff] = S.compute_whisker(nforder);
Tmap = @(z) real(tpoly(Tcoeff, z));
Nflow = @(t,z) tpoly(Ncoeff, z);

%% Approximate inverse to the normal form transformation
[W, Lambda] = eig(R(1:mfddim,1:mfddim));
if strcmpi(invertstyle, 'poly') % faster
    invpoints = W\Xi;
    iT = invpoints/phi(Tmap(invpoints),1:nforder);
    iTmap = @(xi) iT*phi(xi, 1:nforder);
elseif strcmpi(invertstyle, 'mini') % more accurate
    conjtransform = @(rz) reshape([rz(1:end/2)+1i*rz(end/2+1:end), rz(1:end/2)-1i*rz(end/2+1:end)].', [], size(rz,2));
    iTmap = @(eta) conjtransform(fsolve(@(z)Tmap(conjtransform(z))-eta, zeros(size(eta)), optimset('Display','off')));
end

%% Evaluate model
options = struct('isauto', 1, 'isdamped', 1, 'numDigits', 3);
symexp = reduced_dynamics_symbolic(DS.spectrum.Lambda(1:mfddim), Ncoeff, options);
sympref('FloatingPointOutput',true);
fprintf('\\dot{\\rho}_%u =\n\\dot{\\theta_%u} =\n', ceil(0.1:1/2:mfddim/2))
disp(symexp)

[~, zRec] = ode45(Nflow, t(iStart(1):iStart(2)-1), iTmap(Xi(:,1)), odeset('RelTol', 1e-6));
yRec = Mmap(Tmap(zRec.'));

%% Backbone curves
BBCInfo = struct();
if mfddim == 2
yObservable = @(y) max(abs(y(1,:)));
rho = linspace(0,1);
theta = linspace(0,2*pi);
for ir = 1:length(rho)
    zBBC = rho(ir).*exp(1i*[theta;-theta]);
    uBBC(ir) = yObservable(Mmap(Tmap(zBBC)));
end
gamma = compute_gamma(Ncoeff);
damp = @(rho) real([Lambda(1), gamma.'])*phi(rho,0:2:nforder-1);
freq = @(rho) imag([Lambda(1), gamma.'])*phi(rho,0:2:nforder-1);
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
end

%% Subfunction: Multivariate polynomial
function exps=exponents(d,k)
    B = repmat({0:max(k)},1,d);
    A = combvec(B{:}).';
    exps = A(ismember(sum(A,2), k),:);
    [~,ind] = sort(sum(exps,2));
    exps = exps(ind,:);
end

function u = phi(xi, r) % return monomials
    x = reshape(xi, 1, size(xi, 1), []);
    exps = exponents(size(xi, 1),r);
    u = reshape(prod(x.^exps, 2), size(exps, 1), []);
end

function [dXdt, Xtrunc] = ftd(X, t) % finite time difference
    ind = 5:size(X,2)-4; Xtrunc = X(:,ind);
    dX = 4/5*(X(:,ind+1)-X(:,ind-1)) - 1/5*(X(:,ind+2)-X(:,ind-2)) + ...
        4/105*(X(:,ind+3)-X(:,ind-3)) - 1/280*(X(:,ind+4)-X(:,ind-4));
    dXdt = dX./(t(ind+1)-t(ind));
end

%% Subfunction: Express tensor as polynomial function
function xi = tpoly(T, z)
    xi = zeros(size(T(1).coeffs,1), size(z, 2));
    z = reshape(z, 1, size(z, 1), []);
    for o = 1:numel(T)
        if ~isempty(T(o).ind)
            xi = xi + T(o).coeffs*reshape(prod(z.^T(o).ind, 2), size(T(o).ind, 1), []);
        end
    end
end
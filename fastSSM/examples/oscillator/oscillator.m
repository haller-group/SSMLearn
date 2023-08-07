clearvars
close all

%% Setup
N = 2;      % number of masses
mass = 1;
stiffness = 1;
damping = 0.03;

% obsdof = [1,0,0,0]; % observe x1
obsdof = [0,1,0,0]; % observe x2
% obsdof = [0,0,1,1]; % observe sum x1dot+x2dot
% obsdof = [-1,1,0,0]; % observe sum -x1+x2 (not an embedding)
% obsdof = [1,0,0,0;
%           0,1,0,0]; % observe [x1;x2] as a vector
obs = @(x) obsdof*x;

p = 5;  % delay embedding dimension
k = 15; % timelag tau = k*dt

nicered = [0.7,0.1,0.1];
niceblue = [0.1,0.1,0.7];
nicegreen = [0.1,0.9,0.1];
nicegray = [0.6,0.6,0.6];

M = mass*eye(N);
K = stiffness*(2*eye(N) - diag(ones(1,N-1),-1) - diag(ones(1,N-1),1));
C = damping*(2*eye(N) - diag(ones(1,N-1),-1) - diag(ones(1,N-1),1));
A = [zeros(N), eye(N);
     -M\K,      -M\C];
[E,LambdaFull] = eig(A); lambdaFull = diag(LambdaFull);
[~,ind] = sort(-real(lambdaFull));
lambdaFull = lambdaFull(ind); LambdaFull = diag(lambdaFull); E = E(:,ind);
SSMDim = 2;
lambda = lambdaFull(1:SSMDim); Lambda = diag(lambda);
FF = @(x) [zeros(N,1);
            (x(1,:)-x(2,:)).^3-2*x(1,:).^2;
           -(x(1,:)-x(2,:)).^3];
F = @(t,x) A*x - FF(x);
F2 = sptensor(zeros(2*N, 2*N, 2*N));
F2(1,1,1) = -2;
F3 = sptensor(zeros(2*N, 2*N, 2*N, 2*N));
F3(1,1,1,1) = 1;
F3(1,1,1,2) = -3;
F3(1,1,2,2) = 3;
F3(1,2,2,2) = -1;
F3(2,1,1,1) = -1;
F3(2,1,1,2) = 3;
F3(2,1,2,2) = -3;
F3(2,2,2,2) = 1;
fnl = {F2,F3};
x0 = getSSMIC(M, C, K, fnl, 1, 0.225, 2, 0)
endtime = 250;
dt = 0.1;
[t,x] = ode45(F,0:dt:endtime,x0); t = t.'; x = x.';
%% Create string for observable function
for oi = 1:size(obsdof,1)
coord = 1;
obsstring = '';
for o = obsdof(oi,:)
    if o > 0
        if ~isempty(obsstring)
            obsstring = [obsstring, '+'];
        end
        if abs(o) ~= 1
            obsstring = [obsstring, num2str(abs(o))];
        end
        if coord>2
            obsstring = [obsstring, '\dot '];
        end
        obsstring = [obsstring, 'x_',num2str(mod(coord-1,2)+1),'(%s)'];
    elseif o < 0
        obsstring = [obsstring, '-'];
        if abs(o) ~= 1
            obsstring = [obsstring, num2str(abs(o))];
        end
        if coord>2
            obsstring = [obsstring, '\dot '];
        end
        obsstring = [obsstring, 'x_',num2str(mod(coord-1,2)+1),'(%s)'];
    end
    coord = coord + 1;
end
if oi==1
obsstr = @(s) ['$',strrep(obsstring, '%s', s),'$'];
elseif oi==2
obsstr2 = @(s) ['$',strrep(obsstring, '%s', s),'$'];
end
end
xDataFull = {t, x};
xData = {t, obs(x)};
q = size(xData{1,2}, 1)

%% Delay embedding
yData = embedCoordinatesBlockwise(xData, p, k);
tau = k*dt;
V = [exp(lambda.'.*(tau*(0:p-1)'))]; % Vandermonde
if q == 2
    c = [1;1]; % arbitrary nonzero rotation coefficients
    Dobs = obs(E(:,1:2)*diag(c));
    T = [V*diag(Dobs(1,:));V*diag(Dobs(2,:))] % Tangent space
elseif q == 1
    T = V
end
[Mmap, iMmap, Tmap, iTmap, Nflow, yRec] = fastSSMplus(yData, SSMDim, 5, 3, 3, 'poly', [real(T(:,1)),imag(T(:,1))]);
paperFigure('x','time','y',obsstr('t'));
plot(yData{1,2}(1,:), 'k', 'linewidth', 2, 'displayname', 'Simulation')
plot(real(yRec(1,:)), '--', 'linewidth', 2, 'color', nicegreen, 'displayname', 'Prediction')
xlim([0,2000])

%% fastSSM
[MmapX, iMmapX, TmapX, iTmapX, NflowX, xRec] = fastSSMplus(xDataFull,SSMDim,5,3,3,'poly',[real(E(:,1)),imag(E(:,1))]);

%% Plot full phase space
cap = 1;
pind = [1,3,4];
paperFigure('x',['$x_',num2str(pind(1)),'$'],'y',['$x_',num2str(pind(2)),'$'],'z',['$x_',num2str(pind(3)),'$'],'legendcols', 1);
plot2DSSM(pind, iMmapX(xDataFull{1,2}), MmapX, 10, 50, nicegray);
legend('$\mathcal{M}$')
hold on
plot3(x(pind(1),cap:end), x(pind(2),cap:end), x(pind(3),cap:end), 'k', 'linewidth', 1.2, 'displayname', 'Data');
w1 = norm(x(pind,1))*real(E(pind,1))/norm(real(E(pind,1)));
w2 = norm(x(pind,1))*imag(E(pind,1))/norm(imag(E(pind,1)));
quiver3(0,0,0,w1(1),w1(2),w1(3), 'linewidth', 3, 'MaxHeadSize', 1, 'autoscalefactor', 1.3, 'displayname', '$\mathrm{Re}\,E_1$')
quiver3(0,0,0,w2(1),w2(2),w2(3), 'linewidth', 3, 'MaxHeadSize', 1, 'autoscalefactor', 1.3, 'displayname', '$\mathrm{Im}\,E_1$')
[X,Y] = meshgrid(0.4*linspace(-1,1,2),0.4*linspace(-1,1,2));
h = surf(w1(1)*X+w2(1)*Y,w1(2)*X+w2(2)*Y,w1(3)*X+w2(3)*Y, 'displayname', '$T_0\mathcal{M}$');
h.FaceColor = nicered; h.FaceAlpha = 0.5; h.EdgeColor = nicered;
view(-57,15)

%% Plot observable space
if q == 1
    pind = [1,2,3];
    paperFigure('x',obsstr('t'),'y',obsstr(['t+',num2str(k*(pind(2)-1)),'\Delta t']),'z',obsstr(['t+',num2str(k*(pind(3)-1)),'\Delta t']),'legendcols', 1);
else
    pind = [1,2,7];
    paperFigure('x',obsstr('t'),'y',obsstr(['t+',num2str(k*(pind(2)-1)),'\Delta t']),'z',obsstr2(['t+',num2str(k*(pind(3)-1-p)),'\Delta t']),'legendcols', 1);
end
plot2DSSM(pind, iMmap(yData{1,2}), Mmap, 20, 50, nicegray);
legend('$\tilde{\mathcal{M}}$')
hold on
plot3(yData{1,2}(pind(1),cap:end), yData{1,2}(pind(2),cap:end), yData{1,2}(pind(3),cap:end), 'k', 'linewidth', 1.2, 'displayname', 'Data')
v1 = -max(vecnorm(yData{1,2}(pind,:)))*real(T(pind,1))/norm(real(T(pind,1)));
v2 = max(vecnorm(yData{1,2}(pind,:)))*imag(T(pind,1))/norm(imag(T(pind,1)));
quiver3(0,0,0,v1(1),v1(2),v1(3), 'linewidth', 3, 'MaxHeadSize', 1, 'autoscalefactor', 1, 'displayname', '$\mathrm{Re}\,V_1$')
quiver3(0,0,0,v2(1),v2(2),v2(3), 'linewidth', 3, 'MaxHeadSize', 1, 'autoscalefactor', 1, 'displayname', '$\mathrm{Im}\,V_1$')
[X,Y] = meshgrid(0.35*linspace(-1,1,2),0.35*linspace(-1,1,2));
h = surf(v1(1)*X+v2(1)*Y,v1(2)*X+v2(2)*Y,v1(3)*X+v2(3)*Y, 'displayname', '$T_0\tilde{\mathcal{M}}$');
h.FaceColor = nicered; h.FaceAlpha = 0.5; h.EdgeColor = nicered;
view(-61,6)

%% Plot reduced coordinates
paperFigure('x','$\mathrm{Re}(V^\dagger y)$','y','$\mathrm{Im}(V^\dagger y)$','legendcols',0);
zeta = T\yData{1,2};
zeta1 = real(zeta(1,:));
zeta2 = imag(zeta(1,:));
plot(zeta1, zeta2, 'color', nicered, 'linewidth', 1.6)

%%
function yData = embedCoordinatesBlockwise(xData, p, delaySteps)
    q = size(xData{1,2}, 1);
    yData1 = embedCoordinates(xData, p, delaySteps);
    yData = yData1;
    for iTraj = 1:size(yData1, 1)
        for o = 1:q
            yData{iTraj,2}(p*(o-1)+(1:p),:) = yData1{iTraj,2}(o:q:end,:);
        end
    end
end
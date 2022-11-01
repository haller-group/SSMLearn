function FRCSweep = integrateSweep(F, omegaSpan, n, fSweep, observable, varargin)

p = inputParser;
validScalar = @(x) isscalar(x);
addOptional(p, 'step', 0.1, validScalar);
addOptional(p, 'initRevolutions', 200, validScalar);
addOptional(p, 'stepRevolutions', 30, validScalar);
addOptional(p, 'odetol', 1e-4, validScalar);
parse(p, varargin{:});

opts = odeset('RelTol', p.Results.odetol);

w_0 = omegaSpan(1);
Npers = p.Results.initRevolutions;
if size(fSweep, 1) == n
    fSweep = [zeros(n,1); fSweep];
end
F_force = @(t,x,w) F(t,x) + fSweep*cos(w*t);
[t,x] = ode15s(@(t,x) F_force(t,x,w_0),[0 Npers*2*pi/w_0], zeros(2*n,1), opts);
t = t'; x = x';
figure
plot(t,observable(x))
hold on
plot(t,plot(t,observable(fSweep*cos(w_0*t))))

w_sweep = omegaSpan(1):p.Results.step:omegaSpan(2);
u_sweep = zeros(1,length(w_sweep));
phs_sweep = zeros(1,length(w_sweep));
Npers = p.Results.stepRevolutions;
for ii = 1:length(w_sweep)
    [t,x] = ode15s(@(t,x) F_force(t,x,w_sweep(ii)),[0 Npers*2*pi/w_sweep(ii)], x(:,end), opts);
    t = t'; x = x';
    [t,x] = ode15s(@(t,x) F_force(t,x,w_sweep(ii)),[0 1*2*pi/w_sweep(ii)], x(:,end), opts);
    t = t'; x = x';
    u_sweep(ii) = 0.5*abs(max(observable(x)) - min(observable(x)));
    phs_sweep(ii) = acos(observable(x)*vecnorm(fSweep)*cos(w_sweep(ii)*t')/norm(observable(x))/norm(vecnorm(fSweep)*cos(w_sweep(ii)*t)));
    disp(['computed FRC point ', num2str(ii), ' of ', num2str(length(w_sweep)), '. u = ', num2str(u_sweep(ii))])
end

FRCSweep = struct('omega',w_sweep,'amp',u_sweep,'phs',phs_sweep);
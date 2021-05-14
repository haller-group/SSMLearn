function [FRC_NI] = integrateSweep(F, w_span, n, f_sweep, observable, varargin)

p = inputParser;
validScalar = @(x) isscalar(x);
addOptional(p, 'step', 0.1, validScalar);
addOptional(p, 'initRevolutions', 200, validScalar);
addOptional(p, 'stepRevolutions', 30, validScalar);
parse(p, varargin{:});

opts = odeset('AbsTol', 1e-6);

w_0 = w_span(1);
Npers = p.Results.initRevolutions;
F_force = @(t,x,w) F(t,x) + [zeros(n,1); f_sweep*cos(w*t)];
[t_sim,x_sim] = ode15s(@(t,x) F_force(t,x,w_0),[0 Npers*2*pi/w_0], zeros(2*n,1), opts);
t_sim = t_sim'; x_sim = x_sim';
figure
plot(t_sim,observable(x_sim))

w_sweep = w_span(1):p.Results.step:w_span(2);
u_sweep = zeros(1,length(w_sweep));
Npers = p.Results.stepRevolutions;
for ii = 1:length(w_sweep)
    [t_sim,x_sim] = ode15s(@(t,x) F_force(t,x,w_sweep(ii)),[0 Npers*2*pi/w_sweep(ii)], x_sim(:,end), opts);
    t_sim = t_sim'; x_sim = x_sim';
    [t_sim,x_sim] = ode15s(@(t,x) F_force(t,x,w_sweep(ii)),[0 1*2*pi/w_sweep(ii)], x_sim(:,end), opts);
    t_sim = t_sim'; x_sim = x_sim';
    u_sweep(ii) = 0.5*abs(max(observable(x_sim)) - min(observable(x_sim)));
    disp(['computed FRC point ', num2str(ii), ' of ', num2str(length(w_sweep)), '. u = ', num2str(u_sweep(ii))])
end

FRC_NI = struct('omega',w_sweep,'amp',u_sweep);
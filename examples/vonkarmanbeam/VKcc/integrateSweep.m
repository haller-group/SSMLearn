function [FRC_NI] = integrateSweep(F, w_span, n, f_sweep, observable, varargin)

p = inputParser;
validScalar = @(x) isscalar(x);
addOptional(p, 'step', 0.1, validScalar);
addOptional(p, 'initRevolutions', 200, validScalar);
addOptional(p, 'stepRevolutions', 30, validScalar);
parse(p, varargin{:});

w_0 = w_span(1);
Npers = p.Results.initRevolutions;
F_force = @(t,x,w) F(t,x) + [zeros(n,1); f_sweep*cos(w*t)];
[t_sim,x_sim] = ode15s(@(t,x) F_force(t,x,w_0),[0 Npers*2*pi/w_0], zeros(2*n,1));
t_sim = t_sim'; x_sim = x_sim';
figure
plot(t_sim,observable(x_sim))

w_sweep = w_span(1):p.Results.step:w_span(2);
u_sweep = zeros(1,length(w_sweep));
Npers = p.Results.stepRevolutions;
for ii = 1:length(w_sweep)
    [~,x_sim] = ode15s(@(t,x) F_force(t,x,w_sweep(ii)),[0 Npers*2*pi/w_sweep(ii)], x_sim(:,end));
    x_sim = x_sim';
    [~,x_PO] = ode15s(@(t,x) F_force(t,x,w_sweep(ii)),[0 1*2*pi/w_sweep(ii)], x_sim(:,end));
    u_sweep(ii) = max(abs(observable(x_PO.')));
    disp(['computed FRC point ', num2str(ii), ' of ', num2str(length(w_sweep)), '. u = ', num2str(u_sweep(ii))])
end

FRC_NI = struct('omega',w_sweep,'amp',u_sweep);
function [t,pt] = poincareMap(x, odefun, frequency, numIterates, options)
    T = 2*pi/frequency;
    tspan = 0:T:numIterates*T;
    [t, x] = ode45(odefun, tspan, x, options);
    pt = x;
end
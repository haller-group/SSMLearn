function [F,IC] = parabolicSyst(nTraj, ICRadius, mu, omega, lambda, rot, varargin)
% Example dynamical system
% Known to have the invariant manifold z = g(x,y), with optional rotation
% Nonlinear terms can be added to the reduced dynamics x,y as an anonymous function
% Returns the anonynous evolution function and initial conditions

if isempty(varargin)
    nonlinearities = @(t,x) 0;
else
    nonlinearities = varargin{1};
end

R = @(t,x) [mu*x(1,:)+omega*x(2,:); 
            -omega*x(1,:)+mu*x(2,:)] + nonlinearities(t,x);
g = @(t,x) x(1,:).^2 + x(2,:).^2;
Dg = @(t,x) [2*x(1,:)', 2*x(2,:)'];
Q = rotation(rot(1), rot(2), rot(3));
F = @(t,x) Q*[R(t,Q'*x);
            lambda*(Q(:,3:end)'*x-g(t,Q'*x))+sum(Dg(t,Q'*x).*R(t,Q'*x)',2)'];

angle = linspace(0,2*pi,nTraj+1);
angle = angle(1:end-1);
IC = Q*[ICRadius*cos(angle);
          ICRadius*sin(angle);
          g(0, ICRadius*[cos(angle); sin(angle)])];
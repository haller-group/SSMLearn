function [t,x,idxSwitch,l0vect] = odenonsmooth_1switchsurf( f, s, Ds, timeSpan, ...
                                                        x0, varargin)
% Function that integrates mechanical models (either forced or unforced)
% which are piecewise smooth with one switching surface. The switching
% surfaces is assumed to be an hyperplane defined by a single dof.

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
% addParameter(p, 'SigTol', 1e-4, validScalarPosNum);
% addParameter(p, 'RelTol', 1e-5, validScalarPosNum);
addParameter(p, 'SigTol', 1e-6, validScalarPosNum);
addParameter(p, 'RelTol', 1e-6, validScalarPosNum);

addParameter(p, 'AbsTol', 1e-10, validScalarPosNum);
addParameter(p, 'numDofs', length(x0), validScalarPosNum);
addParameter(p, 'switchCoord', find(Ds(0)), validScalarPosNum);
addParameter(p, 'ODEsolver', @ode45);
parse(p, varargin{:});

% Define crossing function and switching dof
cross = @(t,x) ((Ds(x)*f(t,x,1))*(Ds(x)*f(t,x,-1)))>0; 

% Integrate to first switch
t0 = timeSpan(1); t = t0; x = transpose(x0);
if abs(s(x0)) > p.Results.SigTol
    l0 = sign(s(x0));
    [tseg,xseg,switched] = integrateToSwitch( f, s, l0, ...
                                                [t0 timeSpan(end)], x0, p);
    t = [t; tseg(2:end)]; x = [x; xseg(2:end,:)];
    l0vect = l0*ones(1,length(tseg(1:end)));
    if switched == 1; idxSwitch = length(t);  else; idxSwitch = []; end
else
    if sign(s(x0)) == 0
        [~,pos] = min(abs([(Ds(x)*f(t,x.',1)) (Ds(x)*f(t,x.',-1))]));
        if pos == 1; l0 = -1; else; l0 = 1; end
    else
        l0 = sign(s(x0));
    end
    idxSwitch = 1; l0vect = l0;
end
% We are on the switching surface
t0 = t(end); x0 = transpose(x(end,:)); n = length(x0)/2; 
while (t0 < timeSpan(end)) == 1
    if cross(t0,x0) == 1 % We cross the sticking surface
        l0 = -l0; % Switch
        [tseg,xseg,switched] = integrateToSwitch( f, s, l0, ...
                                                [t0 timeSpan(end)], x0, p);
        t = [t; tseg(2:end)]; x = [x; xseg(2:end,:)];
        if switched == 1; idxSwitch = [idxSwitch; length(t)]; end
        t0 = t(end); x0 = transpose(x(end,:)); l0vect = [l0vect l0*ones(1,length(tseg(2:end)))];
    else          % We stick to the switching surface
        [tseg,xseg,switched] = integrateInStick( f, s, Ds, ...
                                                [t0 timeSpan(end)], x0, p); 
        t = [t; tseg(2:end)]; x = [x; xseg(2:end,:)];
        if switched == 1; idxSwitch = [idxSwitch; length(t)]; end
        t0 = t(end); x0 = transpose(x(end,:)); %l0 = sign(x0(1));
        l0 = l0vect(end); l0vect = [l0vect l0*ones(1,length(tseg(2:end)))];
        l0 = -sign(Ds(x0)*f(t0,x0,1));
    end
end
end

function [t,x,switched] = integrateToSwitch( f, s, l0, timeSpan, x0, p)
switched = 0; t0 = timeSpan(1); t1 = timeSpan(end);
t = t0; x = transpose(x0);
while switched == 0 && abs(t1-t0)>0
    opts = odeset('RelTol', p.Results.RelTol, ...
                  'AbsTol', p.Results.AbsTol, ...
                  'Events',@(t,y) eventSwitch(t,y,s,l0));
    [tseg,xseg,te,~,~] = p.Results.ODEsolver(@(t,x) f(t,x,l0), [t0 t1], x0, opts);
    if ~isempty(te) % A switch was detected
        t = [t; tseg(2:end)]; x = [x; xseg(2:end,:)]; switched = 1;
    else % No switch detected
        t = [t; tseg(2:end)]; x = [x; xseg(2:end,:)];
        t0 = tseg(end); x0 = transpose(xseg(end,:)); t1 = t1 + (t0-t1);
    end
end
end

function [value,isterminal,direction] = eventSwitch(t,y,s,l0)
value = s(y); isterminal = 1; direction = -l0;
% value = (-l0 * s(y))>0; isterminal = 1; direction = 0;
end

function [t,x,switched] = integrateInStick( f, s, Ds, timeSpan, x0, p)
ndof = p.Results.numDofs; sdof = p.Results.switchCoord - ndof;
switched = 0; t0 = timeSpan(1); t1 = timeSpan(end); n = length(x0);
t = t0; x = transpose(x0); s0 = s(zeros(n,1)); nSplit = ndof+sdof;
xSwitch = @(y) [y(1:nSplit-1); -s0; y(nSplit:end)];
% l = @(t,y) Ds(x0)*f(t,xSwitch(y),sign(s0));
l = @(t,y) Ds(x0)*f(t,xSwitch(y),1) * Ds(x0)*f(t,xSwitch(y),-1); 
indMat = eye(n); indMat(nSplit,:) = [];
while switched == 0 && abs(t1-t0)>0
    opts = odeset('RelTol', p.Results.RelTol, ...
        'AbsTol', p.Results.AbsTol, ...
        'Events',@(t,y) eventStick(t,y,l));
    lambda_sigma = @(t,y) ((f(t,xSwitch(y),-1) + f(t,xSwitch(y),1)).' * Ds(y).')/((f(t,xSwitch(y),-1) - f(t,xSwitch(y),1)).' * Ds(y).');
%     [tseg,yseg,te,~,~] = p.Results.ODEsolver(@(t,y) indMat*f(t,...
%          xSwitch(y),0), [t0 t1], [x0(1:nSplit-1); x0(nSplit+1:end)], opts);
    [tseg,yseg,te,~,~] = p.Results.ODEsolver(@(t,y) indMat*f(t,...
         xSwitch(y),lambda_sigma(t,y)), [t0 t1], [x0(1:nSplit-1); x0(nSplit+1:end)], opts);
    xseg = [yseg(:,1:nSplit-1) -s0*ones(size(tseg)) yseg(:,nSplit:end)];
    if ~isempty(te) % A switch was detected
        t = [t; tseg(2:end)]; x = [x; xseg(2:end,:)]; switched = 1;
    else % No switch detected
        t = [t; tseg(2:end)]; x = [x; xseg(2:end,:)];
        t0 = tseg(end); x0 = transpose(xseg(end,:)); t1 = t1 + (t0-t1);
    end
end
end

function [value,isterminal,direction] = eventStick(t,y,l)
% value = l(t,y); isterminal = 1; direction = 0;
value = l(t,y)<0; isterminal = 1; direction = 0;

end


function [t,x,idxSwitch,eta] = integratePWSROM(wFunction,...
                 rFunction,vFunction,sFunctionPhys,timeSpan,x0,varargin)
% Only switching between two ROMs

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addParameter(p, 'RelTol', 1e-5, validScalarPosNum);
addParameter(p, 'AbsTol', 1e-10, validScalarPosNum);
addParameter(p, 'ODEsolver', @ode45);
addParameter(p, 'RICmethod', 'fiber');
addParameter(p, 'WeightMatrix', 1);
addParameter(p, 'gammaStoEta', @(s,l) s);
addParameter(p, 'gammaEtatoS', @(eta,l) eta);
parse(p, varargin{:});

% Define initial functions
sFunction =@(t,eta,l) sFunctionPhys(vFunction(t,eta,l));

% Determine if we start in plus or minus
t0 = timeSpan(1); t = t0; x = transpose(x0); idxSwitch = [];
l0 = sign(sFunctionPhys(x0));
if l0 == 0
    eta0 = wFunction(t0,x0,1); Dt = 1e-10; t1 = t0+Dt;
    eta1 = eta0 + Dt*rFunction(t0,eta0,1);
    x1 = vFunction(t1,eta1,1);
    l0 = sign(sFunctionPhys(x1));
end
eta0 = wFunction(t0,x0,l0); eta = transpose(eta0);
% Perform switching
while (t0 < timeSpan(end)) == 1
    [tseg,xseg,switched,etaseg] = integrateToSwitch( wFunction, ...
     rFunction, vFunction, sFunction, l0, [t0 timeSpan(end)], x0, eta0, p);
    t = [t; tseg(2:end)]; x = [x; xseg(2:end,:)]; eta0 = [];
    eta = [eta; etaseg(2:end,:)];
    if switched == 1; idxSwitch = [idxSwitch; length(t)]; end
    t0 = t(end); x0 = transpose(x(end,:)); l0 = -l0;
end
end


function [t,x,switched,eta] = integrateToSwitch( w, r, v, s, l0, ...
                                                     timeSpan, x0, eta0, p)
switched = 0; t0 = timeSpan(1); t1 = timeSpan(end); 
if isempty(eta0) == 1
    eta0 = getReducedInitialCondition(t0,l0,x0,w,v,p);
end
t = t0; x = transpose(x0); eta = transpose(eta0);
while switched == 0 && abs(t1-t0)>0
    opts = odeset('RelTol', p.Results.RelTol, ...
        'AbsTol', p.Results.AbsTol, ...
        'Events',@(t,y) eventSwitch(t,y,s,l0));
    [tseg,etaseg,te,~,~] = p.Results.ODEsolver(@(t,x) r(t,x,l0), ...
                                               [t0 t1], eta0, opts);
    xseg = transpose(v(transpose(tseg),transpose(etaseg),l0));
    t = [t; tseg(2:end)]; x = [x; xseg(2:end,:)]; 
    eta = [eta; etaseg(2:end,:)]; 
    if ~isempty(te) % A switch was detected
        switched = 1;
    else % No switch detected
        t0 = tseg(end); x0 = transpose(xseg(end,:)); t1 = t1 + (t0-t1);
        eta0 = w(t0,x0,l0);
    end
end
end

function [value,isterminal,direction] = eventSwitch(t,y,s,l)
value = s(t,y,l); isterminal = 1; direction = -l;
end

function eta0 = getReducedInitialCondition(t0,l0,x0,w,v,p)
eta0 = w(t0,x0,l0);
if strcmp(p.Results.RICmethod,'physical') == 1
    W = p.Results.WeightMatrix;
    gammaStoEta = p.Results.gammaStoEta;
    deltaS = @(s) x0-v(t0,gammaStoEta(s,l0),l0);
    s0 = fminsearch(@(s) transpose(deltaS(s))*W*deltaS(s),...
                                           p.Results.gammaEtatoS(eta0,l0)); 
    eta0 = gammaStoEta(s0,l0);
end
end
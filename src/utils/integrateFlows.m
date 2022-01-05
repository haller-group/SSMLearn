function etaRec= integrateFlows(flow, etaData)

opts = odeset('RelTol',1e-4);
nTraj = size(etaData,1);
etaRec = cell(nTraj,2);
for iTraj = 1:nTraj
    tStart = etaData{iTraj,1}(1);
    tEnd = etaData{iTraj,1}(end);
    nSamp = length(etaData{iTraj,1});
    [t, x] = ode15s(@(t,y) flow(y), linspace(tStart, tEnd, nSamp), etaData{iTraj,2}(:,1), opts);
    etaRec{iTraj,1} = t.';
    etaRec{iTraj,2} = x.';
end
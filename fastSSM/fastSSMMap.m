function [iMmap, Rflow, yRec, Mmap] = fastSSMMap(yData, mfddim, mfdorder, romorder)
% [iMmap, Rflow, yRec, Mmap = fastSSMMap(yData, mfddim, mfdorder, romorder)
% Same as fastSSMplus, but fits a map instead of a vector field on the SSM.

%% Input, settings
if length(mfdorder) == 1; mfdorder = 1:mfdorder; end
if length(romorder) == 1; romorder = 1:romorder; end
t = horzcat(yData{:,1}); y = horzcat(yData{:,2});
iStart = 1;
for iTraj = 1:size(yData,1); iStart(iTraj+1) = iStart(iTraj)+size(yData{iTraj,1},2); end

%% Dimensionality reduction
[u,s,v] = svds(y, mfddim);
V = (s\u'./max(abs(v'),[],2))';
iMmap = @(y) V'*y;
eta = iMmap(y);
M = y/phi(eta, mfdorder);
Mmap = @(eta) M*phi(eta, mfdorder);

%% Compute reduced dynamics
eta1 = []; etanew = []; timenew = [];
for iTraj = 1:size(yData,1)
    enew = eta(:,iStart(iTraj):iStart(iTraj+1)-2);
    e1 = eta(:,iStart(iTraj)+1:iStart(iTraj+1)-1);
    tnew = t(iStart(iTraj):iStart(iTraj+1)-2);
    eta1 = [eta1, e1]; etanew = [etanew, enew]; timenew = [timenew, tnew];
end

Phieta = phi(etanew, romorder);
R = eta1/Phieta;

Rflow = @(t,eta) R*phi(eta, romorder);

%% Evaluate model
etaRec(:,1) = eta(:,1);
for ii = 2:iStart(2)-1
    etaRec(:,ii) = Rflow(t(:,ii-1), etaRec(:,ii-1));
end

yRec = Mmap(etaRec);
end

%% Multivariate polynomial
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
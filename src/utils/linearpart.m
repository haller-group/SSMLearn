function [W,A,V,lambda,varargout] = linearpart(M,C,K,m,varargin)

% Default gives conservative mode shapes
n = size(M,1);
A = [sparse(n,n) speye(n,n); -M\K -M\C];
if isempty(varargin) == 1 
    % Compute conservative eigenvectors
    m = min(n,m);
    [Vcons,Dcons] = eigs(K,M,m,'smallestabs');
    dfull = diag(Dcons);
    [dfull,pos] = sort(dfull); Vcons = Vcons(:,pos);
    % Phase space modes: definition and sorting
    
    D = diag(Vcons.'*M*Vcons);
    disp(dfull)
    omega = sqrt(dfull); 
    beta = diag(Vcons.'*C*Vcons);
    zeta = (beta./(D.*omega))/2;
    lambda = [omega.*(-zeta + sqrt(zeta.^2 - 1)); omega.*(-zeta - sqrt(zeta.^2 - 1))];
    Vcons = Vcons*diag(1./sqrt(D));
    V = [Vcons sparse(n,m); sparse(n,m) Vcons];
    W = [transpose(Vcons)*M sparse(m,n); ...
        sparse(m,n) transpose(Vcons)*(M)];
%     figure; spy((Wcons*Vcons)>1e-9)
    VO = V; V(:,1:2:end) = VO(:,1:m);
    V(:,2:2:end) = VO(:,m+1:end);
    WO = W; W(1:2:end,:) = WO(1:m,:);
    W(2:2:end,:) = WO(m+1:end,:);
    varargout{1} = dfull;
    varargout{2} = zeta;
else
    % Compute full (damped) eigenvectors
    [V,D] = eig(A);
    dex = diag(D); [~,pos] = sort(abs(imag(dex))); dex = dex(pos);
    lambda = dex;
    V = V(:,pos); 
    VO = V; V = real(V); V(:,2:2:end) = imag(VO(:,1:2:end));
    W = inv(V);
end
end
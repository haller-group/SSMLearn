function [W,A,V] = linearpart(M,C,K,varargin)

% Default gives conservative mode shapes
n = size(M,1);
A = full([zeros(n) eye(n); -M\K -M\C]);

if isempty(varargin) == 1 
    % Compute conservative eigenvectors
    [Vcons,Dcons] = eig(full(K),full(M));
    dfull = diag(Dcons);
    [~,pos] = sort(dfull); Vcons = Vcons(:,pos);
    % Phase space modes: definition and sorting
    V = [Vcons zeros(n); zeros(n) Vcons];
    W = [transpose(Vcons)*full(M) zeros(n); ...
        zeros(n) transpose(Vcons)*full(M)];
%     figure; spy((Wcons*Vcons)>1e-9)
    VO = V; V(:,1:2:end) = VO(:,1:n);
    V(:,2:2:end) = VO(:,n+1:end);
    WO = W; W(1:2:end,:) = WO(1:n,:);
    W(2:2:end,:) = WO(n+1:end,:);
else
    % Compute full (damped) eigenvectors
    [V,D] = eig(A);
    dex = diag(D); [~,pos] = sort(abs(imag(dex))); dex = dex(pos);
    V = V(:,pos); 
    VO = V; V = real(V); V(:,2:2:end) = imag(VO(:,1:2:end));
    W = inv(V);
end
end
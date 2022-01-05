function [f,Df] = fMinimize(z,X,n,k,phi,phi_dim,D_phi_info,L,LI)
% Helper for IMparametrization

% Definition of cost function and its gradient for IMparametrization

% Reshape optimization vector
V = reshape(z(1:n*k),n,k);
H = reshape(z(1+n*k:end),n,phi_dim);

% Compute scalar error function
Y = transpose(V)*X;
Phi = phi(Y);

X_rec = V * Y + H * Phi;
Err = X-X_rec; 
f = (sum(sum((Err.*L).^2)) + sum((H(:).^2).*LI(:)) )/size(X,1)/size(X,2);
if nargout > 1 % gradient required
    Err_L = Err.*(L.^2);
    % Derivative with respect to V
    DFv = Err_L * transpose(Y) + X * transpose(Err_L) * V ; % Matrix: n x k
    for ii = 1:k
        Coeff_i = D_phi_info{ii,1}; ind_i = D_phi_info{ii,2};
        Di_Phi = [ Coeff_i.*[Y; Phi(1:length(Coeff_i)-k,:)]; ...
            zeros(phi_dim-length(Coeff_i),size(Y,2))];
        Di_Phi = Di_Phi(ind_i,:);
        DFv(:,ii) = DFv(:,ii) + X * transpose(sum((H*Di_Phi).*Err_L));
    end
    % Derivative with respect to W
    DFw = transpose( Phi * transpose(Err_L)); % Matrix: n x phi_dim
    Df = -2/size(X,1)/size(X,2) * [reshape(DFv,n*k,1) ; ...
        reshape(DFw,n*phi_dim,1)-H(:).*LI(:)];
end
end
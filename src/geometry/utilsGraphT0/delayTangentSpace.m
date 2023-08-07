function V = delayTangentSpace(tau, p, lambda)
% V = delayTangentSpace(tau, p, lambda)
%
% Compute the eigenvectors of the linear part of a delay embedded system
% given the eigenvalues \lambda_j. 
%
% tau = k*dt

lambda = reshape(lambda, 1, []);

V1 = exp(tau*(0:p-1)'.*lambda);
indImPlus = imag(lambda) > 1e-10;
indImMinus = imag(lambda) < -1e-10;
V = V1;
V(:,indImPlus) = real(V1(:,indImPlus));
V(:,indImMinus) = imag(V1(:,indImPlus));
% V = V./max(abs(V));
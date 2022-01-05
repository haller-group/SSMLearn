function [NCtoC,CtoNC] = rcoordinatesStandardization(Q)
% Standardize coordinates: rotation and normalization according to SVD
% assuming that data Q is centered on the origin
[U,~,~] = svds(Q,size(Q,1));
QR = U'*Q; A = max(QR,[],2); A_1 = A.^(-1);
CtoNC = @(q) diag(A_1)*U'*q; 
NCtoC = @(y) U*diag(A)*y; 
end
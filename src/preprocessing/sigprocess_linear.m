function [omega,lambda,A,t_eig] = sigprocess_linear(X,t,window_length,window_overlap)

X = transpose(X(:,1:end-1));
Y = transpose(X(:,2:end));

lambda = []; A = []; t_eig = [];
ind_end = window_length + 1;
stepping = round(window_length*(1-window_overlap));
while ind_end < size(X,1)
    ind_i  = ind_end-window_length:ind_end;
    X_i = X(ind_i,:);  Y_i = Y(ind_i,:);
    A_i = transpose((transpose(X_i)*X_i) \ (transpose(X_i)*Y_i));
    [~,D_i] = eig(A_i);
    lambda_i = diag(D_i);
    [~,pos] = sort(abs(lambda_i),'ascend');
    lambda = [lambda lambda_i(pos)]; 
    A = [A A_i(:)]; 
    t_eig = [t_eig (t(ind_i(1))+t(ind_i(end)))/2];
    ind_end = ind_end + stepping ;
end
omega = log(lambda)/(t(2)-t(1)); 
end


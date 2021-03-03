function V_M = multivariate_polynomial_lintransf(V,k,M)
% Given phi, a k-variate polynomial of order M, and a transformation x = V*y,
% this code gets the matrix V_M such that phi(x) = phi(V*y) = V_M * phi(y)
Exp_mat = []; phi_sym = [];
phi_info = cell(M,2); u = sym('u',[1 k]);
for ii = 1 : M
Exp_mat_i = multivariante_exponents(k,ii);
phi_sym_i = prod(u.^Exp_mat_i,2);
phi_info{ii,1} = phi_sym_i; 
phi_info{ii,2} = matlabFunction(phi_sym_i,'Vars', {transpose(u)} ); 
phi_info{ii,3} = Exp_mat_i;
end
  
V_M = V; 
for ii = 2:M
phi_i = phi_info{ii,2}; Exp_mat_i = phi_info{ii,3};
phi_Vi = phi_i(V*transpose(Exp_mat_i));
phi_i = phi_i(transpose(Exp_mat_i)); 
V_M_i = phi_Vi/phi_i;
V_M = sparse(blkdiag(V_M,V_M_i));
end

% % Optional Check
% phi = @(x) x;
% for ii = 2:M
% phi_i = phi_info{ii,2}; 
% phi = @(x) [phi(x); phi_i(x)]; 
% end
% % phi = matlabFunction(phi_sym,'Vars', {transpose(u)} );
% x_eval = 10*((rand(k,10000)*2-1)+1i*(rand(k,10000)*2-1));
% a1 = V_M*phi(x_eval); a2 = phi(V*x_eval);
% a3 = sqrt(sum((a1-a2).*conj(a1-a2))); a3_n = sqrt(sum((a2).*conj(a2)));
% a4 = a3./a3_n;
% rel_err_per = mean(a4)*100;
% figure(943); clf; hold on; grid on;
% plot(1:length(a3),a4,'Linewidth',2)
% set(gca,'YScale','log'); 
% xlabel('Sample #')
% ylabel('Error')
end
function l_est_cont = computeEigenvaluesMap(Maps_info, dt)

R_coeff = Maps_info.R.coeff;
SSM_dim = size(R_coeff, 1);
R_coeff_lin = R_coeff(:,1:SSM_dim);
l_est = eig(R_coeff_lin);
l_est_cont = log(l_est)/dt;

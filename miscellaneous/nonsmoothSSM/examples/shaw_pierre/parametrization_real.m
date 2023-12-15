function [z1, z2, coeff_param, coeff_dyn] = parametrization_real(Lambda, V, q0, alpha_num, lambda, A_tilde)
%dare in input D_tilde e V_tilde
    
    lambda_diag = diag(Lambda);
    lambda_1 = lambda_diag(1);
    lambda_1_conj = lambda_diag(3);
    lambda_2 = lambda_diag(2);
    lambda_2_conj = lambda_diag(4);
    
    inv_V = inv(V);
    r = inv_V(:,3);
    p_vect = V(1,:);
    
    %order 2
    P_2 = [2 * real(lambda_1) - real(lambda_2), -imag(lambda_1), 0, - imag(lambda_2), 0, 0;
           2 * imag(lambda_1), 2 * real(lambda_1) - real(lambda_2), -2 * imag(lambda_1), 0, -imag(lambda_2), 0;
           0, imag(lambda_1), 2 * real(lambda_1) - real(lambda_2), 0, 0, -imag(lambda_2);
           imag(lambda_2), 0, 0, 2 * real(lambda_1) - real(lambda_2), -imag(lambda_1), 0;
           0, imag(lambda_2), 0, 2 * imag(lambda_1), 2 * real(lambda_1) - real(lambda_2), -2 * imag(lambda_1);
           0, 0, imag(lambda_2), 0, imag(lambda_1), 2 * real(lambda_1) - real(lambda_2)];

    b_2 = lambda * [r(3) * 3 * alpha_num * q0 * p_vect(1)^2;
                           r(3) * 6 * alpha_num * q0 * p_vect(1) * p_vect(2);
                           r(3) * 3 * alpha_num * q0 * p_vect(2)^2;
                           r(4) * 3 * alpha_num * q0 * p_vect(1)^2;
                           r(4) * 6 * alpha_num * q0 * p_vect(1) * p_vect(2);
                           r(4) * 3 * alpha_num * q0 * p_vect(2)^2];

    coeff_2 = P_2\b_2;

    c_3 = coeff_2(1);
    c_4 = coeff_2(2);
    c_5 = coeff_2(3);
    d_3 = coeff_2(4);
    d_4 = coeff_2(5);
    d_5 = coeff_2(6);
    %order 3
    P_3 = [3 * real(lambda_1) - real(lambda_2), -imag(lambda_1), 0, 0, -imag(lambda_2), 0, 0, 0;
           3 * imag(lambda_1), 3 * real(lambda_1) - real(lambda_2), -2 * imag(lambda_1), 0, 0, -imag(lambda_2), 0, 0;
           0, 2 * imag(lambda_1), 3 * real(lambda_1) - real(lambda_2), -3 * imag(lambda_1), 0, 0, -imag(lambda_2),0;
           0, 0, imag(lambda_1),3 * real(lambda_1) - real(lambda_2),0, 0, 0, -imag(lambda_2);
           imag(lambda_2), 0, 0, 0, 3 * real(lambda_1) - real(lambda_2), -imag(lambda_1), 0, 0;
           0, imag(lambda_2), 0, 0, 3 * imag(lambda_1), 3 * real(lambda_1) - real(lambda_2), -2 * imag(lambda_1), 0;
           0, 0, imag(lambda_2), 0, 0, 2 * imag(lambda_1), 3 * real(lambda_1) - real(lambda_2), -3 * imag(lambda_1);
           0, 0, 0, imag(lambda_2), 0,0,imag(lambda_1), 3 * real(lambda_1) - real(lambda_2)];

    b_3 = [-r(3) * alpha_num * p_vect(1)^3 + lambda * r(3) * (6 * alpha_num * q0 * (c_3 * p_vect(3) + d_3 * p_vect(4)) * p_vect(1)) - lambda * 3 * alpha_num * p_vect(1)^2 * q0 * (2 * c_3 * r(1) + c_4 * r(2));
           -r(3) * 3 * alpha_num * p_vect(1)^2 * p_vect(2) + lambda * r(3) * (3 * alpha_num * q0 * (2 * p_vect(2) * (c_3 * p_vect(3) + d_3 * p_vect(4)) + 2 * p_vect(1) * (c_4 * p_vect(3) + d_4 * p_vect(4)))) - lambda * 3 * alpha_num * p_vect(1)^2 * q0 * (c_4 * r(1) + 2 * c_5 * r(2)) - lambda * 6 * alpha_num * p_vect(1) * p_vect(2) * q0 * (2 * c_3 * r(1) + c_4 * r(2)); 
           -r(3) * 3 * alpha_num * p_vect(1) * p_vect(2)^2 + lambda * r(3) * (3 * alpha_num * q0 * (2 * p_vect(2) * (c_4 * p_vect(3) + d_4 * p_vect(4)) + 2 * p_vect(1) * (c_5 * p_vect(3) + d_5 * p_vect(4)))) - lambda * 3 * alpha_num * p_vect(2)^2 * q0 * (2 * c_3 * r(1) + c_4 * r(2)) - lambda * 6 * alpha_num * p_vect(1) * p_vect(2) * q0 * (c_4 * r(1) + 2 * c_5 * r(2));
           -r(3) * alpha_num * p_vect(2)^3 + lambda * r(3) * (6 * alpha_num * q0 * (c_5 * p_vect(3) + d_5 * p_vect(4)) * p_vect(2)) - lambda * 3 * alpha_num * p_vect(2)^2 * q0 * (c_4 * r(1) + 2 * c_5 * r(2));
           -r(4) * alpha_num * p_vect(1)^3 + lambda * r(4) * (6 * alpha_num * q0 * (c_3 * p_vect(3) + d_3 * p_vect(4)) * p_vect(1)) - lambda * 3 * alpha_num * p_vect(1)^2 * q0 * (2 * d_3 * r(1) + d_4 * r(2));
           -r(4) * 3 * alpha_num * p_vect(1)^2 * p_vect(2) + lambda * r(4) * (3 * alpha_num * q0 * (2 * p_vect(2) * (c_3 * p_vect(3) + d_3 * p_vect(4)) + 2 * p_vect(1) * (c_4 * p_vect(3) + d_4 * p_vect(4)))) - lambda * 3 * alpha_num * p_vect(1)^2 * q0 * (d_4 * r(1) + 2 * d_5 * r(2)) - lambda * 6 * alpha_num * p_vect(1) * p_vect(2) * q0 * (2 * d_3 * r(1) + d_4 * r(2));
           -r(4) * 3 * alpha_num * p_vect(1) * p_vect(2)^2 + lambda * r(4) * (3 * alpha_num * q0 * (2 * p_vect(2) * (c_4 * p_vect(3) + d_4 * p_vect(4)) + 2 * p_vect(1) * (c_5 * p_vect(3) + d_5 * p_vect(4)))) - lambda * 3 * alpha_num * p_vect(2)^2 * q0 * (2 * d_3 * r(1) + d_4 * r(2)) - lambda * 6 * alpha_num * p_vect(1) * p_vect(2) * q0 * (d_4 * r(1) + 2 * d_5 * r(2));
           -r(4) * alpha_num * p_vect(2)^3 + lambda * r(4) * (6 * alpha_num * q0 * (c_5 * p_vect(3) + d_5 * p_vect(4)) * p_vect(2)) - lambda * 3 * alpha_num * p_vect(2)^2 * q0 * (d_4 * r(1) + 2 * d_5 * r(2))];

    coeff_3 = P_3 \ b_3;

    c_6 = coeff_3(1);
    c_7 = coeff_3(2);
    c_8 = coeff_3(3);
    c_9 = coeff_3(4);

    d_6 = coeff_3(5);
    d_7 = coeff_3(6);
    d_8 = coeff_3(7);
    d_9 = coeff_3(8);

    %%%%%%%%%%%%%%%%%%%%% parametrization

    z1 = @(Y) c_3 * Y(1).^2 + c_4 * Y(1).*Y(2) + c_5 * Y(2).^2 + c_6*Y(1).^3 + c_7 * Y(1).^2.*Y(2) + c_8 * Y(1).*Y(2).^2 + c_9 * Y(2).^3;
    z2 = @(Y) d_3 * Y(1).^2 + d_4 * Y(1).*Y(2) + d_5 * Y(2).^2 + d_6*Y(1).^3 + d_7 * Y(1).^2.*Y(2) + d_8 * Y(1).*Y(2).^2 + d_9 * Y(2).^3;
    
   
    coeff_param = [c_3 c_4 c_5 c_6 c_7 c_8 c_9;
                   d_3 d_4 d_5 d_6 d_7 d_8 d_9];
    

    %%%%%%%%%%%%%%%%%%%% coefficients reduced dynamics

    e_1 = A_tilde(1,1);
    e_2 = A_tilde(1,2);
    e_3 = lambda * 3 * r(1) * alpha_num * q0 * p_vect(1)^2;
    e_4 = lambda * 6 * r(1) * alpha_num * q0 * p_vect(1) * p_vect(2);
    e_5 = lambda * 3 * r(1) * alpha_num * q0 * p_vect(2)^2;
    e_6 = (lambda * 6 * r(1) * alpha_num * p_vect(1) * q0 * (c_3 * p_vect(3) + d_3 * p_vect(4)) - r(1) * alpha_num * p_vect(1)^3);
    e_7 = (lambda * 3 * r(1) * alpha_num * q0 * (2 * p_vect(2) * (c_3 * p_vect(3) + d_3 * p_vect(4)) + 2 * p_vect(1) * (c_4 * p_vect(3) + d_4 * p_vect(4))) - 3 * r(1) * alpha_num * p_vect(1)^2 * p_vect(2));
    e_8 = (lambda * 3 * r(1) * alpha_num * q0 * (2 * p_vect(2) * (c_4 * p_vect(3) + d_4 * p_vect(4)) + 2 * p_vect(1) * (c_5 * p_vect(3) + d_5 * p_vect(4))) - 3 * r(1) * alpha_num * p_vect(1) * p_vect(2)^2);
    e_9 = (lambda * 6 * r(1) * alpha_num * p_vect(2) * q0 * (c_5 * p_vect(3) + d_5 * p_vect(4)) - r(1) * alpha_num * p_vect(2)^3);

    f_1 = A_tilde(2,1);
    f_2 = A_tilde(2,2);
    f_3 = lambda * 3 * r(2) * alpha_num * q0 * p_vect(1)^2;
    f_4 = lambda * 6 * r(2) * alpha_num * q0 * p_vect(1) * p_vect(2);
    f_5 = lambda * 3 * r(2) * alpha_num * q0 * p_vect(2)^2;
    f_6 = (lambda * 6 * r(2) * alpha_num * p_vect(1) * q0 * (c_3 * p_vect(3) + d_3 * p_vect(4)) - r(2) * alpha_num * p_vect(1)^3);
    f_7 = (lambda * 3 * r(2) * alpha_num * q0 * (2 * p_vect(2) * (c_3 * p_vect(3) + d_3 * p_vect(4)) + 2 * p_vect(1) * (c_4 * p_vect(3) + d_4 * p_vect(4))) - 3 * r(2) * alpha_num * p_vect(1)^2 * p_vect(2));
    f_8 = (lambda * 3 * r(2) * alpha_num * q0 * (2 * p_vect(2) * (c_4 * p_vect(3) + d_4 * p_vect(4)) + 2 * p_vect(1) * (c_5 * p_vect(3) + d_5 * p_vect(4))) - 3 * r(2) * alpha_num * p_vect(1) * p_vect(2)^2);
    f_9 = (lambda * 6 * r(2) * alpha_num * p_vect(2) * q0 * (c_5 * p_vect(3) + d_5 * p_vect(4)) - r(2) * alpha_num * p_vect(2)^3);


    coeff_dyn = [e_1 e_2 e_3 e_4 e_5 e_6 e_7 e_8 e_9;
                 f_1 f_2 f_3 f_4 f_5 f_6 f_7 f_8 f_9];

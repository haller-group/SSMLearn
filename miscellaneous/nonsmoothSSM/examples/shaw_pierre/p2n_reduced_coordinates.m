function eta_minus = p2n_reduced_coordinates(eta_plus, q0, V_tilde)

y1_plus = eta_plus(1);
y2_plus = eta_plus(2);

x0_plus = [-q0; -q0/2; 0; 0];
x0_minus = [q0; q0/2; 0; 0];

inv_V = inv(V_tilde);

y1_minus = y1_plus + inv_V(1,:) * (x0_plus - x0_minus);
y2_minus = y2_plus + inv_V(2,:) * (x0_plus - x0_minus);

eta_minus = [y1_minus; y2_minus];
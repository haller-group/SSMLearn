function dy = shawpierre(t,x,  c1, c2, k, gamma, amplitude, omega)
    A = [[0, 1, 0, 0]; [-2*k, -c1 - c2, k, c1]; [0, 0,0,1]; [k, c2, -2*k, -c1 - c2]];
    dy = A*x + [0; -gamma*x(1).^3 + amplitude * cos(omega * t); 0 ; 0];
end
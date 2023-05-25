function dy = planar(t, x, a, b, c)
    dy = [x(1)*(x(2)-b); c*x(2)*(x(1)-a)];
end
function w = timeWeighting(c1, c2, t)

w = 1 ./ (1 + c1*exp(-c2*t));
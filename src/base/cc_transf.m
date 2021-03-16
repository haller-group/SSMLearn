function [X] = cc_transf(x)
X = [x; conj(x)];
end
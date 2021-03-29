function [X] = reim_transf(x)
x = x(1:size(x,1)/2,:);
X = [real(x); imag(x)];
end
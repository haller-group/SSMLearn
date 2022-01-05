function X = transformationReIm(x)
x = x(1:size(x,1)/2,:);
X = [real(x); imag(x)];
end
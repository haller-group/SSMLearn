function y = iterateMap(R,iters,y0)
y = zeros(length(y0),iters); y(:,1) = y0;
for ii = 2:iters
y(:,ii) = R(y(:,ii-1));
end
end

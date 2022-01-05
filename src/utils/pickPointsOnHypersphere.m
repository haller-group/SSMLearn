function x = pickPointsOnHypersphere(n, d, varargin)
% x = pickPointsOnHypersphere(n, d)
% x = pickPointsOnHypersphere(n, d, seed)
% Returns n points on a d-1-sphere in R^d
% The algorithm tries to maximize the distance between the points
% The seed for the random number generator can be set in the optional
% argument for consistent results.

if ~isempty(varargin) > 0
    rng(varargin{1})
end

x = normrnd(0,1,[d,max(300,3*n)]);
x = x./vecnorm(x);

while size(x,2) > n
    D = x'*x - eye(size(x,2));
    [~, closestIndex] = max(max(D));
    x(:, closestIndex) = [];
end

% Plot points on a 1- or 2-sphere:
% plot(x(1,:),x(2,:),'ko')
% axis equal

% sphere(200);
% set(findobj(gcf, 'type', 'surface'), 'EdgeColor','none')
% hold on
% plot3(x(1,:),x(2,:),x(3,:),'k.')
% axis equal

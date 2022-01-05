function c = pinkgreen(m)
%pinkgreen    Shades of pink and green color map
%   pinkgreen(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright green, range through shades of green
%    to white, and then through shades of pick to bright pick.
%   pinkgreen, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(pinkgreen)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.
map_base = [208 28 139; 241 182 218; 255 255 255; 184 225 134; 77 172 38]/255;
d_shape = 0.125;
vec_base = [0 0.5-d_shape 0.5 0.5+d_shape 1];
c = interp1(vec_base,map_base,[0:(m-1)]'/(m-1));

% if nargin < 1, m = size(get(gcf,'colormap'),1); end
% if (mod(m,2) == 0)
%     % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
%     m1 = m*0.5;
%     r = (0:m1-1)'/max(m1-1,1);
%     g = r;
%     r = [r; ones(m1,1)];
%     g = [g; flipud(g)];
%     b = flipud(r);
% else
%     % From [0 0 1] to [1 1 1] to [1 0 0];
%     m1 = floor(m*0.5);
%     r = (0:m1-1)'/max(m1,1);
%     g = r;
%     r = [r; ones(m1+1,1)];
%     g = [g; 1; flipud(g)];
%     b = flipud(r);
% end
% c = [r g b]; 
function w = timeWeighting(c1, c2, t, varargin)

w = 1 ./ (1 + c1*exp(-c2*t));

if ~isempty(varargin)
    plot(t, timeWeighting(c1,c2,t), 'r', 'DisplayName', 'weighting')
end
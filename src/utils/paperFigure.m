function fig = paperFigure(varargin)
% fig = customFigure()
% Open a figure with default settings for a nice-looking article.
%
% INPUT
% subPlot     [rows,cols]  Optional subplts.
% OUTPUT
% fig         figure handle
%

fontName = 'times';
fontSize = 30;

p = inputParser;
addOptional(p, 'x', '');
addOptional(p, 'y', '');
addOptional(p, 'z', '');
addOptional(p, 'legendcols', 1);
parse(p, varargin{:});

fig = figure('position', [100,100,560,420]); clf;
hold on; grid on; box on;
set(gca,'defaulttextinterpreter','latex');
set(gca, 'fontname', fontName);
set(gca, 'fontsize', fontSize);
if ~isempty(p.Results.x)
    xlabel(p.Results.x);
end
if ~isempty(p.Results.y)
    ylabel(p.Results.y);
end
if ~isempty(p.Results.z)
    zlabel(p.Results.z);
end
if p.Results.legendcols>0
    legend('location', 'best', 'interpreter', 'latex', 'numcolumns', p.Results.legendcols);
end
function fig = customFigure(varargin)
% fig = customFigure()
% Open a figure with default settings. Is run in each plot function of
% SSMLearn. Update here to change plot settings globally.
%
% INPUT
% figNumber   int>=1       Optional figure number.
% subPlot     [rows,cols]  Optional subplts.
% OUTPUT
% fig         figure handle
%
% Default commands included:

%   clf; hold on; grid on; box on;
%   fontName = 'helvetica';
%   fontSize = 18;

fontName = 'helvetica';
fontSize = 18;

p = inputParser;
addOptional(p, 'figNumber', []);
addOptional(p, 'subPlot', []);
parse(p, varargin{:});

if isempty(p.Results.figNumber) == 1
    fig = figure; clf;
else
    fig = figure(p.Results.figNumber); clf;
end

if isempty(p.Results.subPlot) == 1
    hold on; grid on; box on;
    set(gca, 'fontname', fontName)
    set(gca, 'fontsize', fontSize)
else
    rows = p.Results.subPlot(1);
    cols = p.Results.subPlot(2);
    if rows*cols > 4; fontSize = 12; end
    if rows*cols > 40; fontSize = 6; end
    for idxPlot = 1:rows*cols
        subplot(rows,cols,idxPlot); hold on; grid on; box on;
        set(gca, 'fontname', fontName)
        set(gca, 'fontsize', fontSize)
    end
end
end

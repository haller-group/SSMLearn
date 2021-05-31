function plotFRC(FRC, varargin)

p = inputParser;
validString = @(x) isstring(x)||ischar(x);
validStringOrCell = @(x) validString(x)||iscell(x);
validMode = @(x) strcmp(x,'Amplitude') || strcmp(x,'Phase');
addOptional(p, 'color', [0,0,0]);
addOptional(p, 'datalabel', 'FRC', validStringOrCell);
addParameter(p, 'y', 'Amplitude', validMode);
addParameter(p, 'curves', 1:length(fieldnames(FRC)));
parse(p, varargin{:});
oneLabelOnly = 0;
if validString(p.Results.datalabel)
    datalabel = {p.Results.datalabel};
    oneLabelOnly = 1;
else
    datalabel = p.Results.datalabel;
end

% Plot
hold on; grid on; box on;
for ii = p.Results.curves
    freq_i = FRC.(['F' num2str(ii)]).Freq;
    if strcmp(p.Results.y, 'Amplitude')
        plot_i  = FRC.(['F' num2str(ii)]).Amp;
    elseif strcmp(p.Results.y, 'Phase')
        plot_i  = FRC.(['F' num2str(ii)]).Nf_Phs;
    end
    stab_i = FRC.(['F' num2str(ii)]).Stab;
    [~,pos] = find(abs(diff(stab_i))==1);
    if isempty(pos)==1
%         h_i = plot(freq_i,amp_i,'Color',p.Results.color(ii,:),'Linewidth',2,...
%             'DisplayName', datalabel{ii});
        h_i = plot(freq_i,plot_i,'Color',p.Results.color(ii,:),'Linewidth',2,...
            'DisplayName', datalabel{ii});
    else
        h_i = plot(freq_i(1:pos(1)),plot_i(1:pos(1)),'Color',p.Results.color(ii,:),'Linewidth',2,...
            'DisplayName', datalabel{ii});
        if length(pos)>1
            h_ii = plot(freq_i(pos(1)+1:pos(2)),plot_i(pos(1)+1:pos(2)),'--','Color',p.Results.color(ii,:),'Linewidth',2,...
                'DisplayName', ['FRC unstable -  ', datalabel{ii}]);
            h_iii = plot(freq_i(pos(2)+1:end),plot_i(pos(2)+1:end),'Color',p.Results.color(ii,:),'Linewidth',2);
            h_iii.Annotation.LegendInformation.IconDisplayStyle = 'off';
            h_ii.Annotation.LegendInformation.IconDisplayStyle = 'off';
            
        else
            h_ii = plot(freq_i(pos(1)+1:end),plot_i(pos(1)+1:end),'--','Color',p.Results.color(ii,:),'Linewidth',2,...
                'DisplayName', ['FRC unstable - ',datalabel{ii}]);
        end
    end
    if oneLabelOnly
    if ii~= 1; h_i.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h_ii.Annotation.LegendInformation.IconDisplayStyle = 'off'; end
    end
end

xlabel('$\Omega \, [$rad/s$]$','Interpreter','latex')
ylabel('$u \, [$m$]$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
legend('location','NW')

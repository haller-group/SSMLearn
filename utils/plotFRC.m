function plotFRC(FRC, varargin)

p = inputParser;
validString = @(x) isstring(x)||ischar(x);
addOptional(p, 'color', 'b');
addOptional(p, 'datalabel', '', validString);
parse(p, varargin{:});

% Plot
hold on; grid on; box on;
for ii = 1:length(fieldnames(FRC))
    freq_i = FRC.(['F' num2str(ii)]).Freq;
    amp_i  = FRC.(['F' num2str(ii)]).Amp;
    stab_i = FRC.(['F' num2str(ii)]).Stab;
    [~,pos] = find(abs(diff(stab_i))==1);
    if isempty(pos)==1
%         h_i = plot(freq_i,amp_i,'Color',p.Results.color,'Linewidth',2,...
%             'DisplayName', ['FRC stable - ', p.Results.datalabel]);
        h_i = plot(freq_i,amp_i,'Color',p.Results.color,'Linewidth',2,...
            'DisplayName', ['FRC - ', p.Results.datalabel]);
    else
        h_i = plot(freq_i(1:pos(1)),amp_i(1:pos(1)),'Color',p.Results.color,'Linewidth',2,...
            'DisplayName', ['FRC - ', p.Results.datalabel]);
        if length(pos)>1
            h_ii = plot(freq_i(pos(1)+1:pos(2)),amp_i(pos(1)+1:pos(2)),'--','Color',p.Results.color,'Linewidth',2,...
                'DisplayName', ['FRC unstable -  ', p.Results.datalabel]);
            h_iii = plot(freq_i(pos(2)+1:end),amp_i(pos(2)+1:end),'Color',p.Results.color,'Linewidth',2);
            h_iii.Annotation.LegendInformation.IconDisplayStyle = 'off';
            h_ii.Annotation.LegendInformation.IconDisplayStyle = 'off';
            
        else
            h_ii = plot(freq_i(pos(1)+1:end),amp_i(pos(1)+1:end),'--','Color',p.Results.color,'Linewidth',2,...
                'DisplayName', ['FRC unstable - ',p.Results.datalabel]);
        end
    end
    if ii~= 1; h_i.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h_ii.Annotation.LegendInformation.IconDisplayStyle = 'off'; end
end

xlabel('$\Omega \, [$rad/s$]$','Interpreter','latex')
ylabel('$u \, [$m$]$','Interpreter','latex')
set(gca,'fontname','helvetica')
set(gca,'fontsize',18)
legend('location','NW')

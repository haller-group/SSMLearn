function plotFRC2DSSM(FRC, varargin)

p = inputParser;
validString = @(x) isstring(x)||ischar(x);
validStringOrCell = @(x) validString(x)||iscell(x);
validMode = @(x) strcmp(x,'Amplitude') || strcmp(x,'Phase');
addOptional(p, 'color', [0,0,0]);
addOptional(p, 'datalabel', '', validStringOrCell);
addParameter(p, 'y', 'Amplitude', validMode);
addParameter(p, 'curves', 1:length(fieldnames(FRC)));
addParameter(p, 'freqscale', 1); % pass 2*pi to get hertz
parse(p, varargin{:});
if validString(p.Results.datalabel)
    datalabel = {p.Results.datalabel};
else
    datalabel = p.Results.datalabel;
end

rep = [0; 0];
% Plot
hold on; grid on; box on;
c_idx = 0;
for jj = p.Results.curves
    if size(p.Results.color,1) == 1
        c_idx = 1;
    else
        rep = [0; 0];
        c_idx = c_idx + 1;
    end
    IDnamej = datalabel{c_idx};
    Freq_i = FRC.(['F' num2str(jj)]).Freq;
    if strcmp(p.Results.y, 'Amplitude')
        Plot_i  = FRC.(['F' num2str(jj)]).Amp;
    elseif strcmp(p.Results.y, 'Phase')
        Plot_i  = FRC.(['F' num2str(jj)]).Nf_Phs;
    end
    Stab_i = FRC.(['F' num2str(jj)]).Stab;
    for ll = 1:size(Freq_i,1)
        freq_i = Freq_i(ll,:);
        y_i  = Plot_i(ll,:);
        stab_i = Stab_i(ll,:);
        % Analyze Stab Changes
        [~,pos] = find([1 abs(diff(stab_i))==1]);
        for ii = 1:length(pos)
            if ii == length(pos)
                rep = stableunstableplot(freq_i(pos(ii):end),y_i(pos(ii):end),...
                    p.Results.color(c_idx,:),stab_i(pos(ii)),IDnamej,rep,p);
            else
                rep = stableunstableplot(freq_i(pos(ii):(pos(ii+1)-1)),...
                    y_i(pos(ii):(pos(ii+1)-1)),...
                    p.Results.color(c_idx,:),stab_i(pos(ii)),IDnamej,rep,p);
            end
        end
    end
end

xlabel('$\Omega$','Interpreter','latex')
ylabel('Amplitude','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
legend('location','best')

end

function [rep] = stableunstableplot(freq_i,amp_i,IDcolor,stable,IDname,rep,p)
if stable == 1
    h = plot(freq_i/p.Results.freqscale,amp_i,'Color',IDcolor,'Linewidth',2,...
        'DisplayName', ['FRC stable ' IDname]);
    if rep(1) > 0
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    rep(1) = rep(1)+1;
else
    h = plot(freq_i/p.Results.freqscale,amp_i,'--','Color',IDcolor,'Linewidth',2,...
        'DisplayName', ['FRC unstable ' IDname]);
    if rep(2) > 0
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    rep(2) = rep(2)+1;
end
end
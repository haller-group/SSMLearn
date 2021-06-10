function plotFRC2DSSM(FRC_data,f_red,color,varargin)
IDname = varargin{:};
rep = [0; 0];
for jj = 1:length(f_red)
    if size(color,1) ~= length(f_red)
        c_idx = 1; IDnamej = IDname;
    else
        rep = [0; 0];
        c_idx = jj; IDnamej = [', f = ' num2str(f_red(jj)) IDname];
    end
    Freq_i = FRC_data.(['F' num2str(jj)]).Freq;
    Amp_i  = FRC_data.(['F' num2str(jj)]).Amp;
    Stab_i = FRC_data.(['F' num2str(jj)]).Stab;
    for ll = 1:size(freq_i,1)
        freq_i = Freq_i(ll,:);
        amp_i  = Amp_i(ll,:);
        stab_i = Stab_i(ll,:);
        % Analyze Stab Changes
        [~,pos] = find([1 abs(diff(stab_i))==1]);
        for ii = 1:length(pos)
            if ii == length(pos)
                rep = stableunstableplot(freq_i(pos(ii):end),amp_i(pos(ii):end),...
                    color(c_idx,:),stab_i(pos(ii)),IDnamej,rep);
            else
                rep = stableunstableplot(freq_i(pos(ii):(pos(ii+1)-1)),...
                    amp_i(pos(ii):(pos(ii+1)-1)),...
                    color(c_idx,:),stab_i(pos(ii)),IDnamej,rep);
            end
        end
    end
end
end

function [rep] = stableunstableplot(freq_i,amp_i,IDcolor,stable,IDname,rep)
if stable == 1
    h = plot(freq_i/2/pi,amp_i,'Color',IDcolor,'Linewidth',2,...
        'DisplayName', ['FRC stable' IDname]);
    if rep(1) > 0
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    rep(1) = rep(1)+1;
else
    h = plot(freq_i/2/pi,amp_i,'--','Color',IDcolor,'Linewidth',2,...
        'DisplayName', ['FRC unstable' IDname]);
    if rep(2) > 0
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    rep(2) = rep(2)+1;
end
end
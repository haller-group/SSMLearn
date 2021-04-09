function Errors = ampph_error(x,x_ref,t_ref,varargin)

% Split x_ref in buckets
[~,locs]= findpeaks(abs(x_ref));
if size(locs,2)==1; locs = transpose(locs); end
idx = [1 floor((locs(1:end-1)+locs(2:end))/2) length(x_ref)];
if isempty(varargin) == 0
clf; hold on; grid on; box on;
plot(t_ref,x_ref,'k','Linewidth',1,'DisplayName','Reference Trajectory')
plot(t_ref,x,'b','Linewidth',1,'DisplayName','Trajectory')
plot(t_ref(locs),x_ref(1,locs),'r.','MarkerSize',12,'DisplayName','Peaks')
plot(t_ref(idx),x_ref(idx),'g.','MarkerSize',12,'DisplayName','Intervals')
xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',18)
xlim([t_ref(1) t_ref(end)])
legend
end
N_err = length(idx)-1;
t_err = t_ref(locs);
err_amp = zeros(size(x,1),N_err);
err_amp_rel = zeros(size(x,1),N_err);
err_phase = zeros(size(x,1),N_err);
err_phase_rel = zeros(size(x,1),N_err);
% figure(214); clf; hold on; grid on; box on;
% plot(t_ref,x_ref,'k','Linewidth',1,'DisplayName','Reference Trajectory')
% plot(t_ref,x,'Linewidth',1,'DisplayName','Trajectory')
for ii = 1:N_err
    x_ref_i = x_ref(idx(ii):idx(ii+1)); t_ref_i = t_ref(idx(ii):idx(ii+1));
    x_i = x(:,idx(ii):idx(ii+1));
    [vals,poss] = max(abs(x_i),[],2);
    inds = sub2ind(size(x_i),transpose(1:size(x,1)),poss);
    [val_ref,pos_ref] = max(abs(x_ref_i));
%     plot(t_ref_i(pos_ref),x_ref_i(pos_ref),'r.','MarkerSize',12,'DisplayName','Peaks Ref.')
%     plot(t_ref_i(poss),x_i(inds),'c.','MarkerSize',12,'DisplayName','Peaks Ref.')
%     drawnow;
    err_amp(:,ii) = vals.*sign(x_i(inds)) - val_ref*sign(x_ref_i(pos_ref));
    err_amp_rel(:,ii) = err_amp(:,ii) /val_ref*100;
    err_phase(:,ii) = transpose(t_ref_i(poss)-t_ref_i(pos_ref));
    err_phase_rel(:,ii) = err_phase(:,ii)/(t_ref_i(end)-t_ref_i(1))*100;
end
Errors = struct('t',t_err,'AmpError',err_amp,'AmpErrorR',err_amp_rel,...
                             'PhError',err_phase,'PhErrorR',err_phase_rel);
end
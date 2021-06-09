function surfV(V, embedDim, dim, varargin)

ztext = 'V';
if ~isempty(varargin)
    ztext = varargin{1};
end
    
Vplot = reshape(V(:,dim), [], embedDim);
if size(Vplot, 1) > 1
    figure
    surf(Vplot);
    shading interp
    xlabel('time-dim', 'Interpreter', 'latex')
    ylabel('$x$-dim', 'Interpreter', 'latex')
    zlabel(['$',ztext,'^{',num2str(dim),'}_{t,x}$'], 'Interpreter', 'latex')
    set(gca,'fontname','times')
    set(gca,'fontsize',18)
    
    figure
    plot(1:length(Vplot(:,1)), Vplot(:,1)/norm(Vplot(:,1)), 'LineWidth', 2)
    xlim([1, length(Vplot(:,1))])
    xlabel('$x$-dim', 'Interpreter', 'latex')
    ylabel(['$',ztext,'^{',num2str(dim),'}_{',num2str(1),',x}$'], 'Interpreter', 'latex')
    set(gca,'fontname','times')
    set(gca,'fontsize',18)
end

if size(Vplot, 2) > 1
    figure
    plot(1:embedDim, Vplot(1,:)/norm(Vplot(1,:)), 'LineWidth', 2)
    xlim([1, embedDim])
    xlabel('time-dim', 'Interpreter', 'latex')
    ylabel(['$',ztext,'^{',num2str(dim),'}_{t,',num2str(1),'}$'], 'Interpreter', 'latex')
    set(gca,'fontname','times')
    set(gca,'fontsize',18)
end
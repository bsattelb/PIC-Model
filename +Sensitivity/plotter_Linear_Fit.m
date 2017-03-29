function plotter_Linear_Fit(Nparams, Nsamples, paramNames, QOI, w, output, Xs, varargin)

    if ~exist('Results_LF', 'dir')
        mkdir('Results_LF');
    end
    fig = figure;
    plot(1:Nparams,w,'.-b','MarkerSize',30);
    title(['Global Linear Fit (N = ' int2str(Nsamples) ')'],'Interpreter','latex','Fontsize',16,'FontWeight','bold')
    xlabel('Parameter','Interpreter','latex','Fontsize',14)
    ylabel('Weight','Interpreter','latex','Fontsize',14)
    ylim([-1,1])
    xlim([0,Nparams+1])
    ax = gca;
    ax.XTick = 1:Nparams;
    ax.XTickLabel = paramNames;
    saveas(fig, ['Results_LF/WV.png'])
    clear fig;
    
    if isempty(varargin)
        fig = figure;
        plot(Xs*w, output,'ko');
        p = polyfit(Xs*w, output, 1);
        hold on;
        g = linspace(-2, 2, 100);
        plot(g, g*p(1) + p(2));
        hold off;
        title(['Sufficient Summary Plot (N = ' int2str(Nsamples) ')'],'Interpreter','latex','Fontsize',16,'FontWeight','bold')
        xlabel(['$w^T x_j$'],'Interpreter','latex','FontSize',14)
        ylabel(QOI,'Interpreter','latex','FontSize',14)
        xlim([-2,2])
        axis square;
        grid on;
        saveas(fig, ['Results_LF/SSP.png'])
        clear fig;
    else
        fig = figure;
        errorbar(Xs*w, output, 2*varargin{1}, 'ko');
        p = polyfit(Xs*w, output, 1);
        hold on;
        g = linspace(-2, 2, 100);
        plot(g, g*p(1) + p(2));
        hold off;
        title(['Sufficient Summary Plot (N = ' int2str(Nsamples) ')'],'Interpreter','latex','Fontsize',16,'FontWeight','bold')
        xlabel(['$w^T x_j$'],'Interpreter','latex','FontSize',14)
        ylabel(QOI,'Interpreter','latex','FontSize',14)
        xlim([-2,2])
        axis square;
        grid on;
        saveas(fig, ['Results_LF/SSP.png'])
        clear fig;
    end
    
end
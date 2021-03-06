function plotter_Local_Linear_Model(Nparams, Nsamples, paramNames, QOI, evalues, U, output, Xs, varargin)

    if ~exist('Results_LLRM', 'dir')
        mkdir('Results_LLRM');
    end

    fig = figure;
    semilogy(1:Nparams, evalues, '.-b', 'MarkerSize', 30)
    title(['Approximate Eigenvalues, (N = ' int2str(Nsamples) ')'], 'Interpreter', 'latex','Fontsize', 16,'FontWeight','bold')
    xlim([0,Nparams+1])
    saveas(fig, 'Results_LLRM/evalues.png')
    clear fig;

    for ii = 1:Nparams
        fig = figure;
        plot(1:Nparams, U(:,ii), '.-b', 'MarkerSize',30)
        title(['Weight Vector ' int2str(ii), '(N = ' int2str(Nsamples) ')'],'Interpreter','latex','Fontsize',16,'FontWeight','bold')
        xlabel('Parameters','Interpreter','latex','FontSize',14)
        ylabel('Parameter Weights','Interpreter','latex','FontSize',14)
        xlim([0,Nparams+1])
        ylim([-1,1])
        ax = gca;
        ax.XTick = 1:Nparams;
        ax.XTickLabel = paramNames;
        saveas(fig, ['Results_LLRM/WV' int2str(ii) '.png'])
        clear fig;

		fig = figure;
        if isempty(varargin)
            plot(Xs*U(:,ii), output,'ko');
        else
            errorbar(Xs*U(:,ii), output, 2*varargin{1}, 'ko');
        end
		title(['Sufficient Summary Plot ' int2str(ii) ' (N = ' int2str(Nsamples) ')'],'Interpreter','latex','Fontsize',16,'FontWeight','bold')
		xlabel(['$w^T_' int2str(ii) ' x_j$'],'Interpreter','latex','FontSize',14)
		ylabel(QOI,'Interpreter','latex','FontSize',14)
		xlim([-2,2])
		axis square;
		grid on;
		saveas(fig, ['Results_LLRM/SSP' int2str(ii) '.png'])
		clear fig;
    end

    if Nparams >= 3
        fig = figure;
        scatter3(Xs*U(:,1), Xs*U(:, 2), Xs*U(:, 3), 16, output, 'filled')
        title(['Heat Map (N = ' int2str(Nsamples) ')'],'Interpreter','latex','Fontsize',16,'FontWeight','bold')
        xlabel('$w^T_1 x_j$','Interpreter','latex','FontSize',14)
        ylabel('$w^T_2 x_j$','Interpreter','latex','FontSize',14)
        zlabel('$w^T_3 x_j$','Interpreter','latex','FontSize',14)
        xlim([-2,2])
        ylim([-2,2])
        zlim([-2,2])
        grid on;
        colorbar
        saveas(fig, 'Results_LLRM/HeatMap.png')
        saveas(fig, 'Results_LLRM/HeatMap.fig')
        clear fig;
    end
end
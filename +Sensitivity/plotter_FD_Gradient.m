function plotter_Active_Subspaces(Nparams, Nsamples, paramNames, QOI, evalues, U, output, Xs)

    if ~exist('Results_FDG', 'dir')
        mkdir('Results_FDG');
    end

    fig = figure;
    semilogy(1:Nparams, evalues, '.-b', 'MarkerSize', 30)
    title(['Approximate Eigenvalues, (N = ' int2str(Nsamples) ')'], 'Interpreter', 'latex','Fontsize', 16,'FontWeight','bold')
    xlim([0,Nparams+1])
    saveas(fig, 'Results_FDG/evalues.png')
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
        saveas(fig, ['Results_FDG/WV' int2str(ii) '.png'])
        clear fig;

        fig = figure;
        plot(Xs*U(:,ii), output,'ko');
        title(['Sufficient Summary Plot ' int2str(ii) ' (N = ' int2str(Nsamples) ')'],'Interpreter','latex','Fontsize',16,'FontWeight','bold')
        xlabel(['$w^T_' int2str(ii) ' x_j$'],'Interpreter','latex','FontSize',14)
        ylabel(QOI,'Interpreter','latex','FontSize',14)
        xlim([-2,2])
        axis square;
        grid on;
        saveas(fig, ['Results_FDG/SSP' int2str(ii) '.png'])
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
        saveas(fig, 'Results_FDG/HeatMap.png')
        saveas(fig, 'Results_FDG/HeatMap.fig')
        clear fig;
    end
end
%% Make Plots of the Bounds (value space)
figure
    ii=1;
    whichsquares=[ (n_rectx+1)/2-1  (n_rectx+1)/2 (n_rectx+1)/2+1];
    whichsquares=(n_rectx+1)/2;
    for x = whichsquares
    subplot(1,length(whichsquares),ii)
    
    plot((1:size(Brv,1))/4,Brv(:,x),'.','MarkerSize',10,'color', 'b','HandleVisibility','off');%'linewidth',2,'color', [t_cost/(max(t_costs)+t_cost) t_cost/(max(t_costs)+t_cost) t_cost/(max(t_costs)+t_cost)],'HandleVisibility','off');
    hold on
    %plot((1:size(Blv,1))/4,Blv(:,x),'.','MarkerSize',10,'color', 'r', 'DisplayName', ['y' num2str(1/(t_cost*n_recty))]);%'linewidth',2,'color', [t_cost/(max(t_costs)+t_cost) t_cost/(max(t_costs)+t_cost) t_cost/(max(t_costs)+t_cost)], 'DisplayName', ['y' num2str(1/(t_cost*n_recty))]);
    %plot((1:size(Blv,1))/4,.5*ones(size(Blv,1)),'k', 'HandleVisibility','off');
    %plot((1:size(Blrv,1))/4,Blrv(:,x),'.','MarkerSize',10,'color', 'g', 'DisplayName', ['y' num2str(1/(t_cost*n_recty))]);%'linewidth',2,'color', [t_cost/(max(t_costs)+t_cost) t_cost/(max(t_costs)+t_cost) t_cost/(max(t_costs)+t_cost)], 'DisplayName', ['y' num2str(1/(t_cost*n_recty))]);
    if ii>1
        set(gca,'YTickLabel',[])
    else
        ylabel('Belief (g)')
    end
    set(gca,'fontsize',18)
    set(gca,'linewidth',2)
    box off
    %ylim([0 1])
    xlabel('Time (s)')
    %title(['x = ' ' ' num2str(x-(n_rectx+1)/2)])
    legend('-DynamicLegend', 'Location', 'EastOutside');   
    
    ii=ii+1;
    end
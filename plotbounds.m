%%
figure
    ii=1;
    whichsquares=[ 21 22 23];
    %whichsquares=(n_rectx+1)/2;
    for x_ = whichsquares
    subplot(1,length(whichsquares),ii)
    
    plot((1:size(Br,1))/4,Br(:,x_),'.','color', 'k','HandleVisibility','off');%'linewidth',2,'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)],'HandleVisibility','off');
    hold on
    plot((1:size(Bl,1))/4,Bl(:,x_),'.','color', 'k', 'DisplayName', ['x' num2str(n_rectx)]);%'linewidth',2,'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)], 'DisplayName', ['x' num2str(n_rectx)]);
    plot((1:size(Bl,1))/4,.5*ones(size(Bl,1)),'k', 'HandleVisibility','off');
    if ii>1
        set(gca,'YTickLabel',[])
    else
        ylabel('Belief (g)')
    end
    set(gca,'fontsize',18)
    set(gca,'linewidth',2)
    box off
    ylim([0 1])
    xlabel('Time (s)')
    title(['x = ' ' ' num2str(x_-(n_rectx+1)/2)])
    legend('-DynamicLegend', 'Location', 'EastOutside');   
    
    ii=ii+1;
    end
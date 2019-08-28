%%
mainfig=figure;
set(mainfig,'position',[0 100 1200 500]);
set(mainfig,'color','white')  
difficulty = nanmean(stim,2);
[sorted idx] = sort(difficulty);

  subplot(1,2,1)
  hold on
  whichstim= idx(end-9)';
for i = whichstim
    yyaxis right

    %plot(.25*[1:last(i)],1.02*(stim2(i,1:last(i))+max(stim2(:)))/2, 'k.', 'markersize', 20)

    plot(.25*[1:last(i)],stimuli(i,1:last(i)), 'r-', 'markersize', 20)
    ylabel('Belief')
    ylim([-0.05 1.05])

    yyaxis left

if choice(i)==1
 plot(.25*[1:last(i)],x(i,1:last(i))-x_num-1,'r','linewidth',2)
elseif choice(i)==n_rectx
      plot(.25*[1:last(i)],x(i,1:last(i))-x_num-1,'b','linewidth',2)
end
xlim([0 7])
ylim([-x_num x_num])

xlabel('Time (s)')
ylabel('x')
title('Easy')
set(gca,'fontsize',18)
set(gca,'linewidth',2)
box off
end

subplot(1,2,2)
%% Plot trajectories
figure
difficulty = nanmean(stim2,2);
[sorted idx] = sort(difficulty);
%figure
whichstim=idx(1:6)';%75
whichstim=1001;
C = {[.5 .6 .7],'b','k','g','y','r',[.8 .2 .6]};
for i = whichstim
    hold on
    yyaxis right
    plot(.25*[1:last(i)],1.02*(stim2(i,1:last(i))+max(stim2(:)))/2, 'k.', 'markersize', 20)%'k.', 'markersize', 20)
    hold on
    %plot(.25*[1:length(stim(i,:))],stimuli(i,:), '--r', 'linewidth', 2)%plot(.25*[1:last(i)],stimuli(i,1:last(i)), '--r', 'linewidth', 2)%'r-', 'markersize', 20)
    ylabel('Belief')
    ylim([-0.05 1.05])
    yyaxis left
    hold on
if choice(i)==1
plot(.25*[1:last(i)],x(i,1:last(i))-x_num-1,'-','color','r','linewidth',2)%C{i},'linewidth',2)%'r','linewidth',2)
elseif choice(i)==n_rectx
plot(.25*[1:last(i)],x(i,1:last(i))-x_num-1,'-','color',C{i},'linewidth',2)%'r','linewidth',2)%plot(.25*[1:last(i)],x(i,1:last(i))-x_num-1,'-b','linewidth',2)
end

xlim([0 n_recty])
ylim([-x_num x_num])

xlabel('Time (s)')
ylabel('x')
%title(['(NR-NL)/(NR+NL): ' num2str(nanmean(stim(whichstim,:)))])
%title(['(NR-NL)/(NR+NL): ' num2str(nanmean(stim(whichstim(2),:)))])
set(gca,'fontsize',18)
set(gca,'linewidth',2)
box off
%legend('Trajectory', 'Stimulus','Belief', 'Location', 'southeast')
end 

%% Plot belief
figure
difficulty = nanmean(stim2,2);
[sorted idx] = sort(difficulty);
%figure
whichstim=idx(1:6)';%75
%whichstim=9;41
C = {[.5 .6 .7],'b','k','g','y','r',[.8 .2 .6]};
for i = whichstim
    hold on
    yyaxis right
    plot(.25*[1:length(stim(i,:))],1.02*(stim(i,:)+max(stim(:)))/2, 'k.', 'markersize', 20)%plot(.25*[1:last(i)],1.02*(stim2(i,1:last(i))+max(stim2(:)))/2, 'k.', 'markersize', 20)%'k.', 'markersize', 20)
    hold on
    %CM = jet(i);
    plot(.25*[1:length(stim(i,:))],stimuli(i,:), '--','color',C{i}, 'linewidth', 2)%plot(.25*[1:last(i)],stimuli(i,1:last(i)), '--r', 'linewidth', 2)%'r-', 'markersize', 20)
    ylabel('Belief')
    ylim([-0.05 1.05])
    yyaxis left
    hold on
if choice(i)==1
%plot(.25*[1:last(i)],x(i,1:last(i))-x_num-1,'r','linewidth',2)%'r','linewidth',2)
elseif choice(i)==n_rectx
%plot(.25*[1:last(i)],x(i,1:last(i))-x_num-1,'b','linewidth',2)
end

xlim([0 n_recty])
ylim([-1 1])%ylim([-x_num x_num])

xlabel('Time (s)')
ylabel('Stimulus')%ylabel('x')
%title(['(NR-NL)/(NR+NL): ' num2str(nanmean(stim(whichstim,:)))])
%title(['(NR-NL)/(NR+NL): ' num2str(nanmean(stim(whichstim(2),:)))])
set(gca,'fontsize',18)
set(gca,'linewidth',2)
box off
%legend('Trajectory', 'Stimulus','Belief', 'Location', 'southeast')
end 
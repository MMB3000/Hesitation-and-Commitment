%% Plot Value function
figure
hold on
t_step=4;
plot(gn{t_step},Value{t_step}(:,1),'.r','MarkerSize',20,'DisplayName','Move Left')
plot(gn{t_step},Value{t_step}(:,2),'.k','MarkerSize',20,'DisplayName','Wait')
plot(gn{t_step},Value{t_step}(:,3),'.b','MarkerSize',20,'DisplayName','Move Right')
xlabel('Belief (g)')
ylabel('Value')
title(['Value at t = ' num2str(t_step/4.) ' s'])
legend('-DynamicLegend', 'Location', 'East');   
%%
figure
hold on
times = linspace(1,56,56);
for t_step=1:56
    plot3(t_step/4.*ones(size(gn{t_step})), gn{t_step},Value{t_step}(:,1),'.r')
    plot3(t_step/4.*ones(size(gn{t_step})), gn{t_step},Value{t_step}(:,2),'.k')
    plot3(t_step/4.*ones(size(gn{t_step})), gn{t_step},Value{t_step}(:,3),'.b')
    
    plot3(t_step/4., Bl(t_step,2), Value{t_step}(idivide(int16(t_step),int16(2))+1,1),'.r','MarkerSize',25)
    plot3(t_step/4., Blr(t_step,2), Value{t_step}(idivide(int16(t_step),int16(2))+1,1),'.k','MarkerSize',25)
    plot3(t_step/4., Br(t_step,2), Value{t_step}(idivide(int16(t_step),int16(2))+1,3),'.b','MarkerSize',25)
end
zlim([.1 .9])
view(90,0)
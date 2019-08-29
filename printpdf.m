%%Print pdf of current figure
title_file = 'plots/stimulus_belief_easy.pdf';

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,title_file,'-dpdf','-r0')
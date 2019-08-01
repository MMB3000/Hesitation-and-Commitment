%% Plot 2D Histogram of RT vs Stimulus strength
figure
% ylim([2 14])
% xlim([0 0.7])
[N,Xbins,Ybins] = hist2d(difficulty,last/4,'tile');
hold on
colorbar()
set(gca,'YDir','normal')
xlabel('(NR-NL)/(NR+NL)')
ylabel('Reaction Time (s)')
h = colorbar();
ylabel(h, 'Counts');set(gca,'YDir','normal')
title(['Time Limit = ' ' ' num2str(1/(t_cost*n_recty))])
hold on

%% Plot 2D Histogram of RT vs Absolute Stimulus strength
figure
% ylim([2 14])
% xlim([0 0.7])
[N,Xbins,Ybins] = hist2d(abs(difficulty),last/4,'tile');
hold on
Xbins;
Ybins;
N = flip(N);
matriz = bsxfun(@times,flip(Ybins)',N);
vector1 = sum(matriz,1);
vector2 = sum(N,1);
means = bsxfun(@rdivide,vector1,vector2);
plot(Xbins,means,'.k','MarkerSize',20)
h = colorbar();
ylabel(h, 'Counts');set(gca,'YDir','normal')
xlabel('|(NR-NL)/(NR+NL)|')
ylabel('Reaction Time (s)')
title(['Time Limit = ' ' ' num2str(1/(t_cost*n_recty))])
hold on
%% Plot conditional probability plot
figure
[N,Xbins,Ybins] = hist2d(abs(difficulty),last/4,'tile');
figure
%N = flip(N);
matriz = bsxfun(@times,Ybins',N);%flip(Ybins)',N)
vector1 = sum(matriz,1);
vector2 = sum(N,1);
cond_prob = bsxfun(@rdivide,N,vector2);
means = bsxfun(@rdivide,vector1,vector2);
set(gca,'YDir','normal')
xlabel('|(NR-NL)/(NR+NL)|')
ylabel('Reaction Time (s)')
title(['Time Limit = ' ' ' num2str(1/(t_cost*n_recty))])
imagesc(Xbins,Ybins,cond_prob)
hold on
plot(Xbins,means,'.k','MarkerSize',20)
h = colorbar();
ylabel(h, 'Conditional probability');
set(gca,'YDir','normal')
xlabel('|(NR-NL)/(NR+NL)|')
ylabel('Reaction Time (s)')
title(['Time Limit = ' ' ' num2str(1/(t_cost*n_recty))])
hold on
st=regstats(last/4,abs(difficulty));
if st.tstat.pval(2)<.05
h=    plot(abs(difficulty),st.beta(1)+st.beta(2)*abs(difficulty),'r','linewidth',3);

elseif st.tstat.pval(2)>.05
h=         plot(abs(difficulty),st.beta(1)+st.beta(2)*abs(difficulty),'b','linewidth',3);
end
  legend(h,['\beta = ' ' '  num2str(round(st.tstat.beta(2)*100)/100),'p=' num2str(round(st.tstat.pval(2)*10000)/10000)]);

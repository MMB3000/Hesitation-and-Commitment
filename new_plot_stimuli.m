%Fix number of image (move to end)
%set same ylim for last figures
%set bold font and lines for last figures
%change namesave to y size not tc

maintic = tic;

n_recty = 14;                 %maximum time allowed (must be whole number)
g_num=500;


%% Generate Stimuli
  nstim = 1;
  tic
  stim = GenerateStimulus( nstim, n_recty );    
  toc

%% Transform to Belief with noise
alpha =0.2;

mainfig1=figure;
set(mainfig1,'position',[0 100 800 500]);
set(mainfig1,'color','white')

for repeats = 1:5
  
  tvec = 0:.25:n_recty;
  stim_normalized = stim/max(stim(:)); %normalize so that right=1 left=-1
  
  %flips
  flips = rand(1,size(stim_normalized,2));
  flips(flips>alpha) = 1;
  flips(flips<=alpha) = -1;
  stim_normalized = bsxfun(@times, stim_normalized, flips);
  %end of flips
  
  runs= nan(nstim,length(tvec)-1);
  p= nan(nstim,length(tvec)-1);
  stimuli= nan(nstim,length(tvec)-1);
  tic
  
    for i = 1:size(stim_normalized,1)
            ns=cumsum(ones(size(stim_normalized(i,:))));
            tempvar=[];
            for n = ns
                
                r   = linspace(0,n,n+1);
                g_numerator = bsxfun(@minus,betainc(1-alpha,r+1,n-r+1),betainc(0.5,r+1,n-r+1));
                g_denominator = bsxfun(@minus,betainc(1-alpha,r+1,n-r+1),betainc(alpha,r+1,n-r+1));
                g = bsxfun(@rdivide,g_numerator,g_denominator);
                
               
                r_n = linspace(0,n+1,n+2);
                g_n_numerator = bsxfun(@minus,betainc(1-alpha,r_n+1,n-r_n+2),betainc(0.5,r_n+1,n-r_n+2));
                g_n_denominator = bsxfun(@minus,betainc(1-alpha,r_n+1,n-r_n+2),betainc(alpha,r_n+1,n-r_n+2));
                g_n = bsxfun(@rdivide,g_n_numerator,g_n_denominator);
                
                
                
                rs=sum(max(stim_normalized(i,1:n),0)); %number of rigth counts from stimulus
                
                random_next = rand;
                randomnext(i,n) = random_next;
                if random_next>g(rs+1)
                    p(i,n) = g_n(rs+2);
                else
                    p(i,n) = g_n(rs+1);
                end
                
                stimuli(i,n)=p(i,n);
                
            end
            
    end
    toc
    
   
    figure(1)
    if repeats ==1

    yyaxis left
    hold on
    plot(.25*(1:length(stim(1,:))), stim(1,:)/max(stim(:)), '.k', 'markersize', 20, 'DisplayName', 'Stimulus')
    ylim([-1.02 1.02])
    ylabel('Stimulus')
    yyaxis right
    hold on
    plot(.25*(1:length(p(1,:))), p(1,:), '--', 'color', [0.7 .7 .7], 'linewidth', 3, 'DisplayName', strcat('Belief \alpha =',num2str(alpha)))
    else
    plot(.25*(1:length(p(1,:))), p(1,:), '--', 'color', [0.7 .7 .7], 'linewidth', 3 ,'HandleVisibility','off')
    end

    ylim([-0.01 1.01])
end
    
%% Transform to Belief no noise
alpha =0.;
for repeats = 1%:5
  stim_normalized = stim/max(stim(:));
  runs= nan(nstim,length(tvec)-1);
  %runsBin = nan(nstim,length(tvec)-1);
  p= nan(nstim,length(tvec)-1);
  stimuli= nan(nstim,length(tvec)-1);
  tic
  
    for i = 1:size(stim_normalized,1)
            ns=cumsum(ones(size(stim_normalized(i,:))));
            tempvar=[];
            for n = ns
                
                r   = linspace(0,n,n+1);
                g_numerator = bsxfun(@minus,betainc(1-alpha,r+1,n-r+1),betainc(0.5,r+1,n-r+1));
                g_denominator = bsxfun(@minus,betainc(1-alpha,r+1,n-r+1),betainc(alpha,r+1,n-r+1));
                g = bsxfun(@rdivide,g_numerator,g_denominator);
                
                
                r_n = linspace(0,n+1,n+2);
                g_n_numerator = bsxfun(@minus,betainc(1-alpha,r_n+1,n-r_n+2),betainc(0.5,r_n+1,n-r_n+2));
                g_n_denominator = bsxfun(@minus,betainc(1-alpha,r_n+1,n-r_n+2),betainc(alpha,r_n+1,n-r_n+2));
                g_n = bsxfun(@rdivide,g_n_numerator,g_n_denominator);
                
                
                
                rs=sum(max(stim_normalized(i,1:n),0));
%                 rsn=.5+rs+(randn(size(rs))*(sigma*n));
%                 tempvar=knnsearch(r_n',rsn');
%                 runs(i,n)=tempvar;
                %runsBin(i,n)=runs(i,n);
                
                random_next = rand;
                randomnext(i,n) = random_next;
                if random_next>g(rs+1)
                    p(i,n) = g_n(rs+2);
                else
                    p(i,n) = g_n(rs+1);
                end
                
                stimuli(i,n)=p(i,n);
                
            end
            
    end
    toc    
    
%%
figure(1)
yyaxis right
hold on
if repeats ==1
    plot(.25*(1:length(p(1,:))), p(1,:), '-', 'color', [0 .4 0], 'linewidth', 3, 'DisplayName', 'Belief \alpha = 0')
else
    plot(.25*(1:length(p(1,:))), p(1,:), '-', 'color', [0 .4 0], 'linewidth', 3 ,'HandleVisibility','off')
end
ylabel('Belief')
xlim([0 14])
legend('-DynamicLegend', 'Location', 'best')
end
title('Stimulus vs. Agent Belief')
xlabel('Time (s)')


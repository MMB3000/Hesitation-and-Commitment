%Fix number of image (move to end)
%set same ylim for last figures
%set bold font and lines for last figures
%change namesave to y size not tc

%clear all
clearvars -except easystim hardstim
close all

rng(10)                       %set seed for reproducibility of results
maintic = tic;

n_recty = 14;                 %maximum time allowed (must be whole number)
g_num=500;
alpha = .2;
x_cost = .01;
end_penalty = -1;

%% Generate Stimuli
nstim = 5000;
tic
stim = GenerateStimulus( nstim, n_recty );
%stim = [easystim;easystim;easystim;easystim;easystim;easystim];
%stim = [hardstim;hardstim;hardstim;hardstim;hardstim;hardstim];
toc
stim = stim/max(stim(:));  
%% Transform to Belief
tvec = 0:.25:n_recty;
% 8
runs= nan(nstim,length(tvec)-1);
stimuli= nan(nstim,length(tvec)-1);
tic
rs  = cumsum(max(stim,0),2);
ns  = size(stim,2);
% r_n = nan(ns,g_num+2);
% g_n = nan(ns,g_num+2); 
% for n = 1:ns
%     rnn      = linspace(0,n,g_num+2);
%     r_n(n,:) = rnn;
%     g_n(n,:) = betainc(.5,rnn+1,n+1-rnn,'upper');
% end
% r_n = r_n';
% 
% for n = 1:ns
%     
%     rsn     = rs(:,n) + randn(size(rs(:,n)))*sigma.*sqrt(rs(:,n).*(n-rs(:,n)));
%     tempvar = knnsearch(r_n(:,n),rsn);
%     runs(:,n)   = tempvar;
%     stimuli(:,n)= g_n(n,tempvar);
%     
% end

nmax = size(stim,2);

%flips
flips = rand(nstim,size(stim,2));
flips(flips>alpha) = 1;
flips(flips<=alpha) = -1;
stim_noisy = bsxfun(@times, stim, flips);
%end of flips

for i = 1:size(stim,1)
    %ns=cumsum(ones(size(stim(i,:))));
    tempvar=[];
    for n = 1:ns

        r   = linspace(0,n,n+1);
        g_numerator = bsxfun(@minus,betainc(1-alpha,r+1,n-r+1),betainc(0.5,r+1,n-r+1));
        g_denominator = bsxfun(@minus,betainc(1-alpha,r+1,n-r+1),betainc(alpha,r+1,n-r+1));
        g = bsxfun(@rdivide,g_numerator,g_denominator);


        r_n = linspace(0,n+1,n+2);
        g_n_numerator = bsxfun(@minus,betainc(1-alpha,r_n+1,n-r_n+2),betainc(0.5,r_n+1,n-r_n+2));
        g_n_denominator = bsxfun(@minus,betainc(1-alpha,r_n+1,n-r_n+2),betainc(alpha,r_n+1,n-r_n+2));
        g_n = bsxfun(@rdivide,g_n_numerator,g_n_denominator);



        rs=sum(max(stim_noisy(i,1:n),0));
%                 rsn=.5+rs+(randn(size(rs))*(sigma*n));
%                 tempvar=knnsearch(r_n',rsn');
%                 runs(i,n)=tempvar;
        %runsBin(i,n)=runs(i,n);

        %p(i,n) = g_n(rs+1);
        
        %using belief for next timestep
%         random_next = rand;
%         randomnext(i,n) = random_next;
%         if random_next>g(rs+1)
%             p(i,n) = g_n(rs+2);
%         else
%             p(i,n) = g_n(rs+1);
%         end
%         
%         stimuli(i,n)=p(i,n);
        %
        
        stimuli(i,n)=g(rs+1);

    end
end

toc

%% Loop over time costs
%t_costs = [1/(64*n_recty) 1/(32*n_recty) 1/(16*n_recty) 1/(8*n_recty) 1/(4*n_recty) 1/(2*n_recty) 1/(1*n_recty)];

%n_rectxs = [ 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 ]; %:2:41];% 41 51 55];%21 31];
n_rectxs = 17;%[ 3 5 7 9 11 21 31 41 51 61 71 81 91 101 111];
stim = stim/max(stim(:));
%t_costs = 1/(64*n_recty);
% indx=length(t_costs)+1:length(t_costs)+n_recty;
% for i = 1:n_recty-1  
%t_costs(indx(i))= 1/(4*(n_recty-i));
% end

t_costs = 1/(4*n_recty);
%%

all_pcs1        = nan(length(t_costs),length(n_rectxs));
all_pcs2        = nan(length(t_costs),length(n_rectxs));
all_npresses    = nan(length(t_costs),length(n_rectxs));
all_starts      = nan(length(t_costs),length(n_rectxs));
all_betas       = nan(length(t_costs),length(n_rectxs));
all_poststarts  = nan(length(t_costs),length(n_rectxs));

%pdf printing preferences
matlab.graphics.internal.setPrintPreferences('DefaultPaperPositionMode','auto')


for dim1=1:length(t_costs);
    
    t_cost = t_costs(dim1);
%n_rectxs=17;
tout=1;
%uncomment:
mainfig1=figure;
set(mainfig1,'position',[0 100 800 500]);
set(mainfig1,'color','white')
 
mainfig1=figure;
set(mainfig1,'position',[0 100 800 500]);
set(mainfig1,'color','white')
 
mainfig1=figure;
set(mainfig1,'position',[0 100 800 500]);
set(mainfig1,'color','white')
 
mainfig1=figure;
set(mainfig1,'position',[0 100 1200 500]);
set(mainfig1,'color','white')
 
fig=figure;set(fig,'color','white')
set(fig,'position',[0 100 1200 400]);
 
mainfig1=figure;
set(mainfig1,'position',[0 100 800 500]);
set(mainfig1,'color','white')

mainfig1=figure;
set(mainfig1,'position',[0 100 800 500]);
set(mainfig1,'color','white')
 
mainfig1=figure;
set(mainfig1,'position',[0 100 1000 500]);
set(mainfig1,'color','white')
 
mainfig1=figure;
set(mainfig1,'position',[0 100 800 500]);
set(mainfig1,'color','white')
%end of stuff to uncomment

for dim2 = 1:length(n_rectxs)
    
    n_rectx = n_rectxs(dim2); 
    
clearvars -except nstim maintic dim1 dim2 all_betas all_poststarts all_pcs1 all_pcs2 all_starts all_npresses x_cost x_costs alpha sigma stim stimuli stim_noisy p runs n_recty n_rectx n_rectxs t_cost t_costs pc tout easystim hardstim

%function [] = BinomialAgent() 
%% SET GRID PARAMETERS
%n_rectx = 5;                 %number of button presses to between sides
%n_recty = 14;                 %maximum time allowed (must be whole number)
%Rw=100;                      %Reward for correct responses
%Rp=25;                       %Penalty for incorrect responses
%t_cost = 1/(4*n_recty/4);

%% SET STIMULUS PARAMETERS

stim_duration_sec=n_recty;    %in seconds
ramp_stim_dur_msec=0;         %in miliseconds

chunk_duration_msec=250;      %in miliseconds (includes rise and fall time)
filter_win_msec=30;           %Cosine filtering window in miliseconds
inter_burst_time = 150;       %Stimulus processing time
refractoryT = 0;              %Refractory period on button presses


%% Initialize Value Function

dt = .25;
rng(1);

t_num = n_recty/(chunk_duration_msec/1000);

[ gn, gpn, dgn,dgpn, pggn ] = new_new_belief( alpha, t_num );

% Value = nan(t_num,g_num,n_rectx);
% stop = 0
% for n = 1:t_num
%     
%     Value(n,:,1)    = 1-gn(:,n)';
%     Value(n,:,end)  = gn(:,n)';
%     if stop==0
%         Value(n,:,1)
%         Value(n,:,end)
%         stop = 1;
%     end
%     
% end


Value = {};%nan(t_num,g_num,n_rectx);
stop = 0;
for t = 1:t_num
    Value{t} = nan(t+1,n_rectx); %t+1 because first timestep has 2 possibilities for belief
end

for t = 1:t_num %t:=timestep
    Value{t}(:,1) = 1-gn{t}';
    Value{t}(:,end) = gn{t}';
end

Value{t_num}(:,2:end-1) = end_penalty;  %should be bigger than accumulated time cost at t=t_num

% for n = 1:t_num
%     for n_cell = 1:t_num
%         if n<n_cell+1
%               Value{n_cell}(n,1) = 1-gn{n_cell}(n);
%               Value{n_cell}(n,n_rectx) = 1-gn{n_cell}(n);
%         end
%     end
% end


%Vm = NaN(t_num-1, g_num,n_rectx);  %index of the position of max(E[V]) %why t_num-1
Vm = {}; %check what this is exactly <<<<<<<<
for t = 1:t_num-1
    Vm{t} = nan(t+1,n_rectx); % check whether it should be t or t+1<<<<<<<<
end

xtgrid = ones(1,n_rectx);
xtgrid(1,[1, n_rectx])=[0,0];
TimeCosts = -t_cost*dt*xtgrid; 

 
%% Bellman Equation solution by back-propagation
tic
%back propagation in time (t)  
Values = {};
for t = 1:t_num %t:=timestep
    Values{t} = nan(t+1,n_rectx,n_rectx);
end


for t = 1:t_num
      Values{t}(:,1,1) = squeeze(Value{t}(:,1));
      Values{t}(:,end,end) = squeeze(Value{t}(:,end));
end


for t = t_num-1:-1:1
    % int{p(g(t+dt)|g(t))*V(t+dt,g(t+dt),x(t+dt))dg(t+dt)}
    gg  = pggn{t};
    g   = gn{t};
    gp  = gpn{t};
    %dg  = dgn{t}; %this should not be necessary/valid
    %dgp = dgpn{t}; %^
    
    Evidence = gg*Value{t+1}(:,:); % equation 4.3 thesis
    
    
    % over all non-terminal x (j)
    
    for j = 2:n_rectx-1
        MovementCosts = abs((1:n_rectx)-j);
        MovementCosts(MovementCosts>1)=Inf;
        MovementCosts = -x_cost*MovementCosts;
        ExpectedValues = bsxfun(@plus,Evidence,MovementCosts+TimeCosts);
        [Value{t}(:,j), Vm{t}(:,j)] = max(ExpectedValues+eps*randn(size(ExpectedValues)),[],2); %<<<<<<<<<<<<<<<<<<<<<<<<<<<<< I do not understand this completely
        %we dont need the values, only the location of the max, Vm
        Values{t}(:,j,:) = ExpectedValues;
    end
end
toc


%% Calculate bounds from the Intersections of action Values

Blr     = nan(t_num,n_rectx); %check what this is exactly <<<<<<<< boundary lr g-space?
Blrv    = nan(t_num,n_rectx); %^ boundary lr V-space?
Br      = nan(t_num,n_rectx); %^ 
Brv     = nan(t_num,n_rectx); %^
Bl      = nan(t_num,n_rectx); %^
Blv     = nan(t_num,n_rectx); %^

for t = 1:t_num
    index=1;
    for x = 2:n_rectx-1 %all non-terminal x
        
        vr = Values{t}(:,x,x+1)'; %value of going to the right
        vl = Values{t}(:,x,x-1)'; %value of going to the left
        vw = Values{t}(:,x,x)'; %value of waiting
        
        %comment below
%            subplot(1,3,index)
%            hold on
%            plot3(t/4*ones(size(gn{t})),gn{t},vr,'b')
%            hold on
%            plot3(t/4*ones(size(gn{t})),gn{t},vl,'r')
%            plot3(t/4*ones(size(gn{t})),gn{t},vw,'k')
        %comment above
        
        rlw = InterX([gn{t}'; vl],[gn{t}'; vw]); % intersection (g,V) of value left and value wait? (why r?) <<<<<<
        rwr = InterX([gn{t}'; vw],[gn{t}'; vr]); %^ ...
        rlr = InterX([gn{t}'; vl],[gn{t}'; vr]); %^ ...
                         
        if ~isempty(rlw) &&  ~isempty(rwr) && rlw(2)>rlr(2)
            %plot3(t/4,rlw(1),rlw(2),'.r','markersize',15);%comment
            
            Blv(t,x)    = rlw(2);
            Bl(t,x)     = rlw(1);
            
        end
        
        if ~isempty(rwr) && ~isempty(rlw) && rwr(2)>rlr(2)
            %plot3(t/4,rwr(1),rwr(2),'.b','markersize',15);%comment
            Brv(t,x)    = rwr(2);
            Br(t,x)     = rwr(1);
        end
        if~isempty(rlr)
            %plot3(t/4,rlr(1),rlr(2),'.g','markersize',15);%comment
            Blrv(t,x)   = rlr(2);
            Blr(t,x)    = rlr(1);
        end
        
        tbar=t_num-1-(n_rectx-1)/2;
        if (n_rectx-1)/2 >= t_num-1
            tbar=2;
        end
        if t>tbar-1 && isempty(rlw) && ~isempty(rlr)
            Blv(t,x)    = rlr(2);
            Bl(t,x)     = rlr(1);
        end
        if t>tbar-1 && isempty(rwr)  && ~isempty(rlr)
            Brv(t,x)    = rlr(2);
            Br(t,x)     = rlr(1);
        end
        
%         if Bl(t,x)==0
%             Bl(t,x)=nan;
%         end
%         
%         if Br(t,x)==1
%             Bl(t,x)=nan;
%         end
        
        lbalt           = Bl(1:t,x);
        lbalt(lbalt==0) = [];
        if ~any(~isnan(lbalt)) && isempty(rlw)
            %Blv(t,x)=rlr(2); %comment
            Bl(t,x)     = 0;%rlr(1)
        end
        rbalt           = Br(1:t,x);
        rbalt(rbalt==1) = [];
        if ~any(~isnan(rbalt)) && isempty(rwr)
            %Brv(t,x)=rlr(2);%comment
            Br(t,x)     = 1;%rlr(1);
        end
               
%          zlim([.1 .9])
%          title(['x =' ' ' num2str(x-(n_rectx-1)/2-1) ])
%          xlabel('t')
%          ylabel('g')
%          if index==1
%            zlabel('Value')
%          end
%          index=index+1;
%          hold off
%          view(80,74)
    end
end

for i=1:size(Bl,2)
    if nanmean(Bl(:,i))==1 || nanmean(Bl(:,i))==0
       Bl(:,i)=nan;
    end
    if nanmean(Br(:,i))==1 || nanmean(Br(:,i))==0
       Br(:,i)=nan;
    end    
end
      

% %%
% ts=[];
% xs=[];
% bls=[];
% brs=[];
% bws=[];
% counter=1;
% for i = 1:size(Bl,2)
% bls = [bls; Bl(:,i)];
% brs = [brs; Br(:,i)];
% bws = [bws; Blr(:,i)];
% 
% ts = [ts 1:size(Bl,1)];
% xs = [xs counter*ones(1,size(Bl,1))];
% counter = counter+1;
% end

%% Make Plots of the Bounds
figure(4)
    ii=1;
    whichsquares=[ (n_rectx+1)/2-1  (n_rectx+1)/2 (n_rectx+1)/2+1];
    for x = whichsquares
    subplot(1,length(whichsquares),ii)
    
    plot((1:size(Br,1))/4,Br(:,x),'.','color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)],'HandleVisibility','off');%'linewidth',2,'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)],'HandleVisibility','off');
    hold on
    plot((1:size(Bl,1))/4,Bl(:,x),'.','color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)], 'DisplayName', ['x' num2str(n_rectx)]);%'linewidth',2,'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)], 'DisplayName', ['x' num2str(n_rectx)]);
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
    title(['x = ' ' ' num2str(x-(n_rectx+1)/2)])
    legend('-DynamicLegend', 'Location', 'EastOutside');   
    
    ii=ii+1;
    end
    
%end %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< comment

%% Make Trajectories

        
        x_num = (n_rectx-1)/2;
        x = nan*ones(size(runs,1),size(runs,2));
        y = nan*ones(size(runs,1),size(runs,2));
        x(:,1)=(n_rectx+1)/2;
        
        ts=[];
        xs=[];
        bls=[];
        brs=[];
        bws=[];
        counter=1;
        for i = 1:size(Bl,2)
            bls = [bls; Bl(:,i)];
            brs = [brs; Br(:,i)];
            bws = [bws; Blr(:,i)];
            
            ts = [ts 1:size(Bl,1)];
            xs = [xs counter*ones(1,size(Bl,1))];
            counter = counter+1;
        end
        
        clear last
        for traj = 1:size(stimuli,1)
            
            for in = 1:(length(x(traj,:))) %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< timestep ?
                if x(traj,in)==1 || x(traj,in)==n_rectx ||  isnan(x(traj,in))
                    
                else
                    if stimuli(traj,in)>Br(in,x(traj,in)) || (x(traj,in)-(n_rectx+1)/2>0&& isnan(Br(in,x(traj,in)))) 
                        x(traj,in+1) = x(traj,in)+1;
                    elseif stimuli(traj,in)<Bl(in,x(traj,in)) || (x(traj,in)-(n_rectx+1)/2<0 && isnan(Bl(in,x(traj,in)))) 
                        x(traj,in+1) = x(traj,in)-1;
                    elseif isnan(Br(in,x(traj,in))) && x(traj,in)==(n_rectx+1)/2 
                        x(traj,in+1)=x(traj,in)+(2*(rand(1,1)>.5)-1);
                    elseif (n_rectx-1)/2 >= t_num-1
                         x(traj,in+1)=x(traj,in)+(2*(rand(1,1)>.5)-1);
                                          
                    else 
                        x(traj,in+1)=x(traj,in);

                    end
                    
                end
            end
            if any(find(isnan(x(traj,:)),1,'first'))
            last(traj) = find(isnan(x(traj,:)),1,'first')-1;
            else 
                last(traj) = length(x(traj,:));
            end
            choice(traj)=x(traj,last(traj));
 
        end
 %% Calculate Performance
 
 %stim2 is the experienced stimulus
      stim2=nan(size(stim));
      for ind = 1:size(stim2,1)
          stim2(ind,1:last(ind))=stim(ind,1:last(ind));
      end
      %stim2(stim2==0)=-1;
  
        
  means = mean(stim,2);
  for i=1:length(last)
      x(i,last(i))
      if means(i)<0 && x(i,last(i))==1
          responses(i)=1;
      elseif means(i)>0 && x(i,last(i))==1
          responses(i)=2;
      elseif means(i)<0 && x(i,last(i))==n_rectx
          responses(i)=2;
      elseif means(i)>0 && x(i,last(i))==n_rectx
          responses(i)=1;
      end
  end
          
  pc(tout)=length(find(responses==1))/length(responses);  
  
%   rts = last; %unique(last);
%   hold on
%   pcFixed=[];
%   pcFixed2=[];
%   for i=1:length(rts)
%    inds=find((last==rts(i))==1);
% 
%   
%    traces = x(inds,:);
%    lasttrace = traces(:,rts(i));
%    numerator=length(find((lasttrace==1 & nanmean(stim2(inds,:),2)<0) | ((lasttrace==n_rectx & mean(stim2(inds,:),2)>0)) ));
%              
%    pcFixed(i)=numerator;
%    pcFixed2(i)=length(inds);
%   
%   end
  

  means_e = nanmean(stim2,2);
  for i=1:length(last)
      if means_e(i)<0 && x(i,last(i))==1
          responses_e(i)=1;
      elseif means_e(i)>0 && x(i,last(i))==1
          responses_e(i)=2;
      elseif means_e(i)<0 && x(i,last(i))==n_rectx
          responses_e(i)=2;
      elseif means_e(i)>0 && x(i,last(i))==n_rectx
          responses_e(i)=1;
      end
  end
          
  pc_e(tout)=length(find(responses_e==1))/length(responses_e);
  tout=tout+1;

%% Make Plots

figure(5)
subplot(1,2,1)
hold on
plot(n_rectx,pc(tout-1),'.','markersize',40,'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)], 'DisplayName', ['x' num2str(n_rectx)]);
title('%Correct (p)')
xlabel('Arena Width')
ylabel('% Correct')
ylim([0 1])
set(gca,'linewidth',2)
set(gca, 'fontsize',18)
box off
legend('-DynamicLegend', 'Location', 'EastOutside');   

all_pcs1(dim1,dim2)=pc(tout-1);

subplot(1,2,2)
hold on
plot(n_rectx,pc_e(tout-1),'.','markersize',40,'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)], 'DisplayName', ['x' num2str(n_rectx)]);
title('%Correct (Experienced)')
xlabel('Arena Width')

all_pcs2(dim1,dim2)=pc_e(tout-1);

ylim([0 1])
set(gca,'linewidth',2)
set(gca, 'fontsize',18)
box off
legend('-DynamicLegend', 'Location', 'best');   %'EastOutside');   

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

    plot(.25*[1:last(i)],1.02*(stim2(i,1:last(i))+max(stim2(:)))/2, 'k.', 'markersize', 20)

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
whichstim=idx(end-(37))';
for i = whichstim
    hold off
    yyaxis right
    plot(.25*[1:last(i)],1.02*(stim2(i,1:last(i))+max(stim2(:)))/2, 'k.', 'markersize', 20)
    hold on
    plot(.25*[1:last(i)],stimuli(i,1:last(i)), 'r-', 'markersize', 20)
    ylabel('Belief')
    ylim([-0.05 1.05])
    yyaxis left
    hold off
if choice(i)==1
plot(.25*[1:last(i)],x(i,1:last(i))-x_num-1,'r','linewidth',2)
elseif choice(i)==n_rectx
plot(.25*[1:last(i)],x(i,1:last(i))-x_num-1,'b','linewidth',2)
end

xlim([0 n_recty])
ylim([-x_num x_num])

xlabel('Time (s)')
ylabel('x')
title('Hard')
set(gca,'fontsize',18)
set(gca,'linewidth',2)
box off
legend('Trajectory', 'Stimulus','Belief', 'Location', 'southeast')
end 


% fig=figure;
% set(fig,'color','white')      
start = [];
npresses=[];
for i = 1:size(x,1)
start(i)=find(abs(diff(x(i,:)))>0,1,'first');
npresses(i)=nansum(abs(diff(x(i,:)-x(i,1))));
end
% plot(difficulty, npresses/((n_rectx-1)/2),'.','markersize',40)
% hold on
% st=regstats(npresses/((n_rectx-1)/2),difficulty);
% if st.tstat.pval(2)<.05
% h=    plot(difficulty,st.beta(1)+st.beta(2)*difficulty,'r','linewidth',3);
% 
% elseif st.tstat.pval(2)>.05
% h=         plot(difficulty,st.beta(1)+st.beta(2)*difficulty,'b','linewidth',3);
% end
% %  legend(h,['\beta = ' ' '  num2str(round(st.tstat.beta(2)*100)/100),'p=' num2str(round(st.tstat.pval(2)*10000)/10000)],'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)]);
% 
% xlabel('|(NR-NL)/(NR+NL)|')
% ylabel('Number of Presses/Required')
% title(['N Presses = ' ' ' num2str((n_rectx-1)/2)])
% set(gca,'linewidth',2)
% set(gca, 'fontsize',18)
% box off




figure(1)
hold all
plot(n_rectx, mean(start/4),'.','markersize',40,'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)], 'DisplayName', ['x' num2str(n_rectx)]);
errorbar(n_rectx, mean(start/4),std(start/4,[],2),'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)],'linewidth',2,'HandleVisibility','off');

xlabel('Arena Width')
ylabel('Start Time (s)')
title('Start Moving Time vs Width')
set(gca,'linewidth',2)
set(gca, 'fontsize',18)
box off
legend('-DynamicLegend', 'Location', 'EastOutside');   

all_starts(dim1,dim2) = mean(start/4);

figure(2)
hold all
plot(n_rectx, mean(npresses/((n_rectx-1)/2)),'.','markersize',40,'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)], 'DisplayName', ['x' num2str(n_rectx)]);
errorbar(n_rectx,mean(npresses/((n_rectx-1)/2)),std(npresses/((n_rectx-1)/2)),'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)],'HandleVisibility','off');

all_npresses(dim1,dim2) = mean(npresses/((n_rectx-1)/2));

xlabel('Arena Width')
ylabel('Number of Presses/Required')
title('Hesitation vs Arena Width')
set(gca,'linewidth',2)
set(gca, 'fontsize',18)
box off
legend('-DynamicLegend', 'Location', 'EastOutside');   


figure(3)
hold all

st=regstats(npresses./((n_rectx-1)/2),difficulty);
plot(n_rectx, round(st.tstat.beta(2)*100)/100,'.','markersize',40,'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)], 'DisplayName', ['x' num2str(n_rectx)]);

xlabel('Arena Width')
ylabel('\beta')
title('\beta vs Width')
set(gca,'linewidth',2)
set(gca, 'fontsize',18)
box off

all_betas(dim1,dim2) = round(st.tstat.beta(2)*100)/100;
legend('-DynamicLegend', 'Location', 'EastOutside');   

figure(6)
hold all      
      
x2=x;
for i = 1:size(x2,1)
x2(i,1:start(i))=nan;
end
x2=[nan(size(x2,1),1) x2];

x2=diff(x2,1,2);
Zero_idxs=[];
for i = 1:size(x2,1)
    idx = (x2(i,1:end-1) == 0) & (x2(i,2:end) ~= 0 );
    N = sum(idx);
    if (x2(end-1) ~= 0) & (x2(end) == 0)
        N = N + 1;
    end
    Zero_idxs(i)=N;
end
Zero_idxs(Zero_idxs==0)=nan;
  plot(n_rectx, nanmean(Zero_idxs),'.','markersize',40,'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)], 'DisplayName', ['x' num2str(n_rectx)]);
    errorbar(n_rectx, nanmean(Zero_idxs),nanstd(Zero_idxs,[],2),'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)],'linewidth',2,'HandleVisibility','off');

    xlabel('Arena Width')
    ylabel('Average # Post-Start Waits')
    title('Average # of Post-Start Waits given a Wait')
      set(gca,'linewidth',2)
  set(gca, 'fontsize',18)
  box off
      legend('-DynamicLegend', 'Location', 'EastOutside');   

  figure(7)
  hold all
  Zero_idxs(isnan(Zero_idxs))=0;
  plot(n_rectx, mean(Zero_idxs),'.','markersize',40,'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)], 'DisplayName', ['x' num2str(n_rectx)]);
    errorbar(n_rectx, mean(Zero_idxs),std(Zero_idxs,[],2),'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)],'linewidth',2,'HandleVisibility','off');

    xlabel('Arena Width')
    ylabel('Average # of Post-Start Waits')
    title('Post-Start Waiting vs Arena Width')
      set(gca,'linewidth',2)
  set(gca, 'fontsize',18)
  box off
    legend('-DynamicLegend', 'Location', 'EastOutside');   
  
  all_poststarts(dim1,dim2) = mean(Zero_idxs);
  
  
mainfig=figure;
%set(mainfig,'position',[0 100 1200 500]);
set(mainfig,'color','white') 
plot(abs(difficulty), start/4,'.','markersize',40)
hold on
st=regstats(start/4,abs(difficulty));
if st.tstat.pval(2)<.05
h=    plot(abs(difficulty),st.beta(1)+st.beta(2)*abs(difficulty),'r','linewidth',3);

elseif st.tstat.pval(2)>.05
h=         plot(abs(difficulty),st.beta(1)+st.beta(2)*abs(difficulty),'b','linewidth',3);
end
legend(h,['\beta = ' ' '  num2str(round(st.tstat.beta(2)*100)/100),'p=' num2str(round(st.tstat.pval(2)*10000)/10000)]);

xlabel('|(NR-NL)/(NR+NL)|')
ylabel('Start Time (s)')
title(['Arena Width = ' ' ' num2str(n_rectx)])
set(gca,'linewidth',2)
set(gca, 'fontsize',18)
box off    

  
mainfig=figure;
%set(mainfig,'position',[0 100 1200 500]);
set(mainfig,'color','white') 
hold off
plot(abs(difficulty), last/4,'.','markersize',40)
hold on
st=regstats(last/4,abs(difficulty));
if st.tstat.pval(2)<.05
h=    plot(abs(difficulty),st.beta(1)+st.beta(2)*abs(difficulty),'r','linewidth',3);

elseif st.tstat.pval(2)>.05
h=         plot(abs(difficulty),st.beta(1)+st.beta(2)*abs(difficulty),'b','linewidth',3);
end
  legend(h,['\beta = ' ' '  num2str(round(st.tstat.beta(2)*100)/100),'p=' num2str(round(st.tstat.pval(2)*10000)/10000)]);

xlabel('|(NR-NL)/(NR+NL)|')
ylabel('Reaction Time (s)')
title(['Arena Width = ' ' ' num2str(n_rectx)])
set(gca,'linewidth',2)
set(gca, 'fontsize',18)
box off    

mainfig=figure;
%set(mainfig,'position',[0 100 1200 500]);
set(mainfig,'color','white') 
hold off
plot(difficulty, last/4,'.','markersize',40)
hold on
st=regstats(last/4,difficulty);
% if st.tstat.pval(2)<.05
% h=    plot(difficulty,st.beta(1)+st.beta(2)*difficulty,'r','linewidth',3);
% 
% elseif st.tstat.pval(2)>.05
% h=         plot(difficulty,st.beta(1)+st.beta(2)*difficulty,'b','linewidth',3);
% end
%  legend(h,['\beta = ' ' '  num2str(round(st.tstat.beta(2)*100)/100),'p=' num2str(round(st.tstat.pval(2)*10000)/10000)],'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)]);
xlim([-1 1])
xlabel('(NR-NL)/(NR+NL)')
ylabel('Reaction Time (s)')
title(['Arena Width = ' ' ' num2str(n_rectx)])
set(gca,'linewidth',2)
set(gca, 'fontsize',18)
box off    

  
mainfig=figure;
%set(mainfig,'position',[0 100 1200 500]);
set(mainfig,'color','white') 
hold off
plot(abs(difficulty), npresses,'.','markersize',40)
hold on
st=regstats(npresses,abs(difficulty));
if st.tstat.pval(2)<.05
h=    plot(abs(difficulty),st.beta(1)+st.beta(2)*abs(difficulty),'r','linewidth',3);

elseif st.tstat.pval(2)>.05
h=         plot(abs(difficulty),st.beta(1)+st.beta(2)*abs(difficulty),'b','linewidth',3);
end
legend(h,['\beta = ' ' '  num2str(round(st.tstat.beta(2)*100)/100),'p=' num2str(round(st.tstat.pval(2)*10000)/10000)]);

xlabel('|(NR-NL)/(NR+NL)|')
ylabel('Number of Presses')
title(['Arena Width = ' ' ' num2str(n_rectx)])
set(gca,'linewidth',2)
set(gca, 'fontsize',18)
box off    


choices = choice;
choices(choices == 1) = 2;

choices(choices == n_rectx) = 1;
choices = choices-1;

figure(8)
%plot((difficulty), choices,'.','markersize',40)
hold all
st = glmfit((difficulty),choices','binomial','link','logit');
z =  st(1) + ((difficulty) * st(2));
z = 1 ./ (1 + exp(-z));

h =    plot((difficulty),z,'.','markersize',15, 'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)],'DisplayName',[ ['x' num2str(n_rectx)] ' \beta = ' ' '  num2str(round(st(2)*100)/100)]);

legend('-DynamicLegend', 'Location', 'EastOutside');

%legend(h,['\beta = ' ' '  num2str(round(st(2)*100)/100)]);
xlim([-1 1])

xlabel('(NR-NL)/(NR+NL)')
ylabel('Proportion Left')
title(['Arena Width = ' ' ' num2str(n_rectx)])
set(gca,'linewidth',2)
set(gca, 'fontsize',18)
box off    


figure(9)
hold all
plot(n_rectx, mean(last/4),'.','markersize',40,'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)],'DisplayName',[ ['x' num2str(n_rectx)]]);
errorbar(n_rectx, mean(last/4),std(last/4,[],2),'color', [n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1) n_rectx/(max(n_rectxs)+1)],'linewidth',2,'HandleVisibility','off');

xlabel('Arena Width')
ylabel('Reaction Time (s)')
title('Reaction Time vs Arena Width')
set(gca,'linewidth',2)
set(gca, 'fontsize',18)
box off
legend('-DynamicLegend', 'Location', 'EastOutside');   

end
%%
h = get(0,'children');
%h = sort(h);
for j=1:length(h)
  saveas(h(j), ['AllAgentx' num2str(n_recty) '_tc_'   num2str(t_cost) '_' 'xc_' num2str(x_cost) '_' num2str(j) 'alpha_' num2str(alpha) '.fig']);
  saveas(h(j), ['AllAgentnx' num2str(n_recty) '_tc_'   num2str(t_cost) '_' 'xc_' num2str(x_cost) '_' num2str(j) 'alpha_' num2str(alpha) '.pdf']);
  saveas(h(j), ['AllAgentnx' num2str(n_recty) '_tc_'  num2str(t_cost) '_' 'xc_' num2str(x_cost) '_' num2str(j) 'alpha_' num2str(alpha) '.png']);
end

close all
end %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<uncomment

computation_time = toc(maintic);
toc(maintic)

%%
[X,Y]=meshgrid(n_rectxs,([ 16*14 8*14 4*14 2*14 14  ]));
fig=figure;
set(fig,'color','white')
surf(X,Y,all_npresses(1:5,:),'FaceAlpha',0.5)
zlim([1 2])
xlabel('x')
ylabel('y')
zlabel('N presses/Required')

fig=figure;
set(fig,'color','white')
surf(X,Y,all_starts(1:5,:),'FaceAlpha',0.5)
zlim([0 2.5])

xlabel('x')
ylabel('y')
zlabel('Start')


fig=figure;
set(fig,'color','white')
surf(X,Y,all_betas(1:5,:),'FaceAlpha',0.5)
zlim([-.2 0])

xlabel('x')
ylabel('y')
zlabel('\beta')


fig=figure;
set(fig,'color','white')
surf(X,Y,all_poststarts(1:5,:),'FaceAlpha',0.5)
zlim([0 2])

xlabel('x')
ylabel('y')
zlabel('Post Start Wait')

fig=figure;
set(fig,'color','white')
surf(X,Y,all_pcs1(1:5,:),'FaceAlpha',0.5)

xlabel('x')
ylabel('y')
zlabel('%Correct')
zlim([.6 .9])


[X,Y]=meshgrid(n_rectxs,[13:-1:1]);% 14 2*14 4*14 8*14 16*14]);
fig=figure;
set(fig,'color','white')
surf(X,Y,all_npresses(6:end,:),'FaceAlpha',0.5)
zlim([1 2])

xlabel('x')
ylabel('y')
zlabel('N presses/Required')

fig=figure;
set(fig,'color','white')
surf(X,Y,all_starts(6:end,:),'FaceAlpha',0.5)
zlim([0 2.5])

xlabel('x')
ylabel('y')
zlabel('Start')

fig=figure;
set(fig,'color','white')
surf(X,Y,all_betas(6:end,:),'FaceAlpha',0.5)
zlim([-.2 0])

xlabel('x')
ylabel('y')
zlabel('\beta')

fig=figure;
set(fig,'color','white')
surf(X,Y,all_poststarts(6:end,:),'FaceAlpha',0.5)
zlim([0 2])

xlabel('x')
ylabel('y')
zlabel('Post Start Wait')

fig=figure;
set(fig,'color','white')
surf(X,Y,all_pcs1(6:end,:),'FaceAlpha',0.5)

xlabel('x')
ylabel('y')
zlabel('%Correct')
zlim([.6 .9])

h = get(0,'children');
%h = sort(h);
for j=1:length(h)
  saveas(h(j), ['AllAgentx_arena' num2str(n_recty) '_tc_'   num2str(t_cost) '_' 'xc_' num2str(x_cost) '_' num2str(j) 'alpha_' num2str(alpha) '.fig']);
  print(h(j), ['AllAgentnx_arena' num2str(n_recty) '_tc_'   num2str(t_cost) '_' 'xc_' num2str(x_cost) '_' num2str(j) 'alpha_' num2str(alpha) '.pdf'], '-fillpage', '-dpdf');
  saveas(h(j), ['AllAgentnx_arena' num2str(n_recty) '_tc_'  num2str(t_cost) '_' 'xc_' num2str(x_cost) '_' num2str(j) 'alpha_' num2str(alpha) '.png']);
end
  save(['AllAgentnx_arenaTrajs' num2str(n_recty) '_tc_'  num2str(t_cost) '_' 'xc_' num2str(x_cost) '_' num2str(j) 'alpha_' num2str(alpha) '.mat'],'x','stim','stim2','stimuli','p','Bl','Br');

  
  close all

[X,Y]=meshgrid(n_rectxs,[16*14 8*14 4*14 2*14 14  13:-1:1]);% 14 2*14 4*14 8*14 16*14]);
fig=figure;
set(fig,'color','white')
surf(X,Y,all_npresses(1:end,:),'FaceAlpha',0.5)
zlim([1 2])

xlabel('x')
ylabel('y')
zlabel('N presses/Required')

fig=figure;
set(fig,'color','white')
surf(X,Y,all_starts(1:end,:),'FaceAlpha',0.5)
zlim([0 2.5])

xlabel('x')
ylabel('y')
zlabel('Start')

fig=figure;
set(fig,'color','white')
surf(X,Y,all_betas(1:end,:),'FaceAlpha',0.5)
zlim([-.2 0])

xlabel('x')
ylabel('y')
zlabel('\beta')

fig=figure;
set(fig,'color','white')
surf(X,Y,all_poststarts(1:end,:),'FaceAlpha',0.5)
zlim([0 2])

xlabel('x')
ylabel('y')
zlabel('Post Start Wait')

fig=figure;
set(fig,'color','white')
surf(X,Y,all_pcs1(1:end,:),'FaceAlpha',0.5)

xlabel('x')
ylabel('y')
zlabel('%Correct')
zlim([.6 .9])

h = get(0,'children');
%h = sort(h);cl
for j=1:length(h)
  saveas(h(j), ['arena' num2str(n_recty) '_tc_'   num2str(t_cost) '_' 'xc_' num2str(x_cost) '_' num2str(j) 'alpha_' num2str(alpha) '.fig']);
  saveas(h(j), ['arena' num2str(n_recty) '_tc_'   num2str(t_cost) '_' 'xc_' num2str(x_cost) '_' num2str(j) 'alpha_' num2str(alpha) '.pdf']);
  saveas(h(j), ['arena' num2str(n_recty) '_tc_'  num2str(t_cost) '_' 'xc_' num2str(x_cost) '_' num2str(j) 'alpha_' num2str(alpha) '.png']);
end


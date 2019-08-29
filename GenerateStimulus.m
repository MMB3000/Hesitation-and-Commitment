function [ stim ] = GenerateStim( nruns, n_recty )
%sample stimuli from our task

                   %Penalty for incorrect responses

%% SET STIMULUS PARAMETERS
ILDSo=[-6,6];       %in dB
%prob_values = [(22/28-8/28) (18/28-12/28) (16/28-14/28)];
prob_values = [.2  .2 .2];
mean_shifts = [ILDSo(1)*prob_values(1), ILDSo(1)*prob_values(2), ILDSo(1)*prob_values(3),...
               ILDSo(2)*prob_values(3),ILDSo(2)*prob_values(2),ILDSo(2)*prob_values(1)];
AverageBinaural_Level=60;     %in dB4+4+4++
neighbor = 0;

stim_duration_sec=n_recty;    %in seconds
ramp_stim_dur_msec=0;         %in miliseconds

chunk_duration_msec=250;      %in miliseconds (includes rise and fall time)
filter_win_msec=30;           %Cosine filtering window in miliseconds
inter_burst_time = 150;       %Stimulus processing time
refractoryT = 0;              %Refractory period on button presses

%Fs=44100;                       %Sample Rate in Hz

% %% SET EXPERIMENT CONDITIONS
% cts = [0];
% n_rectys = [n_recty n_recty n_recty];
% 
% trials_per_condition=12;           %must be a multiple of 3
%  
% %Conditions are predetermined and shuffled
% blocknumber='1';
% seed_is=str2num(blocknumber);
% if seed_is>5  &&  seed_is<11
% 
%     seed_is = seed_is-5;
% elseif seed_is>10 && seed_is<16
%     
%    seed_is = seed_is-10;
%    
% elseif seed_is>15
%    seed_is=seed_is-15;
%    
% end
% rng(seed_is);    
% cond_1=ones(1,trials_per_condition);
% cond_2=2*ones(1,trials_per_condition);
% cond_3=3*ones(1,trials_per_condition);
% cond_4=4*ones(1,trials_per_condition);
% cond_5=5*ones(1,trials_per_condition);
% cond_6=6*ones(1,trials_per_condition);
% % cond_7=7*ones(1,trials_per_condition);
% % cond_8=8*ones(1,trials_per_condition);
% 
% conditions_are=horzcat(cond_1,cond_2,cond_3,cond_4,cond_5,cond_6);%,cond_7,cond_8);
% i_rand=randperm(length(conditions_are));
% conditions_shuffled=conditions_are(i_rand);
%     
% repeat_3versions=length(conditions_are)/3;
% unshuffled_versions=[];
% for num=1:repeat_3versions
%     uno_23=[1 2 3];
%     unshuffled_versions=[unshuffled_versions uno_23];
% end
% shuffled_3=unshuffled_versions(i_rand); %specifies 1 of the 4 versions
% presets_order=shuffled_3;  %saving name for shuffled_4


%% 4 DIFFERENT VERSIONS OF ILD SEQUENCE
% %ILD sequence condition 1
ILDS=ILDSo;%+mean_shifts(1);
% good=0;
% while good==0
%     ILD_preset_1 = any_chunks_gen_uniform(ILDS,mean_shifts(1), AverageBinaural_Level, chunk_duration_msec, stim_duration_sec, 1, filter_win_msec, neighbor, 1);
%     if abs(mean_shifts(1)-mean(ILD_preset_1))<eps
%         good=1;
%     end
% end
% 
% 
% %ILD sequence condition 2
% ILDS=ILDSo;%+mean_shifts(2);
% good=0;
% while good==0 
%     ILD_preset_2 = any_chunks_gen_uniform(ILDS,mean_shifts(2), AverageBinaural_Level, chunk_duration_msec, stim_duration_sec, 1, filter_win_msec, neighbor, 1);
%     if  abs(mean_shifts(2)-mean(ILD_preset_2)) < eps
%         good=1;
%     end
% end 
% 
% 
% %ILD sequence condition 3
% ILDS=ILDSo;%+mean_shifts(3);
% good=0;
% while good==0  
%     ILD_preset_3 = any_chunks_gen_uniform(ILDS, mean_shifts(3), AverageBinaural_Level, chunk_duration_msec, stim_duration_sec, 1, filter_win_msec, neighbor,1);
%     if  abs(mean_shifts(3)-mean(ILD_preset_3)) < eps
%         good=1;
%     end
% end 


% %ILD sequence condition 4
% ILD_preset_4=-ILD_preset_3;
% %ILD sequence condition 5
% ILD_preset_5=-ILD_preset_2;
% %ILD sequence condition 6
% ILD_preset_6=-ILD_preset_1;
% 
% number_of_chunks=length(ILD_preset_6);
% 
% ILD_presets_1=vertcat(ILD_preset_1(1,:),ILD_preset_2(1,:), ILD_preset_3(1,:),ILD_preset_4(1,:), ILD_preset_5(1,:), ILD_preset_6(1,:));%, ILD_preset_7(1,:), ILD_preset_8(1,:));
% ILD_presets_1 = (ILD_presets_1/6+1)/2;
% 
% 
% 
% %sigma =.01;
% tnum=length(ILD_presets_1);

%%
tnum=n_recty*1/(chunk_duration_msec/1000);%any_chunks_gen_uniform(ILDS, mean_shifts(3), AverageBinaural_Level, chunk_duration_msec, stim_duration_sec, 1, filter_win_msec, neighbor,1);

stim = nan(nruns,tnum);


for i = 1:nruns
    
    imod = mod(i,6);
    if imod ==0
        imod = 6;  
    end
    
    if imod<4
        stim(i,:)=ILDS(1)*(2*discreternd(tnum,[.54 .46])-3);
    else
        stim(i,:)=ILDS(1)*(2*discreternd(tnum,[.46 .54])-3);
    end
%     stim(i,:)=gcurr(randperm(lg));
%     stim(i,:)=any_chunks_gen_uniform(ILDS, mean_shifts(3), AverageBinaural_Level, chunk_duration_msec, stim_duration_sec, 1, filter_win_msec, neighbor,1);

end

end


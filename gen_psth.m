function [sdf,sdf_final]= gen_psth(post_process, working_dir)
%%
close all
for i=1:length(post_process.rast.time_stamps)
    combo(i).time_stamps=str2double(post_process.rast.time_stamps(i,:));
    combo(i).grasp=char(post_process.rast.grasp(i,:));
    combo(i).movetype=char(post_process.rast.type(i,:));
    combo(i).hit=char(post_process.rast.hit(i,:));
    combo(i).block=str2double(post_process.rast.block(i,:));
    combo(i).slc_times=str2double(post_process.rast.slc_times(i,:));
end
%% SELECT RELEVANT TRIALS
clear i
j=1;
for i = 1:length(combo)
    combo(i).adj_time=combo(i).time_stamps-combo(i).time_stamps(3);
    if combo(i).adj_time(end)>1.5 && strcmp(combo(i).movetype,'close') && ...
            combo(i).adj_time(1)<-3
        psth(j)=combo(i);
        j=j+1;
    else
        %disp(['skip trial ' num2str(i)])
    end
end
clear i j

%% USE ONLY THE RELVEANT TRIALS
clear i a block_num gauss_smooth j nsp_data nsp_st p1 p2 sigma tstep
sdf_raw={};
sdf={};
grasp_time=(-1.999:.001:1);

for i =1:length(psth) %i is the trial number
    %only use channels with enough info on hand closing
    block_num=psth(i).block;
    nsp_data=post_process.nsp.dataset(block_num).sct;

    %sample for first channel
    for j=1:128 %j is the channel number
        nsp_spikes_all=nsp_data{j,1};
        beg=round(psth(i).time_stamps(1),3);
        cue=round(psth(i).time_stamps(2),3);
        go=round(psth(i).time_stamps(3),3);
        fin=round(psth(i).time_stamps(end),3);
         
        idx=find(nsp_spikes_all>=beg & nsp_spikes_all<=fin);
        nsp_st=round(nsp_spikes_all(idx),3);
        binary_spikes=[];
        
        if isempty(nsp_st)
            sdf_raw{j,i}=zeros([1,3000]);
        else
            %Gaussian Smoothing within the trial
            binary_spikes=[];
            tstep=.001; % Resolution for SDF [s]
            sigma=.05; % Width of gaussian/window [s] (50ms)
            t_gauss=(-.5:tstep:.5);
            mu=0;
            p1=(-0.5)*((t_gauss-mu)/sigma).^2;
            p2=((sigma^2)*2*pi);
            gauss_func=exp(p1)./p2;
            time=beg:tstep:fin; 
            time=round(time,3);

            for k= 1:length(time)
                if find(time(k)==nsp_st)
                    binary_spikes(k)=1;
                else
                    binary_spikes(k)=0;
                end
            end
            a=conv(binary_spikes,gauss_func,'same');
            hold_time=find(time==cue);
            pre_time=find(time==cue);
            move_time=find(time==go);
            
            %select out three second intervals
            hold_epoch=a((hold_time-1000):(hold_time-1));
            pre_epoch=a(pre_time:(pre_time+999));
            move_epoch=a(move_time:(move_time+999));
            total=[hold_epoch, pre_epoch, move_epoch];
            
            raster{j,i}=binary_spikes;
            raster{129,i}=psth(i).grasp;
            sdf_raw{j,i}=total;
            sdf_raw{129,i}=psth(i).grasp;

%             if j ==83
%             subplot(2,2,1)
%             hold on
%             xline(beg,'Linewidth',2,'Color','r')
%             xline(cue,'Linewidth',2,'Color','r')
%             xline(go,'Linewidth',2,'Color','r')
%             xline(fin,'Linewidth',2,'Color','r')
%             area([cue-1, cue+1],[1.05 1.05],'FaceColor','#D4E6F1')
%             area([go, go+1],[1.05 1.05],'FaceColor','#D4E6F1')
%             plot(time,binary_spikes,'LineWidth',1,'Color','k')
%             title('Sample Trial Raster Plot')
%             ylim([.99 1])
% 
%             subplot(2,2,2)
%             title('Sample Trial Smoothed Data')
%             hold on          
%             area([cue-1, cue+1],[350, 350],'FaceColor','#D4E6F1')
%             area([go, go+1],[350, 350],'FaceColor','#D4E6F1')
%             xline(beg,'Linewidth',1,'Color','r')
%             xline(cue,'Linewidth',1,'Color','r')
%             xline(go,'Linewidth',1,'Color','r')
%             xline(fin,'Linewidth',1,'Color','r')
%             plot(time, a, 'Color','k','Linewidth',1.5)  
% 
% 
%             subplot(2,2,3)
%             title('Sample Trial 3 Second Selection')
%             hold on
%             area([-2, -1],[350 350],'FaceColor','#E8DAEF')
%             xline(-1,'Linewidth',1,'Color','r')
%             xline(0,'Linewidth',1,'Color','r')
%             xline(-2,'Linewidth',1,'Color','r')
%             xline(1,'Linewidth',1,'Color','r')
%             plot(grasp_time,total,'Color','k','Linewidth',1.5)
% 
%             subplot(2,2,4)
%             plot(grasp_time,z_trial,'Color','k','Linewidth',1.5)
%             title('Sample Trial Z-Scored')
%             hold on
%             xline(-1,'Linewidth',1,'Color','r')
%             xline(0,'Linewidth',1,'Color','r')
%             xline(-2,'Linewidth',1,'Color','r')
%             xline(1,'Linewidth',1,'Color','r')
%             %spikes{j,i}=binary_spikes;
%             %full_trial{j,i}=a;
%             disp('saving figure')
%             saveas(gcf,[working_dir filesep 'TrialDemo.jpg'])
%             end
% 
%             figure(3)
%             hold on
%             xline(beg,'Linewidth',2,'Color','r')
%             xline(cue,'Linewidth',2,'Color','r')
%             xline(go,'Linewidth',2,'Color','r')
%             xline(fin,'Linewidth',2,'Color','r')
%             area([cue-1, cue+1],[1.05 1.05],'FaceColor','#D4E6F1')
%             area([go, go+1],[1.05 1.05],'FaceColor','#D4E6F1')
%             plot(time,binary_spikes*j,'|','Color','k')
%             title('Sample Trial Raster Plot')


        end
    end
end
close all
clear i a block_num gauss_smooth j nsp_data nsp_st p1 p2 sigm...
    ans beg binary_spikes cue fin gauss_func go hold_epoch...
    hold_time idx k move_epoch move_time mu nsp_spikes_all...
    t_gauss time total tstep pre_time pre_epoch

%% Z SCORE BY SURROUNDING FOUR HOLD EPOCHS
fail_chan=[];
fail_trial=[];
fail_data=[];
num_trials=size(sdf_raw,2);
for i =1:num_trials %i is the trial number
    for j=1:128
        if i<=2 %if you're at the beginning use just the next one
            trial_select=cell2mat(sdf_raw(j,1:5)');
        elseif i>=(num_trials-1) %if you're at the end use just the previous one
            trial_select=cell2mat(sdf_raw(j,(num_trials-4):num_trials)');
        else
            trial_select=cell2mat(sdf_raw(j,(i-2):(i+2))'); 
            %select one trial on either side
        end
        curr_trial=sdf_raw{j,i};
        avg_by=trial_select(:,1:1000); %select hold epoch
        if find(avg_by) %if the selected data has any spikes surrounding
            z_mean=mean(mean(avg_by));
            z_std=std(mean(avg_by,1));
            z_trial=(curr_trial-z_mean)/z_std;
        else %if not
            z_trial=curr_trial;
        end

%         if isnan(z_trial)
%             fail_chan=[fail_chan, j];
%             fail_trial=[fail_trial, i];
%             fail_data=[fail_data; avg_by];
%         end

        sdf{j,i}=z_trial;

        sdf{129,i}=sdf_raw{129,i};
    end
end

%% GROUP BY GRASP
close all
clear avg_grasp chan_num colors empt full_trial g gnames grasp grasp_avg...
    grasp_matrix grasp_out grasp_time hold_time hold_epoch i j mems sigma...
    z_grasp_avg z_mean z_std
gnames={'power' 'palmer' 'pinch' 'disc' 'lateral'};
colors=jet(length(gnames))*0.8;
grasp_time=(-1.999:.001:1);

for i = 1:128
    chan_num=i;
    for g= 1:length(gnames)
        grasp_matrix=[];
        grasp=gnames(g);
        empt={};
        for k = 1:length(psth)
            empt{k}=psth(k).grasp;
        end
        mems=find(ismember(empt,grasp));

        grasp_trials=sdf(i,mems);
        grasp_matrix=[];
        for j = 1:length(grasp_trials)
            grasp_matrix(j,:)=grasp_trials{1,j};
        end

%         figure(3)
%         sgtitle('Trials by Grasp')
%         subplot(5,1,g)
%         for j = 1:length(grasp_trials)
%             plot(grasp_time,grasp_matrix(j,:),'Color',colors(g,:))
%             title(grasp)
%             hold on
%         end

        grasp_avg=mean(grasp_matrix,1);

        %matrix needed to run anova
        sdf_final{i,g}=grasp_trials;

        %use this to actually plot the psth
        grasp_out{i,g}=grasp_avg;

% 
%         if i<65
%            figure(4)
%            subplot(8,8,chan_num);                      
%            plot(grasp_time,grasp_avg,'Color',colors(g,:))
%            hold on
%            xline(0,'Linewidth',1,'Color','k')
%            xline(-1,'Linewidth',1,'Color','k')
%            title(num2str(chan_num))
%        elseif i>64
%            figure(5)
%            subplot(8,8,chan_num-64);                       
%            plot(grasp_time,grasp_avg,'Color',colors(g,:))
%            hold on
%            xline(0,'Linewidth',1,'Color','k')
%            xline(-1,'Linewidth',1,'Color','k')
%            title(num2str(chan_num))
%         end
    end
end
% figure(4)
% set(gcf, 'Position', get(0, 'Screensize'));
% figure(5)
% set(gcf, 'Position', get(0, 'Screensize'));
% saveas(figure(4),[working_dir filesep 'PSTH1.jpg'])
% saveas(figure(5),[working_dir filesep 'PSTH2.jpg'])
% disp('saving figures')
% %% INDIVIDUAL PSTH
% chan_num=83;
% grasp_time=(-1.999:.001:1);
% figure(chan_num)
% for i = 1:5
%     plot(grasp_time,grasp_out{chan_num,i},'Color',colors(i,:),'Linewidth',2)
%     hold on
%     title('Sample Channel PSTH')
%     %[~,~,m,~] = normfit(grasp_out{chan_num,i});
%     %fHandle=errorPatch(grasp_time', m', colors(i,:), 0.2 );
% end
% xline(-1,'Linewidth',1,'Color','r')
% xline(0,'Linewidth',1,'Color','r')
% legend([gnames 'cue' 'movement'],'Location','eastoutside')
% disp('saving figure')
% saveas(figure(chan_num),[working_dir filesep 'ExamplePSTH.jpg'])

end
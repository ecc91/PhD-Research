function [post_process, working_dir]=sort_trials(ol_blocks, cl_blocks,...
    filt_used, box_folder ,nsp_num, loop_id)
% this function is meant to sort trials in each block and output their
% spike times, slc data, nsp data, and other relevant parameters

path=['C:\Users\econo\Box\ReHAB_RP1\Session Data' filesep box_folder];
blocks=[ol_blocks cl_blocks];
nsp_names={'NSP1' 'NSP2' 'NSP3'};
arrayChannelSets = {1:128, 129:192, 193:256, 257:384};
arrayNames = {'Sensory','Parietal','IFG','Motor'};
mkdir([box_folder filesep date filesep loop_id])

working_dir=[box_folder filesep date filesep loop_id];
%% Set the blocks to be analyzed (cl blocks or ol blocks)
if strcmp(loop_id,'open')
    new_blocks=ol_blocks;
elseif strcmp(loop_id,'closed')
    new_blocks=cl_blocks;
end

%% filter information
if strcmp(loop_id,'closed')
    num_filt=max(filt_used);
    rootPath = [path filesep 'Data' filesep 'NCS Data'];
    tmp = dir([rootPath filesep 'Filters' '*(' num2str(num_filt) ')*.mat']);
    fileName = [rootPath filesep tmp.name];
    filt=load(fileName);
    
    figure()
    sgtitle('Features by Filter')
    for i = 1:num_filt
        [~,y,z]=find(filt.sFILTERS(i).sBUILD.filter.K);
        subplot(num_filt,1,i)
        plot(y,z,'.','Markersize',20)
        title(['Filter' num2str(i)])
        xlim([1 768])
        hold on
        xline(1,'-','Sensory - nctx', 'LabelOrientation','horizontal')
        xline(128,'-', 'AIP - nctx','LabelOrientation','horizontal')
        xline(193, '-','IFG - nctx','LabelOrientation','horizontal')
        xline(128*2,'-', 'Motor - nctx','LabelOrientation','horizontal')
        xline(128*3,'-','Sensory - sbp','LabelOrientation','horizontal')
        xline(128*4,'-', 'AIP - sbp','LabelOrientation','horizontal')
        xline(577, '-','IFG - sbp','LabelOrientation','horizontal')
        xline(128*5,'-', 'Motor - sbp','LabelOrientation','horizontal')
        
        Weights=z;
        Channels=y;
        table(Weights,Channels)
        wchan(i).data=[z y];
    end
    disp('saving figure')
    saveas(gcf,[working_dir filesep 'FilterInfo.jpg'])
else
    wchan=[];
end
clear i x y z tmp Weights rootPath num_filt fileName ans nsp2

%% neural data spike threshold crossings
wchan=wchan();
nspName=char(nsp_names(nsp_num));
disp('loading ns5 and slc files')

for i =  1:length(new_blocks)
    % Load .ns5 file and extract spiking times
    blockNum=new_blocks(i)
    fileName = getNS5FileName_RP1(path, nspName, blockNum);
    threshMultiplier = [4.5];

    %car first 64 channels
    carChans = [1:64];
    [spikeTimes1, timeAxis] = getfeatures_RP1( fileName, threshMultiplier, carChans );
    timeWindow = [timeAxis(1) timeAxis(end)];

    %car last 64 channels
    carChans = [64:128];
    [spikeTimes2, timeAxis] = getfeatures_RP1( fileName, threshMultiplier, carChans );

    spikeTimes=[spikeTimes1(1:64) ; spikeTimes2(65:128)];

    nsp.dataset(i).sct=spikeTimes;
    nsp.dataset(i).window=timeWindow;

    clear spikeTimes timeAxis fr fileName threshMultiplier carChans...
        counts timeWindow

    num_filt=filt_used(i);
    rootPath = [path filesep 'Data' filesep 'SLC Data'];
    tmp = dir([rootPath filesep 'SLCdata' '*(' num2str(blockNum) ')*.mat']);
    fileName = [rootPath filesep tmp.name];
    slc(i).data=load(fileName);
    clear x y tmp rawdata

end

disp('finished inputting nsp data')
%% clear unneccesary variables
disp('clearing variables')
clear ap_acutal ap_goal carChans Channels i j k  rootPath ...
    threshMultiplier a act ap b ch_list chg combo dat dh disc ep first...
    g goal grasp l last lateral minl move_start next_ap next_ep palmer ...
    pinch power pre_act pre_goal pre_start rast set_to_zero start ...
    start_idx start_row stop t temp_array trial_end trial_wind tris...
    b blockNum combo dat dh fin first_trial grasp h hits hold_start i j move_start...
    pre_start rng start start_idx temp_array trial_end

%% Divide information by grasp and trial
clear ap_acutal ap_goal carChans Channels i j k  rootPath ...
    threshMultiplier a act ap b ch_list chg combo dat dh disc ep first...
    g goal grasp l last lateral minl move_start next_ap next_ep palmer ...
    pinch power pre_act pre_goal pre_start rast set_to_zero start ...
    start_idx start_row stop t temp_array trial_end trial_wind tris combo


if nsp_num==1
    clock_name='nspClock';
elseif nsp_num==2
    clock_name='nsp2Clock';
elseif nsp_num==3
    clock_name='nsp3Clock';
end

combo=[];

for i = 1:length(new_blocks)
    disp(['block' num2str(new_blocks(i))])

    ap=[];
    ap.actual=slc(i).data.task.auxiliary.values(:,9);
    ap.goal=slc(i).data.task.auxiliary.values(:,43);
    targ=slc(i).data.task.auxiliary.values(:,54);
    hits=slc(i).data.task.auxiliary.values(:,53);
    ap.epoch(1,1:4)="base";

    
    for j = 1:length(ap.actual) %find the epoch diliniations 
        act=ap.actual(j);
        goal=ap.goal(j);
        ap.epoch(j,3)=targ(j); % 0-trial blank, 1-premovement, 2-move
        ap.epoch(j,1)=act;
        ap.epoch(j,2)=goal;
        %find which grasp is being done in each trial
        if goal==0 || goal==100 || goal==52
            ap.epoch(j,4)="power";
        elseif goal==1 || goal==99 || goal==51
            ap.epoch(j,4)="lateral";
        elseif goal==2 || goal==98 || goal==50
            ap.epoch(j,4)="palmer";
        elseif goal==3 || goal==97 || goal==49
            ap.epoch(j,4)="pinch";
        elseif goal==4 || goal==96 || goal==48
            ap.epoch(j,4)="disc";
        end
    end
    clear j  targ act goal


    grasp(i).info=ap.epoch;
    grasp(i).info(1,5)=1;
    grasp(i).info(1,6)=3;
    grasp(i).info(:,8)=hits;

    for j = 1:length(grasp(i).info)-1
        targ=ap.epoch(j,3);
        next_targ=ap.epoch(j+1,3);
        if targ=="2" && next_targ=="0"
             grasp(i).info(j+1,5)=str2double(grasp(i).info(j,5))+1;
        else
            grasp(i).info(j+1,5)=str2double(grasp(i).info(j,5));
        end
        
        if targ=="0" && next_targ=="1" %blank to premove
            grasp(i).info(j+1,6)=1;
        elseif targ=="1" && next_targ=="2" %premove to move
            grasp(i).info(j+1,6)=2;
        elseif targ=="2" && next_targ=="0" %move to blank
            grasp(i).info(j+1,6)=3;
        end

        if str2num(grasp(i).info(j,2))>6
             grasp(i).info(j,7)="close";
        else
            grasp(i).info(j,7)="open";
        end
    end

    clear j

    
%limit to trials that only have hold, premove, move
    first_trial=find(str2double(grasp(i).info(:,3))==0,1);
    
    a=find(grasp(i).info(:,6)=="1");
    b=find(grasp(i).info(:,6)=="2");
    c=find(grasp(i).info(:,6)=="3");
    a=a(a>=first_trial); %only find the ones that occur when the first trial starts
    b=b(b>=first_trial);
    c=c(c>=first_trial);
    if first_trial == c(1)
        c=c;
    else
        c=[first_trial; c];
    end %set the first transition to when the trial starts
    l(1)=length(a);
    l(2)=length(b);
    l(3)=length(c);
    minl=min(l); 

    for j=1:minl
        bounds=[c(j), a(j), b(j)];
        hold_start= c(j);
        pre_start=a(j);
        move_start=b(j);
        grasp(i).time_stamps(j,:)= [hold_start, pre_start, move_start,...
            1, grasp(i).info(pre_start,4),grasp(i).info(pre_start,7),...
            grasp(i).info(pre_start,5)];
    end

    for k =2:minl
        grasp(i).time_stamps(k-1,4)=c(k)-1;
    end
    grasp(i).time_stamps(end,4)=length(grasp(i).info);

    %check to see if there's a hit in that range
    for h=1:size(grasp(i).time_stamps,1)
        fin=str2double(grasp(i).time_stamps(h,4));
        start=fin-24;
        rng=str2double(grasp(i).info(start:fin,8));
        if all(rng==1)
            grasp(i).time_stamps(h,8)="hit";
        else
            grasp(i).time_stamps(h,8)="miss";   
        end
    end
    
    clear j
    clear j k a ans ap ap_actual b bounds c e f hold_start minl...
        move_start next_targ targ rawdata pre_start l num_filt...
        chan_num start_idx
    %find the start idx for each nsp (in this case only nsp2)
    slc_dat=slc(i).data.clocks;
    nm=fieldnames(slc_dat);
    ind=contains(nm, clock_name);
    dat=slc_dat.(nm{ind});
    ns5StartRecordingIdx=find(diff(dat)<0,1);

    if isempty(ns5StartRecordingIdx)
        ns5StartRecordingIdx = dat(1);
    end

    start_idx(i)=double(ceil(ns5StartRecordingIdx*1000));
    for j=1:length(grasp(i).time_stamps)
        grasp(i).time_stamps(j,9)=i;
        %timestamps of 1:4 is where the SLC time stamps are
        grasp(i).time_stamps(j,10)=(str2double(grasp(i).time_stamps(j,1))...
            *20+start_idx(i))/1000;
        grasp(i).time_stamps(j,11)=(str2double(grasp(i).time_stamps(j,2))...
            *20+start_idx(i))/1000;
        grasp(i).time_stamps(j,12)=(str2double(grasp(i).time_stamps(j,3))...
            *20+start_idx(i))/1000;
        grasp(i).time_stamps(j,13)=(str2double(grasp(i).time_stamps(j,4))...
            *20+start_idx(i))/1000;

       combo=[combo; grasp(i).time_stamps(j,:)];
    end    
end
%%
rast.time_stamps(:,:)=combo(:,10:13);
rast.grasp=combo(:,5);
rast.type=combo(:,6);
rast.hit=combo(:,8);
rast.block=combo(:,9);
rast.slc_times=combo(:,1:4);
clear blockNum dat dh first_trial  h hits hold_start i j move_start...
    pre_start rng start start_idx temp_array trial_end

%% SAVE YOUR DATA

post_process.rast=rast;
post_process.nsp=nsp;
post_process.slc=slc;
post_process.path=path;
post_process.blocks.cl=cl_blocks;
post_process.blocks.ol=ol_blocks;
post_process.blocks.all=blocks;
post_process.filt.chan=wchan;
post_process.filt.used=filt_used;
post_process.name=nspName;
post_process.array.sets=arrayChannelSets;
post_process.array.names=arrayNames;
post_process.loop_id=loop_id;

%%
disp('saving post processed data')
save([working_dir filesep 'processed_data.mat'], 'post_process')
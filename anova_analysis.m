function [mult_comp]=anova_analysis(sdf, time)
%does an anova for the first length ms of the premovement and movement
%epochs (cannot exceed 1 second)

time=1000*time; %inputted as seconds and convert to ms
anova_names={'power' 'palmer' 'pinch' 'disc' 'lateral'};
pre_sig_chans=[];
move_sig_chans=[];

for i = 1:128
    anova_pre=[];
    pre_names=[];
    anova_move=[];
    move_names=[];
    for j=1:5
        a=sdf{i,j};
        for k =1:length(a)
           b= a{1,k};
           if isempty(b)
               b=zeros(3000,1);
           end
           h=b(1:1000)';
           p=b(1001:(1001+(time-1)))';
           m=b(2001:(2001+(time-1)))';

           hold=mean(h);
           pre=mean(p);
           move=mean(m);

           anova_pre=[anova_pre, hold, pre];
           pre_names=[pre_names, {'hold'}, {anova_names{j}}];
           anova_move=[anova_move, hold, move];
           move_names=[move_names, {'hold'}, {anova_names{j}}];
        end

    end

    [ppre,~,spre]=anova1(anova_pre,pre_names,'off'); %outputs p, table, and stats
    %[ppre,~,spre]=anova1(anova_pre,pre_names);
    c=multcompare(spre);
    [pmove,~,smove]=anova1(anova_move, move_names,'off');
    %[pmove,~,smove]=anova1(anova_move, move_names);
    d=multcompare(smove);
    mult_comp{i,1}=c;
    mult_comp{i,2}=d;
    mult_comp{i,3}=ppre;
    mult_comp{i,4}=pmove;
    mult_comp{i,5}=anova_pre;
    mult_comp{i,6}=pre_names;
    mult_comp{i,7}=anova_move;
    mult_comp{i,8}=move_names;
    close all
    
    if ppre<.05
        pre_sig_chans=[pre_sig_chans , i];
        comp_data=c;
            comp_groups=[];
            for j = 1:15 %only look at comparisons 
                p_val=comp_data(j,6);
                if p_val<.05
                    comp_groups=[comp_groups, {[num2str(comp_data(j,1)) '/' num2str(comp_data(j,2))]}];
                end
            end
        mult_comp{i,9}=pre_sig_chans;
        mult_comp{i,10}=comp_groups;

    elseif pmove<.05
        move_sig_chans=[move_sig_chans , i];
        comp_data=d;
            comp_groups=[];
            for j = 1:15 %only look at comparisons 
                p_val=comp_data(j,6);
                if p_val<.05
                    comp_groups=[comp_groups, {[num2str(comp_data(j,1)) '/' num2str(comp_data(j,2))]}];
                end
            end
        mult_comp{i,11}=move_sig_chans;
        mult_comp{i,12}=comp_groups;
    else
        mult_comp{i,9}=[];
        mult_comp{i,10}=[];
        mult_comp{i,11}=[];
        mult_comp{i,12}=[];
    end
end

end
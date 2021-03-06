function [sdf,mult_comp, pre_dpca, move_dpca]=full_analysis(ol_blocks, cl_blocks,...
    filt_used, box_folder ,nsp_num, loop_id)
%the full analysis consists of three steps:
%   1. sorting the trials
%   2. generating psths
%   3. sorting through an anova for 
%   4. apply dPCA

[post_process, working_dir] = sort_trials(ol_blocks, cl_blocks,...
    filt_used, box_folder ,nsp_num, loop_id);
close all
[sdf, anova_sdf] = gen_psth(post_process, working_dir);
close all
[mult_comp]=anova_analysis(anova_sdf, .5);
close all
[pre_dpca, move_dpca]=prep_dpca(sdf, working_dir, 1:64);
[pre_dpca, move_dpca]=prep_dpca(sdf, working_dir, 65:128);

end
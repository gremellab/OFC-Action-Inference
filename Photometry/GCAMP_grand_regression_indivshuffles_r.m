function [Model] = GCAMP_grand_regression_indivshuffles_r(Grouped_GCAMP)
%Last updated 5/02/22 by Christian Cazares
%This function takes the Grouped_GCAMP data structure
%and specifically the individual lme table table
%and plops all the individual animal data together into one big table
%which is used to build the simple and complex LMEs predicting activity at
%different time points

tic
Grand_Regression_Cell =[];
Shuf_distributions = [];
Shuf_n_back_Lengths_distribution ={};

%lump all the mice together into one big cell
for i = 1:length(Grouped_GCAMP.Mice)
    Grand_Regression_Cell = [Grand_Regression_Cell; Grouped_GCAMP.Mice{1,i}.regression_cell];
    Shuf_distributions = [Shuf_distributions; Grouped_GCAMP.Mice{1,i}.Shuffled_Durs_Distribtuion ]; 
    Shuf_n_back_Lengths_distribution = [Shuf_n_back_Lengths_distribution Grouped_GCAMP.Mice{1,i}.shuf_n_back_Lengths_distribution];      
end

%get column names from the individual mouse data, since the vars dont
%change
column_names = Grouped_GCAMP.Mice{1,1}.data_table_variables;
Grand_Regression_Table = cell2table(Grand_Regression_Cell, 'VariableNames',column_names);

%add in an n zero reward, for use with post offset graphs
 Grand_Regression_Table.n_zero_reward = Grand_Regression_Table.n_minus_zero_Durations_All>=Grand_Regression_Table.criteria_indicator;

%% shuffled lp order regressions
%create vars for shuffled coefficients

shuf_lme_interp_AUC_lponly_n0_coef_ALL_simple = [];
shuf_lme_interp_AUC_lponly_n1_coef_ALL_simple = [];

shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL_simple = [];
shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_simple = [];

shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_simple = [];
shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_simple = [];


shuf_lme_interp_AUC_lponly_n0_coef_ALL_complex = [];
shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex = [];

shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL_complex = [];
shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex = [];

shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_complex = [];
shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex = [];

%Find Presses over 10s, use this idx to exclude them from the model
longidx =find(Grand_Regression_Table.n_minus_zero_Durations_All>=10000);
%should be redundant, removed long presses earlier in pipeline


for j =  1:size(Shuf_n_back_Lengths_distribution,1)
    %need to erase the i_shuffled array every time, this will be used to create
    %the table for regression
    i_shuffled =[];
    for i = 1:size(Shuf_n_back_Lengths_distribution,2)
        %first shuffle for all mice
        i_shuffled = [i_shuffled; Shuf_n_back_Lengths_distribution{j,i}'];
    end
    
    %add in the appropriate shuffled n lp (ie, the shuffled dist that built the
    %n-back array)
    i_shuffled = [i_shuffled Shuf_distributions(:,j)];
    
    T_shuffled = array2table(i_shuffled,'VariableNames',{'shuf_n_minus_one_Durations_All',...
        'shuf_n_minus_two_Durations_All', 'shuf_n_minus_three_Durations_All', 'shuf_n_minus_four_Durations_All',...
        'shuf_n_minus_five_Durations_All', 'shuf_n_minus_six_Durations_All', 'shuf_n_minus_seven_Durations_All',...
        'shuf_n_minus_eight_Durations_All', 'shuf_n_minus_nine_Durations_All', 'shuf_n_minus_ten_Durations_All',...
        'shuf_LP_Durations_All'});
    
    T_Both = [T_shuffled Grand_Regression_Table];
    
    %% shuffle each individual lag in the interpolated models
    %e.g., only n-0 duration is shuffled
    % Simple
    shuf_modelspec_interp_AUC_lponly_n0_simple = 'interp_all_AUC ~ shuf_LP_Durations_All + interp_all_AUC_n1 +  n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + lp_on_times + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
    shuf_lme_interp_AUC_lponly_n0_simple = fitlme(T_Both,shuf_modelspec_interp_AUC_lponly_n0_simple,'Exclude',longidx);
    shuf_lme_interp_AUC_lponly_n0_coef_simple = shuf_lme_interp_AUC_lponly_n0_simple.Coefficients.Estimate;
    shuf_lme_interp_AUC_lponly_n0_coef_ALL_simple = [shuf_lme_interp_AUC_lponly_n0_coef_ALL_simple shuf_lme_interp_AUC_lponly_n0_coef_simple];
    
    shuf_modelspec_interp_AUC_lponly_n1_simple = 'interp_all_AUC ~ n_minus_zero_Durations_All + interp_all_AUC_n1 + shuf_n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + lp_on_times + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
    shuf_lme_interp_AUC_lponly_n1_simple = fitlme(T_Both,shuf_modelspec_interp_AUC_lponly_n1_simple,'Exclude',longidx);
    shuf_lme_interp_AUC_lponly_n1_coef_simple = shuf_lme_interp_AUC_lponly_n1_simple.Coefficients.Estimate;
    shuf_lme_interp_AUC_lponly_n1_coef_ALL_simple = [shuf_lme_interp_AUC_lponly_n1_coef_ALL_simple shuf_lme_interp_AUC_lponly_n1_coef_simple];
    
    % Complex
    shuf_modelspec_interp_AUC_lponly_n0_complex = 'interp_all_AUC ~ shuf_LP_Durations_All + interp_all_AUC_n1 +  n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + n_minus_one_Durations_All:n_minus_one_reward_All + n_minus_one_Durations_All:HE_n1_Indicator + n_minus_one_Durations_All:ipi1 + lp_on_times + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
    shuf_lme_interp_AUC_lponly_n0_complex = fitlme(T_Both,shuf_modelspec_interp_AUC_lponly_n0_complex,'Exclude',longidx);
    shuf_lme_interp_AUC_lponly_n0_coef_complex = shuf_lme_interp_AUC_lponly_n0_complex.Coefficients.Estimate;
    shuf_lme_interp_AUC_lponly_n0_coef_ALL_complex = [shuf_lme_interp_AUC_lponly_n0_coef_ALL_complex shuf_lme_interp_AUC_lponly_n0_coef_complex];
    
    shuf_modelspec_interp_AUC_lponly_n1_complex = 'interp_all_AUC ~ n_minus_zero_Durations_All + interp_all_AUC_n1 + shuf_n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + shuf_n_minus_one_Durations_All:n_minus_one_reward_All + shuf_n_minus_one_Durations_All:HE_n1_Indicator + shuf_n_minus_one_Durations_All:ipi1 + lp_on_times + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
    shuf_lme_interp_AUC_lponly_n1_complex = fitlme(T_Both,shuf_modelspec_interp_AUC_lponly_n1_complex,'Exclude',longidx);
    shuf_lme_interp_AUC_lponly_n1_coef_complex = shuf_lme_interp_AUC_lponly_n1_complex.Coefficients.Estimate;
    shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex = [shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex shuf_lme_interp_AUC_lponly_n1_coef_complex];

    
    %% shuffle each lag in the pre onset
    % Simple
    shuf_modelspec_pre_onset_base_raw_1s_n0_simple  =  'base_neg_1_to_onset_AUC_raw ~ shuf_LP_Durations_All + base_neg_1_to_onset_AUC_raw_n1 + n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + lp_on_times + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
    shuf_lme_pre_onset_base_raw_1s_n0_simple = fitlme(T_Both,shuf_modelspec_pre_onset_base_raw_1s_n0_simple,'Exclude',longidx);
    shuf_lme_pre_onset_base_raw_1s_n0_coef_simple = shuf_lme_pre_onset_base_raw_1s_n0_simple.Coefficients.Estimate;
    shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL_simple = [shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL_simple shuf_lme_pre_onset_base_raw_1s_n0_coef_simple];
    
    shuf_modelspec_pre_onset_base_raw_1s_n1_simple  = 'base_neg_1_to_onset_AUC_raw ~ n_minus_zero_Durations_All + base_neg_1_to_onset_AUC_raw_n1 + shuf_n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + lp_on_times + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
    shuf_lme_pre_onset_base_raw_1s_n1_simple = fitlme(T_Both,shuf_modelspec_pre_onset_base_raw_1s_n1_simple,'Exclude',longidx);
    shuf_lme_pre_onset_base_raw_1s_n1_coef_simple = shuf_lme_pre_onset_base_raw_1s_n1_simple.Coefficients.Estimate;
    shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_simple = [shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_simple shuf_lme_pre_onset_base_raw_1s_n1_coef_simple];
    
    % Complex
    shuf_modelspec_pre_onset_base_raw_1s_n0_complex  =  'base_neg_1_to_onset_AUC_raw ~ shuf_LP_Durations_All + base_neg_1_to_onset_AUC_raw_n1 + n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + n_minus_one_Durations_All:n_minus_one_reward_All + n_minus_one_Durations_All:HE_n1_Indicator + n_minus_one_Durations_All:ipi1 + lp_on_times + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
    shuf_lme_pre_onset_base_raw_1s_n0_complex = fitlme(T_Both,shuf_modelspec_pre_onset_base_raw_1s_n0_complex,'Exclude',longidx);
    shuf_lme_pre_onset_base_raw_1s_n0_coef_complex = shuf_lme_pre_onset_base_raw_1s_n0_complex.Coefficients.Estimate;
    shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL_complex = [shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL_complex shuf_lme_pre_onset_base_raw_1s_n0_coef_complex];
    
    shuf_modelspec_pre_onset_base_raw_1s_n1_complex  = 'base_neg_1_to_onset_AUC_raw ~ n_minus_zero_Durations_All + base_neg_1_to_onset_AUC_raw_n1 + shuf_n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + shuf_n_minus_one_Durations_All:n_minus_one_reward_All + shuf_n_minus_one_Durations_All:HE_n1_Indicator + shuf_n_minus_one_Durations_All:ipi1 + lp_on_times + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
    shuf_lme_pre_onset_base_raw_1s_n1_complex = fitlme(T_Both,shuf_modelspec_pre_onset_base_raw_1s_n1_complex,'Exclude',longidx);
    shuf_lme_pre_onset_base_raw_1s_n1_coef_complex = shuf_lme_pre_onset_base_raw_1s_n1_complex.Coefficients.Estimate;
    shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex = [shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex shuf_lme_pre_onset_base_raw_1s_n1_coef_complex];
    

    %% post offset 0 to 1s
    % Simple
    shuf_modelspec_post_offset_base_raw_1s_n0_simple =  ' base_offset_to1_AUC_raw ~ shuf_LP_Durations_All + base_offset_to1_AUC_raw_n1 +  n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + lp_on_times + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
    shuf_lme_post_offset_base_raw_1s_n0_simple = fitlme(T_Both,shuf_modelspec_post_offset_base_raw_1s_n0_simple,'Exclude',longidx);
    shuf_lme_post_offset_base_raw_1s_n0_coef_simple = shuf_lme_post_offset_base_raw_1s_n0_simple.Coefficients.Estimate;
    shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_simple = [shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_simple shuf_lme_post_offset_base_raw_1s_n0_coef_simple];
    
    shuf_modelspec_post_offset_base_raw_1s_n1_simple =  ' base_offset_to1_AUC_raw ~  n_minus_zero_Durations_All + base_offset_to1_AUC_raw_n1 +  shuf_n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + lp_on_times + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
    shuf_lme_post_offset_base_raw_1s_n1_simple = fitlme(T_Both,shuf_modelspec_post_offset_base_raw_1s_n1_simple,'Exclude',longidx);
    shuf_lme_post_offset_base_raw_1s_n1_coef_simple = shuf_lme_post_offset_base_raw_1s_n1_simple.Coefficients.Estimate;
    shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_simple = [shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_simple shuf_lme_post_offset_base_raw_1s_n1_coef_simple];
    
    % Complex
    shuf_modelspec_post_offset_base_raw_1s_n0_complex =  ' base_offset_to1_AUC_raw ~ shuf_LP_Durations_All + base_offset_to1_AUC_raw_n1 +  n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + n_minus_one_Durations_All:n_minus_one_reward_All + n_minus_one_Durations_All:HE_n1_Indicator + n_minus_one_Durations_All:ipi1 + lp_on_times + n_zero_reward + shuf_LP_Durations_All*n_zero_reward + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
    shuf_lme_post_offset_base_raw_1s_n0_complex = fitlme(T_Both,shuf_modelspec_post_offset_base_raw_1s_n0_complex,'Exclude',longidx);
    shuf_lme_post_offset_base_raw_1s_n0_coef_complex = shuf_lme_post_offset_base_raw_1s_n0_complex.Coefficients.Estimate;
    shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_complex = [shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_complex shuf_lme_post_offset_base_raw_1s_n0_coef_complex];
    
    shuf_modelspec_post_offset_base_raw_1s_n1_complex =  ' base_offset_to1_AUC_raw ~  n_minus_zero_Durations_All + base_offset_to1_AUC_raw_n1 +  shuf_n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + shuf_n_minus_one_Durations_All:n_minus_one_reward_All + shuf_n_minus_one_Durations_All:HE_n1_Indicator + shuf_n_minus_one_Durations_All:ipi1 + lp_on_times + n_zero_reward + n_minus_zero_Durations_All*n_zero_reward + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
    shuf_lme_post_offset_base_raw_1s_n1_complex = fitlme(T_Both,shuf_modelspec_post_offset_base_raw_1s_n1_complex,'Exclude',longidx);
    shuf_lme_post_offset_base_raw_1s_n1_coef_complex = shuf_lme_post_offset_base_raw_1s_n1_complex.Coefficients.Estimate;
    shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex = [shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex shuf_lme_post_offset_base_raw_1s_n1_coef_complex];


end

%% first create a purely beahvioral LME predicting N duration
%n - 1 to n - 10 onlt
modelspec_behavior_lponly = 'n_minus_zero_Durations_All ~ n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + n_minus_one_Durations_All:n_minus_one_reward_All + n_minus_one_Durations_All:HE_n1_Indicator + n_minus_one_Durations_All:ipi1 + lp_on_times + (1|mouse_indicator) + (1|day_indicator) + (1|criteria_indicator)';
lme_behavior_lponly = fitlme(Grand_Regression_Table,modelspec_behavior_lponly,'Exclude',longidx);
lme_behavior_lponly_coef = lme_behavior_lponly.Coefficients.Estimate;
lme_behavior_lponly_se = lme_behavior_lponly.Coefficients.SE;


%% interpolated activity during the press

%simple relation to durations
modelspec_interp_all_n6_act = 'interp_all_AUC ~ n_minus_zero_Durations_All + interp_all_AUC_n1 + n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + lp_on_times + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
lme_interp_all_n6_act = fitlme(Grand_Regression_Table,modelspec_interp_all_n6_act,'Exclude',longidx);
Model.Interp.lme_interp_all_n6_act = lme_interp_all_n6_act;
lme_interp_all_n6_act_coef = lme_interp_all_n6_act.Coefficients.Estimate;
Model.Interp.lme_interp_all_n6_act_coef = lme_interp_all_n6_act_coef;
Model.Interp.lme_interp_all_n6_act_se = lme_interp_all_n6_act.Coefficients.SE;
Model.Interp.lme_interp_all_n6_act_pval = lme_interp_all_n6_act.Coefficients.pValue;
Model.Interp.lme_interp_all_n6_act_name = lme_interp_all_n6_act.Coefficients.Name;
lme_interp_all_n6_act_anova = anova(lme_interp_all_n6_act);
Model.Interp.lme_interp_all_n6_act_anova = lme_interp_all_n6_act_anova;
Model.Interp.lme_interp_all_n6_act_fstat = lme_interp_all_n6_act_anova.FStat;
Model.Interp.lme_interp_all_n6_act_fpval = lme_interp_all_n6_act_anova.pValue;
Model.Interp.lme_interp_all_n6_act_upper = lme_interp_all_n6_act.Coefficients.Upper;
Model.Interp.lme_interp_all_n6_act_lower = lme_interp_all_n6_act.Coefficients.Lower;

%perm test comparing to shuffled coefs
n_p_interp_n0_simple = sum(abs(shuf_lme_interp_AUC_lponly_n0_coef_ALL_simple(find(contains(shuf_lme_interp_AUC_lponly_n0_simple.Coefficients.Name,'shuf_LP_Durations_All')),:)) > abs(lme_interp_all_n6_act_coef(find(contains(lme_interp_all_n6_act.Coefficients.Name,'n_minus_zero_Durations_All')))))/size(shuf_lme_interp_AUC_lponly_n0_coef_ALL_simple,2);
n_p_interp_n1_simple = sum(abs(shuf_lme_interp_AUC_lponly_n1_coef_ALL_simple(find(contains(shuf_lme_interp_AUC_lponly_n1_simple.Coefficients.Name,'shuf_n_minus_one_Durations_All')),:)) > abs(lme_interp_all_n6_act_coef(find(contains(lme_interp_all_n6_act.Coefficients.Name,'n_minus_one_Durations_All')))))/size(shuf_lme_interp_AUC_lponly_n1_coef_ALL_simple,2);

n_p_interp_n0_to_n6_simple = [n_p_interp_n0_simple; n_p_interp_n1_simple];

%shuffled mean/sem
interp_n0_shuf_mean = mean(shuf_lme_interp_AUC_lponly_n0_coef_ALL_simple(find(contains(shuf_lme_interp_AUC_lponly_n0_simple.Coefficients.Name,'shuf_LP_Durations_All')),:));
interp_n0_shuf_sem = std(shuf_lme_interp_AUC_lponly_n0_coef_ALL_simple(find(contains(shuf_lme_interp_AUC_lponly_n0_simple.Coefficients.Name,'shuf_LP_Durations_All')),:))/sqrt(size(shuf_lme_interp_AUC_lponly_n0_coef_ALL_simple,2));
interp_n1_shuf_mean = mean(shuf_lme_interp_AUC_lponly_n1_coef_ALL_simple(find(contains(shuf_lme_interp_AUC_lponly_n1_simple.Coefficients.Name,'shuf_n_minus_one_Durations_All')),:));
interp_n1_shuf_sem = std(shuf_lme_interp_AUC_lponly_n1_coef_ALL_simple(find(contains(shuf_lme_interp_AUC_lponly_n1_simple.Coefficients.Name,'shuf_n_minus_one_Durations_All')),:))/sqrt(size(shuf_lme_interp_AUC_lponly_n1_coef_ALL_simple,2));

interp_n0_n6_mean_simple = [interp_n0_shuf_mean; interp_n1_shuf_mean];
interp_n0_n6_sem_simple = [interp_n0_shuf_sem; interp_n1_shuf_sem];

Model.Interp.Shuff_Coeff_Name_Simple = [{'n-0 Duration'}; {'n-1 Duration'}];
Model.Interp.Shuff_pval_Simple = n_p_interp_n0_to_n6_simple;
Model.Interp.Shuff_Coeff_Simple = interp_n0_n6_mean_simple;
Model.Interp.Shuff_SE_Simple = interp_n0_n6_sem_simple;

%now add in the cannonical beh model terms in addition
modelspec_interp_all_n6_act_andints = 'interp_all_AUC ~ n_minus_zero_Durations_All + interp_all_AUC_n1 + n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + n_minus_one_Durations_All:n_minus_one_reward_All + n_minus_one_Durations_All:HE_n1_Indicator + n_minus_one_Durations_All:ipi1 + lp_on_times + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
lme_interp_all_n6_act_andints = fitlme(Grand_Regression_Table,modelspec_interp_all_n6_act_andints,'Exclude',longidx);
Model.Interp.lme_interp_all_n6_act_andints = lme_interp_all_n6_act_andints;
Model.Interp.lme_interp_all_n6_act_andints_names = lme_interp_all_n6_act_andints.Coefficients.Name;
lme_interp_all_n6_act_andints_coef = lme_interp_all_n6_act_andints.Coefficients.Estimate;
Model.Interp.lme_interp_all_n6_act_andints_coef = lme_interp_all_n6_act_andints_coef;
Model.Interp.lme_interp_all_n6_act_andints_pvals = lme_interp_all_n6_act_andints.Coefficients.pValue;	
Model.Interp.lme_interp_all_n6_act_andints_se = lme_interp_all_n6_act_andints.Coefficients.SE;	
Model.Interp.lme_interp_all_n6_act_andints_df = lme_interp_all_n6_act_andints.Coefficients.DF;	
Model.Interp.lme_interp_all_n6_act_andints_upper = lme_interp_all_n6_act_andints.Coefficients.Upper;	
Model.Interp.lme_interp_all_n6_act_andints_lower = lme_interp_all_n6_act_andints.Coefficients.Lower;	
lme_interp_all_n6_act_andints_anova = anova(lme_interp_all_n6_act_andints);
Model.Interp.lme_interp_all_n6_act_andints_anova = lme_interp_all_n6_act_andints_anova;
Model.Interp.lme_interp_all_n6_act_andints_anova_Fval = lme_interp_all_n6_act_andints_anova.FStat;
Model.Interp.lme_interp_all_n6_act_andints_anova_Fpval = lme_interp_all_n6_act_andints_anova.pValue;


%perm test comparing to shuffled coefs
n_p_interp_n0_complex = sum(abs(shuf_lme_interp_AUC_lponly_n0_coef_ALL_complex(find(contains(shuf_lme_interp_AUC_lponly_n0_complex.Coefficients.Name,'shuf_LP_Durations_All')),:)) > abs(lme_interp_all_n6_act_andints_coef(find(contains(lme_interp_all_n6_act_andints.Coefficients.Name,'n_minus_zero_Durations_All')))))/size(shuf_lme_interp_AUC_lponly_n0_coef_ALL_complex,2);
n_p_interp_n1_complex = sum(abs(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex(2,:)) > abs(lme_interp_all_n6_act_andints_coef(3)))/size(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex,2);
n_p_interp_n1_outcome_complex = sum(abs(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex(find(contains(shuf_lme_interp_AUC_lponly_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:n_minus_one_reward_All_true')),:)) > abs(lme_interp_all_n6_act_andints_coef(find(contains(lme_interp_all_n6_act_andints.Coefficients.Name,'n_minus_one_Durations_All:n_minus_one_reward_All_true')))))/size(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex,2);
n_p_interp_n1_he_complex =  sum(abs(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex(find(contains(shuf_lme_interp_AUC_lponly_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:HE_n1_Indicator_1')),:)) > abs(lme_interp_all_n6_act_andints_coef(find(contains(lme_interp_all_n6_act_andints.Coefficients.Name,'n_minus_one_Durations_All:HE_n1_Indicator_1')))))/size(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex,2);
n_p_interp_n1_ipi_complex = sum(abs(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex(find(contains(shuf_lme_interp_AUC_lponly_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:ipi1')),:)) > abs(lme_interp_all_n6_act_andints_coef(find(contains(lme_interp_all_n6_act_andints.Coefficients.Name,'n_minus_one_Durations_All:ipi1')))))/size(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex,2);

n_p_interp_n0_to_n6_complex = [n_p_interp_n0_complex; n_p_interp_n1_complex; n_p_interp_n1_outcome_complex; n_p_interp_n1_he_complex; n_p_interp_n1_ipi_complex];

%shuffled mean/sem
interp_n0_shuf_mean = mean(shuf_lme_interp_AUC_lponly_n0_coef_ALL_complex(find(contains(shuf_lme_interp_AUC_lponly_n0_complex.Coefficients.Name,'shuf_LP_Durations_All')),:));
interp_n0_shuf_sem = std(shuf_lme_interp_AUC_lponly_n0_coef_ALL_complex(find(contains(shuf_lme_interp_AUC_lponly_n0_complex.Coefficients.Name,'shuf_LP_Durations_All')),:))/sqrt(length(shuf_lme_interp_AUC_lponly_n0_coef_ALL_complex));
interp_n1_shuf_mean = mean(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex(2,:));
interp_n1_shuf_sem = std(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex(2,:))/sqrt(length(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex));

interp_n1_outcome_shuf_mean = mean(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex(find(contains(shuf_lme_interp_AUC_lponly_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:n_minus_one_reward_All_true')),:));
interp_n1_outcome_shuf_sem = std(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex(find(contains(shuf_lme_interp_AUC_lponly_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:n_minus_one_reward_All_true')),:))/sqrt(size(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex,2));
interp_n1_he_shuf_mean = mean(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex(find(contains(shuf_lme_interp_AUC_lponly_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:HE_n1_Indicator_1')),:));
interp_n1_he_shuf_sem = std(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex(find(contains(shuf_lme_interp_AUC_lponly_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:HE_n1_Indicator_1')),:))/sqrt(size(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex,2));
interp_n1_ipi_shuf_mean = mean(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex(find(contains(shuf_lme_interp_AUC_lponly_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:ipi1')),:));
interp_n1_ipi_shuf_sem = std(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex(find(contains(shuf_lme_interp_AUC_lponly_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:ipi1')),:))/sqrt(size(shuf_lme_interp_AUC_lponly_n1_coef_ALL_complex,2));


interp_n0_n6_mean_complex = [interp_n0_shuf_mean; interp_n1_shuf_mean; interp_n1_outcome_shuf_mean; interp_n1_he_shuf_mean; interp_n1_ipi_shuf_mean];
interp_n0_n6_sem_complex = [interp_n0_shuf_sem; interp_n1_shuf_sem; interp_n1_outcome_shuf_sem; interp_n1_he_shuf_sem; interp_n1_ipi_shuf_sem];

Model.Interp.Shuff_Coeff_Name_Complex = [{'n-0 Duration'}; {'n-1 Duration'}; {'n-1 Duration * n-1 outcome'}; {'n-1 Duration * n-1 head entry'}; {'n-1 Duration * n-1 ipi'}];
Model.Interp.Shuff_pval_Complex = n_p_interp_n0_to_n6_complex;
Model.Interp.Shuff_Coeff_Complex = interp_n0_n6_mean_complex;
Model.Interp.Shuff_SE_Complex = interp_n0_n6_sem_complex;

%How well do the simple/complex models predict activity?
[yhat yhatCI] =predict(lme_interp_all_n6_act,Grand_Regression_Table,'Simultaneous',true);
correctish_pred =Grand_Regression_Table.interp_all_AUC <=yhatCI(:,2) &  Grand_Regression_Table.interp_all_AUC >=yhatCI(:,1);
correctCI_prop = sum(correctish_pred)/length(Grand_Regression_Table.interp_all_AUC)
Model.Interp.simple_correctish_pred = correctish_pred;
Model.Interp.simple_correctCI_prop = correctCI_prop;

 
[yhat yhatCI] =predict(lme_interp_all_n6_act_andints,Grand_Regression_Table,'Simultaneous',true);
correctish_pred = Grand_Regression_Table.interp_all_AUC <=yhatCI(:,2) &  Grand_Regression_Table.interp_all_AUC >=yhatCI(:,1);
correctCI_prop = sum(correctish_pred)/length(Grand_Regression_Table.interp_all_AUC)
Model.Interp.complex_correctish_pred = correctish_pred;
Model.Interp.complex_correctCI_prop = correctCI_prop;


%% pre onset 1s 
% Simple
modelspec_pre_onset_base_raw_1s_n6_act = 'base_neg_1_to_onset_AUC_raw ~ n_minus_zero_Durations_All + base_neg_1_to_onset_AUC_raw_n1 + n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + lp_on_times + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
lme_pre_onset_base_raw_1s_n6_act = fitlme(Grand_Regression_Table,modelspec_pre_onset_base_raw_1s_n6_act,'Exclude',longidx);
Model.Onset.lme_pre_onset_base_raw_1s_n6_act = lme_pre_onset_base_raw_1s_n6_act;
lme_pre_onset_base_raw_1s_n6_act_coef = lme_pre_onset_base_raw_1s_n6_act.Coefficients.Estimate;
Model.Onset.lme_pre_onset_base_raw_1s_n6_act_coef = lme_pre_onset_base_raw_1s_n6_act_coef;
Model.Onset.lme_pre_onset_base_raw_1s_n6_act_se = lme_pre_onset_base_raw_1s_n6_act.Coefficients.SE;
Model.Onset.lme_pre_onset_base_raw_1s_n6_act_pval = lme_pre_onset_base_raw_1s_n6_act.Coefficients.pValue;
Model.Onset.lme_pre_onset_base_raw_1s_n6_act_name = lme_pre_onset_base_raw_1s_n6_act.Coefficients.Name;
lme_pre_onset_base_raw_1s_n6_act_anova = anova(lme_pre_onset_base_raw_1s_n6_act);
Model.Onset.lme_pre_onset_base_raw_1s_n6_act_anova = lme_pre_onset_base_raw_1s_n6_act_anova;
Model.Onset.lme_pre_onset_base_raw_1s_n6_act_fstat = lme_pre_onset_base_raw_1s_n6_act_anova.FStat;
Model.Onset.lme_pre_onset_base_raw_1s_n6_act_fpval = lme_pre_onset_base_raw_1s_n6_act_anova.pValue;
Model.Onset.lme_pre_onset_base_raw_1s_n6_act_upper = lme_pre_onset_base_raw_1s_n6_act.Coefficients.Upper;
Model.Onset.lme_pre_onset_base_raw_1s_n6_act_lower = lme_pre_onset_base_raw_1s_n6_act.Coefficients.Lower;

avg_preonset_AUC = nanmean(Grand_Regression_Table.base_neg_1_to_onset_AUC_raw);
med_preonset_AUC = nanmedian(Grand_Regression_Table.base_neg_1_to_onset_AUC_raw);


%perm test comparing to shuffled coefs
n_p_pre_n0_simple = sum( abs(shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL_simple(find(contains(shuf_lme_pre_onset_base_raw_1s_n0_simple.Coefficients.Name,'shuf_LP_Durations_All')),:)) > abs(lme_pre_onset_base_raw_1s_n6_act_coef(find(contains(lme_pre_onset_base_raw_1s_n6_act.Coefficients.Name,'n_minus_zero_Durations_All')))))/size(shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL_simple,2);
n_p_pre_n1_simple = sum( abs(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_simple(find(contains(shuf_lme_pre_onset_base_raw_1s_n1_simple.Coefficients.Name,'shuf_n_minus_one_Durations_All')),:)) > abs(lme_pre_onset_base_raw_1s_n6_act_coef(find(contains(lme_pre_onset_base_raw_1s_n6_act.Coefficients.Name,'n_minus_one_Durations_All')))))/size(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_simple,2);

n_p_pre_n0_to_n6_simple = [n_p_pre_n0_simple; n_p_pre_n1_simple];

%shuffled mean/sem
pre_n0_shuf_mean = mean(shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL_simple(find(contains(shuf_lme_pre_onset_base_raw_1s_n0_simple.Coefficients.Name,'shuf_LP_Durations_All')),:));
pre_n0_shuf_sem = std(shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL_simple(find(contains(shuf_lme_pre_onset_base_raw_1s_n0_simple.Coefficients.Name,'shuf_LP_Durations_All')),:))/sqrt(length(shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL_simple));
pre_n1_shuf_mean = mean(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_simple(find(contains(shuf_lme_pre_onset_base_raw_1s_n1_simple.Coefficients.Name,'shuf_n_minus_one_Durations_All')),:));
pre_n1_shuf_sem = std(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_simple(find(contains(shuf_lme_pre_onset_base_raw_1s_n1_simple.Coefficients.Name,'shuf_n_minus_one_Durations_All')),:))/sqrt(length(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_simple));

pre_n0_n6_mean_simple = [pre_n0_shuf_mean; pre_n1_shuf_mean];
pre_n0_n6_sem_simple = [pre_n0_shuf_sem; pre_n1_shuf_sem];

Model.Onset.Shuff_Coeff_Name_Simple = [{'n-0 Duration'}; {'n-1 Duration'}];
Model.Onset.Shuff_pval_Simple = n_p_pre_n0_to_n6_simple;
Model.Onset.Shuff_Coeff_Simple = pre_n0_n6_mean_simple;
Model.Onset.Shuff_SE_Simple = pre_n0_n6_sem_simple;


% Complex
%pre onset 1s n1n6 ints
modelspec_preonset_all_n6_act_andints = 'base_neg_1_to_onset_AUC_raw ~ n_minus_zero_Durations_All + base_neg_1_to_onset_AUC_raw_n1 + n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + n_minus_one_Durations_All:n_minus_one_reward_All + n_minus_one_Durations_All:HE_n1_Indicator + n_minus_one_Durations_All:ipi1 + lp_on_times + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
lme_preonset_all_n6_act_andints = fitlme(Grand_Regression_Table,modelspec_preonset_all_n6_act_andints,'Exclude',longidx);
Model.Onset.lme_preonset_all_n6_act_andints = lme_preonset_all_n6_act_andints;
lme_preonset_all_n6_act_andints_coef = lme_preonset_all_n6_act_andints.Coefficients.Estimate;
Model.Onset.lme_preonset_all_n6_act_andints_coef = lme_preonset_all_n6_act_andints_coef;
Model.Onset.lme_preonset_all_n6_act_andints_se = lme_preonset_all_n6_act_andints.Coefficients.SE;
Model.Onset.lme_preonset_all_n6_act_andints_pvals = lme_preonset_all_n6_act_andints.Coefficients.pValue;	
Model.Onset.lme_preonset_all_n6_act_andints_names = lme_preonset_all_n6_act_andints.Coefficients.Name;
Model.Onset.lme_preonset_all_n6_act_andints_df = lme_preonset_all_n6_act_andints.Coefficients.DF;	
Model.Onset.lme_preonset_all_n6_act_andints_upper = lme_preonset_all_n6_act_andints.Coefficients.Upper;	
Model.Onset.lme_preonset_all_n6_act_andints_lower = lme_preonset_all_n6_act_andints.Coefficients.Lower;	
lme_preonset_all_n6_act_andints_anova = anova(lme_preonset_all_n6_act_andints);
Model.Onset.lme_preonset_all_n6_act_andints_anova = lme_preonset_all_n6_act_andints_anova;
Model.Onset.lme_preonset_all_n6_act_andints_anova_Fval = lme_preonset_all_n6_act_andints_anova.FStat;
Model.Onset.lme_preonset_all_n6_act_andints_anova_Fpval = lme_preonset_all_n6_act_andints_anova.pValue;

%perm test comparing to shuffled coefs
n_p_pre_n0_complex = sum(abs(shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL_complex(find(contains(shuf_lme_pre_onset_base_raw_1s_n0_complex.Coefficients.Name,'shuf_LP_Durations_All')),:)) > abs(lme_preonset_all_n6_act_andints_coef(find(contains(lme_preonset_all_n6_act_andints.Coefficients.Name,'n_minus_zero_Durations_All')))))/size(shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL_complex,2);
n_p_pre_n1_complex = sum(abs(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex(2,:)) > abs(lme_preonset_all_n6_act_andints_coef(3)))/size(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex,2);
n_p_pre_n1_outcome_complex = sum(abs(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_pre_onset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:n_minus_one_reward_All_true')),:)) > abs(lme_preonset_all_n6_act_andints_coef(find(contains(lme_preonset_all_n6_act_andints.Coefficients.Name,'n_minus_one_Durations_All:n_minus_one_reward_All_true')))))/size(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex,2);
n_p_pre_n1_he_complex =  sum(abs(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_pre_onset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:HE_n1_Indicator_1')),:)) > abs(lme_preonset_all_n6_act_andints_coef(find(contains(lme_preonset_all_n6_act_andints.Coefficients.Name,'n_minus_one_Durations_All:HE_n1_Indicator_1')))))/size(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex,2);
n_p_pre_n1_ipi_complex = sum(abs(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_pre_onset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:ipi1')),:)) > abs(lme_preonset_all_n6_act_andints_coef(find(contains(lme_preonset_all_n6_act_andints.Coefficients.Name,'n_minus_one_Durations_All:ipi1')))))/size(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex,2);

n_p_pre_n0_to_n6_complex = [n_p_pre_n0_complex; n_p_pre_n1_complex; n_p_pre_n1_outcome_complex; n_p_pre_n1_he_complex; n_p_pre_n1_ipi_complex];

%shuffled mean/sem
pre_n0_shuf_mean = mean(shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL_complex(find(contains(shuf_lme_pre_onset_base_raw_1s_n0_complex.Coefficients.Name,'shuf_LP_Durations_All')),:));
pre_n0_shuf_sem = std(shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL_complex(find(contains(shuf_lme_pre_onset_base_raw_1s_n0_complex.Coefficients.Name,'shuf_LP_Durations_All')),:))/sqrt(size((shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL_complex),2));
pre_n1_shuf_mean = mean(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex(2,:));
pre_n1_shuf_sem = std(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex(2,:))/sqrt(size((shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex),2));

pre_n1_outcome_shuf_mean = mean(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_pre_onset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:n_minus_one_reward_All_true')),:));
pre_n1_outcome_shuf_sem = std(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_pre_onset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:n_minus_one_reward_All_true')),:))/sqrt(size((shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex),2));
pre_n1_he_shuf_mean = mean(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_pre_onset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:HE_n1_Indicator_1')),:));
pre_n1_he_shuf_sem = std(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_pre_onset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:HE_n1_Indicator_1')),:))/sqrt(size((shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex),2));
pre_n1_ipi_shuf_mean = mean(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_pre_onset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:ipi1')),:));
pre_n1_ipi_shuf_sem = std(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_pre_onset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:ipi1')),:))/sqrt(size((shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL_complex),2));


pre_n0_n6_mean_complex = [pre_n0_shuf_mean; pre_n1_shuf_mean; pre_n1_outcome_shuf_mean; pre_n1_he_shuf_mean; pre_n1_ipi_shuf_mean];
pre_n0_n6_sem_complex = [pre_n0_shuf_sem; pre_n1_shuf_sem; pre_n1_outcome_shuf_sem; pre_n1_he_shuf_sem; pre_n1_ipi_shuf_sem];

Model.Onset.Shuff_Coeff_Name_Complex = [{'n-0 Duration'}; {'n-1 Duration'}; {'n-1 Duration * n-1 outcome'}; {'n-1 Duration * n-1 head entry'}; {'n-1 Duration * n-1 ipi'}];
Model.Onset.Shuff_pval_Complex = n_p_pre_n0_to_n6_complex;
Model.Onset.Shuff_Coeff_Complex = pre_n0_n6_mean_complex;
Model.Onset.Shuff_SE_Complex = pre_n0_n6_sem_complex;


%How well do the simple/complex models predict activity?
[yhat yhatCI] =predict(lme_pre_onset_base_raw_1s_n6_act,Grand_Regression_Table,'Simultaneous',true);
correctish_pred =Grand_Regression_Table.base_neg_1_to_onset_AUC_raw <=yhatCI(:,2) &  Grand_Regression_Table.base_neg_1_to_onset_AUC_raw >=yhatCI(:,1);
correctCI_prop = sum(correctish_pred)/length(Grand_Regression_Table.base_neg_1_to_onset_AUC_raw)
Model.Onset.simple_correctish_pred = correctish_pred;
Model.Onset.simple_correctCI_prop = correctCI_prop;

 
[yhat yhatCI] =predict(lme_preonset_all_n6_act_andints,Grand_Regression_Table,'Simultaneous',true);
correctish_pred =Grand_Regression_Table.base_neg_1_to_onset_AUC_raw <=yhatCI(:,2) &  Grand_Regression_Table.base_neg_1_to_onset_AUC_raw >=yhatCI(:,1);
correctCI_prop = sum(correctish_pred)/length(Grand_Regression_Table.base_neg_1_to_onset_AUC_raw)
Model.Onset.complex_correctish_pred = correctish_pred;
Model.Onset.complex_correctCI_prop = correctCI_prop;

%% post offset 1s
modelspec_post_offset_base_raw_to1_n6_act = 'base_offset_to1_AUC_raw ~ n_minus_zero_Durations_All + base_offset_to1_AUC_raw_n1 + n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + lp_on_times + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
lme_post_offset_base_raw_to1_n6_act = fitlme(Grand_Regression_Table,modelspec_post_offset_base_raw_to1_n6_act,'Exclude',longidx);
Model.Offset.lme_post_offset_base_raw_to1_n6_act = lme_post_offset_base_raw_to1_n6_act;
Model.Offset.lme_post_offset_base_raw_to1_n6_act_name = lme_post_offset_base_raw_to1_n6_act.Coefficients.Name;
lme_post_offset_base_raw_to1_n6_act_coef = lme_post_offset_base_raw_to1_n6_act.Coefficients.Estimate;
Model.Offset.lme_post_offset_base_raw_to1_n6_act_coef = lme_post_offset_base_raw_to1_n6_act_coef;
Model.Offset.lme_post_offset_base_raw_to1_n6_act_se = lme_post_offset_base_raw_to1_n6_act.Coefficients.SE;
Model.Offset.lme_post_offset_base_raw_to1_n6_act_pval = lme_post_offset_base_raw_to1_n6_act.Coefficients.pValue;
Model.Offset.lme_post_offset_base_raw_to1_n6_act_name = lme_post_offset_base_raw_to1_n6_act.Coefficients.Name;
lme_post_offset_base_raw_to1_n6_act_anova = anova(lme_post_offset_base_raw_to1_n6_act);
Model.Offset.lme_post_offset_base_raw_to1_n6_act_anova = lme_post_offset_base_raw_to1_n6_act_anova;
Model.Offset.lme_post_offset_base_raw_to1_n6_act_fstat = lme_post_offset_base_raw_to1_n6_act_anova.FStat;
Model.Offset.lme_post_offset_base_raw_to1_n6_act_fpval = lme_post_offset_base_raw_to1_n6_act_anova.pValue;
Model.Offset.lme_post_offset_base_raw_to1_n6_act_upper = lme_post_offset_base_raw_to1_n6_act.Coefficients.Upper;
Model.Offset.lme_post_offset_base_raw_to1_n6_act_lower = lme_post_offset_base_raw_to1_n6_act.Coefficients.Lower;

avg_post_offset_AUC = nanmean(Grand_Regression_Table.base_offset_to1_AUC_raw);
med_post_offset_AUC = nanmedian(Grand_Regression_Table.base_offset_to1_AUC_raw);

%perm test comparing to shuffled coefs
n_p_post_n0_simple = sum( abs(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_simple(find(contains(shuf_lme_post_offset_base_raw_1s_n0_simple.Coefficients.Name,'shuf_LP_Durations_All')),:)) > abs(lme_post_offset_base_raw_to1_n6_act_coef(find(contains(lme_pre_onset_base_raw_1s_n6_act.Coefficients.Name,'n_minus_zero_Durations_All')))))/size(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_simple,2);
n_p_post_n1_simple = sum( abs(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_simple(find(contains(shuf_lme_post_offset_base_raw_1s_n1_simple.Coefficients.Name,'shuf_n_minus_one_Durations_All')),:)) > abs(lme_post_offset_base_raw_to1_n6_act_coef(find(contains(lme_pre_onset_base_raw_1s_n6_act.Coefficients.Name,'n_minus_one_Durations_All')))))/size(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_simple,2);

n_p_post_n0_to_n6_simple = [n_p_post_n0_simple; n_p_post_n1_simple];

%shuffled mean/sem
post_n0_shuf_mean = mean(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_simple(find(contains(shuf_lme_post_offset_base_raw_1s_n0_simple.Coefficients.Name,'shuf_LP_Durations_All')),:));
post_n0_shuf_sem = std(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_simple(find(contains(shuf_lme_post_offset_base_raw_1s_n0_simple.Coefficients.Name,'shuf_LP_Durations_All')),:))/sqrt(length(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_simple));
post_n1_shuf_mean = mean(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_simple(find(contains(shuf_lme_post_offset_base_raw_1s_n1_simple.Coefficients.Name,'shuf_n_minus_one_Durations_All')),:));
post_n1_shuf_sem = std(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_simple(find(contains(shuf_lme_post_offset_base_raw_1s_n1_simple.Coefficients.Name,'shuf_n_minus_one_Durations_All')),:))/sqrt(length(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_simple));

post_n0_n6_mean_simple = [post_n0_shuf_mean; post_n1_shuf_mean];
post_n0_n6_sem_simple = [post_n0_shuf_sem; post_n1_shuf_sem];

Model.Offset.Shuff_Coeff_Name_Simple = [{'n-0 Duration'}; {'n-1 Duration'}];
Model.Offset.Shuff_pval_Simple = n_p_post_n0_to_n6_simple;
Model.Offset.Shuff_Coeff_Simple = post_n0_n6_mean_simple;
Model.Offset.Shuff_SE_Simple = post_n0_n6_sem_simple;

% add in n_zero_reward (i.e., was the just finished lp rewarded) and
% interaction
% n_minus_zero_Durations_All*n_zero_reward +
% moving_average_lp_length_n7andback:n_zero_reward
modelspec_postoffset_all_n6_act_andints_n0 = 'base_offset_to1_AUC_raw ~ n_minus_zero_Durations_All + base_offset_to1_AUC_raw + n_minus_one_Durations_All + HE_n1_Indicator + ipi1 + n_minus_one_reward_All + n_minus_one_Durations_All:n_minus_one_reward_All + n_minus_one_Durations_All:HE_n1_Indicator + n_minus_one_Durations_All:ipi1 + lp_on_times + n_zero_reward + n_minus_zero_Durations_All*n_zero_reward + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
lme_postoffset_all_n6_act_andints_n0 = fitlme(Grand_Regression_Table,modelspec_postoffset_all_n6_act_andints_n0,'Exclude',longidx);
Model.Offset.lme_postoffset_all_n6_act_andints_n0 = lme_postoffset_all_n6_act_andints_n0;
lme_postoffset_all_n6_act_andints_n0_coef = lme_postoffset_all_n6_act_andints_n0.Coefficients.Estimate;
Model.Offset.lme_postoffset_all_n6_act_andints_n0_coef = lme_postoffset_all_n6_act_andints_n0_coef;
Model.Offset.lme_postoffset_all_n6_act_andints_n0_pvals = lme_postoffset_all_n6_act_andints_n0.Coefficients.pValue;	
Model.Offset.lme_postoffset_all_n6_act_andints_n0_se = lme_postoffset_all_n6_act_andints_n0.Coefficients.SE;	
Model.Offset.lme_postoffset_all_n6_act_andints_n0_names = lme_postoffset_all_n6_act_andints_n0.Coefficients.Name;
Model.Offset.lme_postoffset_all_n6_act_andints_n0_df = lme_postoffset_all_n6_act_andints_n0.Coefficients.DF;	
Model.Offset.lme_postoffset_all_n6_act_andints_n0_upper = lme_postoffset_all_n6_act_andints_n0.Coefficients.Upper;	
Model.Offset.lme_postoffset_all_n6_act_andints_n0_lower = lme_postoffset_all_n6_act_andints_n0.Coefficients.Lower;	
lme_postoffset_all_n6_act_andints_n0_anova = anova(lme_postoffset_all_n6_act_andints_n0);
Model.Offset.lme_postoffset_all_n6_act_andints_n0_anova = lme_postoffset_all_n6_act_andints_n0_anova ;
Model.Offset.lme_postoffset_all_n6_act_andints_n0_anova_Fval = lme_postoffset_all_n6_act_andints_n0_anova.FStat;
Model.Offset.lme_postoffset_all_n6_act_andints_n0_anova_Fpval = lme_postoffset_all_n6_act_andints_n0_anova.pValue;


%perm test comparing to shuffled coefs
n_p_post_n0_complex = sum( abs(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_complex(2,:)) > abs(lme_postoffset_all_n6_act_andints_n0_coef(2)))/size(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_complex,2);
n_p_post_n1_complex = sum( abs(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex(2,:)) > abs(lme_postoffset_all_n6_act_andints_n0_coef(3)))/size(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex,2);

n_p_pre_n0_reward_complex = sum( abs(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_complex(find(contains(shuf_lme_post_offset_base_raw_1s_n0_complex.Coefficients.Name,'shuf_LP_Durations_All:n_zero_reward_1')),:)) > abs(lme_postoffset_all_n6_act_andints_n0_coef(find(contains(lme_postoffset_all_n6_act_andints_n0.Coefficients.Name,'n_minus_zero_Durations_All:n_zero_reward_1')))))/size(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_complex,2);
n_p_pre_n1_outcome_complex = sum(abs(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_post_offset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:n_minus_one_reward_All_true')),:)) > abs(lme_postoffset_all_n6_act_andints_n0_coef(find(contains(lme_postoffset_all_n6_act_andints_n0.Coefficients.Name,'n_minus_one_Durations_All:n_minus_one_reward_All_true')))))/size(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex,2);
n_p_pre_n1_he_complex =  sum(abs(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_post_offset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:HE_n1_Indicator_1')),:)) > abs(lme_postoffset_all_n6_act_andints_n0_coef(find(contains(lme_postoffset_all_n6_act_andints_n0.Coefficients.Name,'n_minus_one_Durations_All:HE_n1_Indicator_1')))))/size(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex,2);
n_p_pre_n1_ipi_complex = sum(abs(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_post_offset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:ipi1')),:)) > abs(lme_postoffset_all_n6_act_andints_n0_coef(find(contains(lme_postoffset_all_n6_act_andints_n0.Coefficients.Name,'n_minus_one_Durations_All:ipi1')))))/size(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex,2);

n_p_post_n0_to_n6_complex = [n_p_post_n0_complex; n_p_post_n1_complex; n_p_pre_n0_reward_complex; n_p_pre_n1_outcome_complex; n_p_pre_n1_he_complex; n_p_pre_n1_ipi_complex];

%shuffled mean/sem
post_n0_shuf_mean = mean(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_complex(2,:));
post_n0_shuf_sem = std(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_complex(2,:))/sqrt(size((shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_complex),2));
post_n1_shuf_mean = mean(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex(2,:));
post_n1_shuf_sem = std(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex(2,:))/sqrt(size((shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex),2));

post_n0_n0_reward_shuf_mean = mean(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_complex(find(contains(shuf_lme_post_offset_base_raw_1s_n0_complex.Coefficients.Name,'shuf_LP_Durations_All:n_zero_reward_1')),:));
post_n0_n0_reward_shuf_sem = std(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_complex(find(contains(shuf_lme_post_offset_base_raw_1s_n0_complex.Coefficients.Name,'shuf_LP_Durations_All:n_zero_reward_1')),:))/sqrt(size((shuf_lme_post_offset_base_raw_1s_n0_coef_ALL_complex),2));

post_n1_outcome_shuf_mean = mean(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_post_offset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:n_minus_one_reward_All_true')),:));
post_n1_outcome_shuf_sem = std(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_post_offset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:n_minus_one_reward_All_true')),:))/sqrt(size((shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex),2));
post_n1_he_shuf_mean = mean(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_post_offset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:HE_n1_Indicator_1')),:));
post_n1_he_shuf_sem = std(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_post_offset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:HE_n1_Indicator_1')),:))/sqrt(size((shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex),2));
post_n1_ipi_shuf_mean = mean(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_post_offset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:ipi1')),:));
post_n1_ipi_shuf_sem = std(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex(find(contains(shuf_lme_post_offset_base_raw_1s_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:ipi1')),:))/sqrt(size((shuf_lme_post_offset_base_raw_1s_n1_coef_ALL_complex),2));

post_n0_n6_mean_complex = [post_n0_shuf_mean; post_n1_shuf_mean; post_n0_n0_reward_shuf_mean; post_n1_outcome_shuf_mean; post_n1_he_shuf_mean; post_n1_ipi_shuf_mean];
post_n0_n6_sem_complex =[post_n0_shuf_sem; post_n1_shuf_sem; post_n0_n0_reward_shuf_sem; post_n1_outcome_shuf_sem; post_n1_he_shuf_sem; post_n1_ipi_shuf_sem];

Model.Offset.Shuff_Coeff_Name_complex = [{'n-0 Duration'}; {'n-1 Duration'}; {'n-0 Duration * n-0 outcome'}; {'n-1 Duration * n-1 outcome'}; {'n-1 Duration * n-1 head entry'}; {'n-1 Duration * n-1 ipi'}];
Model.Offset.Shuff_pval_complex = n_p_post_n0_to_n6_complex;
Model.Offset.Shuff_Coeff_complex = post_n0_n6_mean_complex;
Model.Offset.Shuff_SE_complex = post_n0_n6_sem_complex;


%How well do the simple/complex models predict activity?
[yhat yhatCI] =predict(lme_post_offset_base_raw_to1_n6_act,Grand_Regression_Table,'Simultaneous',true);
correctish_pred =Grand_Regression_Table.base_offset_to1_AUC_raw <=yhatCI(:,2) &  Grand_Regression_Table.base_offset_to1_AUC_raw >=yhatCI(:,1);
correctCI_prop = sum(correctish_pred)/length(Grand_Regression_Table.base_offset_to1_AUC_raw)
Model.Offset.simple_correctish_pred = correctish_pred;
Model.Offset.simple_correctCI_prop = correctCI_prop;

 
[yhat yhatCI] =predict(lme_postoffset_all_n6_act_andints_n0,Grand_Regression_Table,'Simultaneous',true);
correctish_pred =Grand_Regression_Table.base_offset_to1_AUC_raw <=yhatCI(:,2) &  Grand_Regression_Table.base_offset_to1_AUC_raw >=yhatCI(:,1);
correctCI_prop = sum(correctish_pred)/length(Grand_Regression_Table.base_offset_to1_AUC_raw)
Model.Offset.complex_correctish_pred = correctish_pred;
Model.Offset.complex_correctCI_prop = correctCI_prop;

toc
end


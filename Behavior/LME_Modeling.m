%% Purpose
%last edited 5/02/22 by Christian Cazares
%The main purpose of this script is to create Linear Mixed Effect Models (LMEs)
%to predict the duration of lever presses (n) given subjective experiential
%information. This includes things like n-back press durations, n-back
%headentries (HE), n-back reward, interpressinterval (IPI), time in
%session, % of presses over criteria.
%These are mixed effect models becuase they will also include random
%effects of day of training and mouse, to account for the repeated nature
%of the datasets. There are sections here for creating relatively simple
%LMEs, as well as more complex LMEs that are arrived at via a priori
%hypothesis/questions of interest combined with model selection.

%an example simple model would seek to predict n duration given prior press
%durations, whether those presses were rewarded, and interactions between
%those terms
%a = regression coeffecient
%E = error term
% C =  press duration
%R = reward or no
%Cn = ao +  + a1C(n-1) + a2C(n-2) + a3C(n-3)

%           + a4R(n-1) + a5R(n-2) + a6R(n-3)
%           + a7C(n-1)R(n-1) + a8C(n-2)R(n-2) +
%            a9C(n-3)R(n-3) + E(t)
 
%R will be coded as 0 for failed press, or 1 for rewarded
%C will be duration in ms

%There is also some other code for smaller questions that will still
%utilize the Data structure obtained from the
%MEDPC_Behavior_Extract_For_Regressions_r
%This includes things like generating within session performance graphs,
%histograms of press durations, and raw duration
%values before, during, or after optogenetic stimulation (if applicable).

%For this, you need to have a Data file generated from MEDPC_Behavior_Extract_For_Regressions
%with this Data structure being loaded under "load" (you must be in the
%folder with the data structure).
%The Data File is Organized First by Day of Training: Data.Day
%The Next level is mouse: Data.Day(x).Mouse
%Within mouse, we have a large structure that includes all the variables for
%all the mice on a given day. Some of these are single numbers (e.g., total number of lever presses)
%while others are cell arrays (e.g. all the individual lever presses a
%mouse made).
%This table goes in "mouse order". So whichever mouse was first
%(numerically) in the original folder will be first here. 
%EACH ROW IS ONE MOUSE.
%The Data.Day(x).Mouse.Name variable will label each mouse.


%% Load the Data and create variables
tic
nback_iter = 3; % how many lever presses back to look at
total_mice = size(Data.Day(2).Mouse,2);
%get number of days
day_number = size(Data.Day,2);
%Initialize Variables
LP_Durations_Cell = {};
LP_Timestamps_Cell ={};
Rewarded_Durations_Cell = {};
Unrewarded_Durations_Cell ={};
IPI_Cell ={};
Shuffled_LP_Duration_Cell ={};
HE_Indicator_Cell ={};
mouse_indicator = [];
day_indicator=[];
ID_indicator=[];
Logcial_LP_Cell ={};
cell_length_var = 0;
Shufffled_Logical_Lever_Press_Cell = {};
HE_n_1_Indicator_Cell ={};
HE_n_2_Indicator_Cell = {};
HE_n_3_Indicator_Cell = {};
IPI_nback_Cell ={};
criteria_indicator_Cell={};
criteria_percent_indicator_Cell={};
up_state_idx_n1_Cell ={};
up_state_idx_n2_Cell ={};
up_state_idx_n3_Cell ={};
avg_duration_indicator_Cell ={};
reward_indicator_Cell ={};
reward_indicator_n1_Cell ={};
reward_indicator_n2_Cell ={};
reward_indicator_n3_Cell ={};
%used with opto stim
stim_indicator_cell = {};
stim_indicator_n1_cell = {};
stim_indicator_n2_cell = {};
stim_indicator_n3_cell = {};
histogram_counts_Cell ={};
histogram_bins_Cell={};
histogam_indexs_Cell={};
%indiv lme predictions
indiv_correctCIprop_criteria_Cell ={};
indiv_correctCIprop_criteria_complex_Cell ={};
indiv_lme_10_indiv_se_Cell ={};
indiv_lme_10_indiv_coef_Cell ={};
indiv_lme_reduced_se_Cell ={};
indiv_lme_reduced_coef_Cell ={};
%individual lme r^2s
r_squared_criteria_Cell={};
r_squared_adjusted_criteria_Cell = {};
r_squared_complex_criteria_Cell = {};
r_squared_adjusted_complex_criteria_Cell = {};
mouse_number_for_corrs =[];
day_number_for_corrs =[];
%% Loop through data structure to get data input into Cells
for i = 1:day_number %outer loop goes across days  
    for j = 1:total_mice %inner loop goes across mice within a day        
    
     IPI_Cell = [IPI_Cell Data.Day(i).Mouse(j).Lever_Press_IPI];
     mouse_number_for_corrs = [mouse_number_for_corrs; j];
     day_number_for_corrs = [day_number_for_corrs; i];
     indiv_lme_10_indiv_se_Cell = [indiv_lme_10_indiv_se_Cell Data.Day(i).Mouse(j).lme_10_indiv_se];
     indiv_lme_10_indiv_coef_Cell = [indiv_lme_10_indiv_coef_Cell   Data.Day(i).Mouse(j).lme_10_indiv_coef];
     indiv_lme_reduced_se_Cell =[indiv_lme_reduced_se_Cell Data.Day(i).Mouse(j).lme_moderate_remove_se];
     indiv_lme_reduced_coef_Cell =[indiv_lme_reduced_coef_Cell Data.Day(i).Mouse(j).lme_moderate_remove_coef];
     indiv_correctCIprop_criteria_Cell =[indiv_correctCIprop_criteria_Cell Data.Day(i).Mouse(j).correctCIprop_criteria];
     indiv_correctCIprop_criteria_complex_Cell =[indiv_correctCIprop_criteria_complex_Cell Data.Day(i).Mouse(j).correctCIprop_criteria_complex];
     r_squared_criteria_Cell=[r_squared_criteria_Cell  Data.Day(i).Mouse(j).r_squared_criteria];
     r_squared_adjusted_criteria_Cell = [r_squared_adjusted_criteria_Cell Data.Day(i).Mouse(j).r_squared_adjusted_criteria];
     r_squared_complex_criteria_Cell = [r_squared_complex_criteria_Cell Data.Day(i).Mouse(j).r_squared_complex_criteria ];
     r_squared_adjusted_complex_criteria_Cell = [r_squared_adjusted_complex_criteria_Cell Data.Day(i).Mouse(j).r_squared_adjusted_complex_criteria ];  
     LP_Durations_Cell = [LP_Durations_Cell Data.Day(i).Mouse(j).Lever_Press_duration];
     Logcial_LP_Cell = [Logcial_LP_Cell Data.Day(i).Mouse(j).Logical_LPs];
     LP_Timestamps_Cell = [LP_Timestamps_Cell Data.Day(i).Mouse(j).Lever_Press_ts];
     Rewarded_Durations_Cell = [Rewarded_Durations_Cell Data.Day(i).Mouse(j).Lever_Press_Rewarded_Lengths];
     Unrewarded_Durations_Cell = [Unrewarded_Durations_Cell Data.Day(i).Mouse(j).Lever_Press_Not_Rewarded_Lengths];
     %create indicator variables that will give us an array to mark every LP
     %with the mouse and day it happened. These will be the Random Effects
     cell_length_var = cell_length_var + 1;
     %we multiple the ones here by j, the loop iteration variable for mice.
     %thus mice, won't have the same numbers in the folders, but will
     %instead be labeled 1:total_mice based on their order in the Data
     %Structure. Ditto for day_indicator, where we multiple an array of
     %ones by i (the day iteration variable). Of note, these indicators are
     %made to be as long as the Lever press durations, so that EACH LP is
     %tagged both with mouse and day.
     mouse_indicator = [mouse_indicator; (ones(size(LP_Durations_Cell{:,cell_length_var}))*j)];
     ID = Data.Day(i).Mouse(j).Name(1:4);
     C    = cell(size(LP_Durations_Cell{:,cell_length_var}));
     C(:) = {ID};
     ID_indicator = [ID_indicator; C];
     day_indicator = [day_indicator;     (ones(size(LP_Durations_Cell{:,cell_length_var}))*i)];
     %these indicators already have the appropriate size
     criteria_indicator_Cell = [criteria_indicator_Cell Data.Day(i).Mouse(j).criteria_indicator];
     criteria_percent_indicator_Cell = [criteria_percent_indicator_Cell  Data.Day(i).Mouse(j).crit_percent_indicator];
     avg_duration_indicator_Cell = [avg_duration_indicator_Cell  Data.Day(i).Mouse(j).avg_duration_indicator];
     reward_indicator_Cell =[reward_indicator_Cell Data.Day(i).Mouse(j).RE_Indicator];
     reward_indicator_n1_Cell =[reward_indicator_n1_Cell Data.Day(i).Mouse(j).RE_n_1_Indicator];
     reward_indicator_n2_Cell =[reward_indicator_n2_Cell Data.Day(i).Mouse(j).RE_n_2_Indicator];
     reward_indicator_n3_Cell =[reward_indicator_n3_Cell Data.Day(i).Mouse(j).RE_n_3_Indicator];
     stim_indicator_cell = [stim_indicator_cell Data.Day(i).Mouse(j).stim_ind_all];
     stim_indicator_n1_cell =[stim_indicator_n1_cell Data.Day(i).Mouse(j).stim_ind_all_n1];
     stim_indicator_n2_cell = [stim_indicator_n2_cell Data.Day(i).Mouse(j).stim_ind_all_n2 ];
     stim_indicator_n3_cell  = [stim_indicator_n3_cell Data.Day(i).Mouse(j).stim_ind_all_n3 ];
     Shuffled_LP_Duration_Cell =[Shuffled_LP_Duration_Cell Data.Day(i).Mouse(j).Shuffled_Lever_Press_lengths];
     HE_Indicator_Cell = [HE_Indicator_Cell Data.Day(i).Mouse(j).HE_Indicator];
     Shufffled_Logical_Lever_Press_Cell =[Shufffled_Logical_Lever_Press_Cell Data.Day(i).Mouse(j).Shuffled_Logical_Lever_press];
     HE_n_1_Indicator_Cell = [HE_n_1_Indicator_Cell Data.Day(i).Mouse(j).HE_n_1_Indicator];
     HE_n_2_Indicator_Cell = [HE_n_2_Indicator_Cell Data.Day(i).Mouse(j).HE_n_2_Indicator];
     HE_n_3_Indicator_Cell = [HE_n_3_Indicator_Cell Data.Day(i).Mouse(j).HE_n_3_Indicator];
     histogram_counts_Cell =[histogram_counts_Cell  Data.Day(i).Mouse(j).histogram_counts];
     histogram_bins_Cell=[histogram_bins_Cell Data.Day(i).Mouse(j).histogram_edges];
     histogam_indexs_Cell=[histogam_indexs_Cell Data.Day(i).Mouse(j).histogram_bin_indexs];    
    end %mouse loop end
end %day loop end

%% loop through the cell to plop all animals together in a single matrix, rather than a cell.
%we will use this long matrices to create a table for LMEs
LP_Durations_All = [];
Logical_LP_All =[];
LP_Timestamps_All =[];
Rewarded_Durations_All = [];
Unrewarded_Durations_All =[];
n_minus_one_Rewarded_Durations_All = [];
n_minus_two_Rewarded_Durations_All =[];
n_minus_three_Rewarded_Durations_All =[];
n_minus_four_Rewarded_Durations_All =[];
IPI_All =[];
Shuffled_Durations_All =[];
Shufffled_Logical_Lever_Press_All =[];
HE_Indicator_All =[];
HE_n_1_Indicator_All = []; 
HE_n_2_Indicator_All = []; 
HE_n_3_Indicator_All = [];
criteria_indicator_All =[];
criteria_percent_indicator_All =[];
avg_duration_indicator_All =[];
reward_indicator_All =[];
reward_indicator_n1_All =[];
reward_indicator_n2_All =[];
reward_indicator_n3_All =[];
stim_indicator_All=[];
stim_indicator_n1_All=[];
stim_indicator_n2_All=[];
stim_indicator_n3_All=[];
histogram_counts_All =[];
histogram_bins_All=[];
histogam_indexs_All=[];       
%individual LME predictions
correctCIprop_criteria_All =[];
correctCIprop_criteria_complex_All =[];
indiv_lme_10_indiv_se_All =[];
indiv_lme_10_indiv_coef_All =[];
indiv_lme_reduced_se_All =[];
indiv_lme_reduced_coef_All =[];
r_squared_criteria_All=[];
r_squared_adjusted_criteria_All = [];
r_squared_complex_criteria_All = [];
r_squared_adjusted_complex_criteria_All = [];
for i = 1:length(LP_Durations_Cell)    
     correctCIprop_criteria_All =[correctCIprop_criteria_All; indiv_correctCIprop_criteria_Cell{i}];
     correctCIprop_criteria_complex_All =[correctCIprop_criteria_complex_All; indiv_correctCIprop_criteria_complex_Cell{i}];
     indiv_lme_10_indiv_se_All =[indiv_lme_10_indiv_se_All indiv_lme_10_indiv_se_Cell{i}];
     indiv_lme_10_indiv_coef_All =[indiv_lme_10_indiv_coef_All indiv_lme_10_indiv_coef_Cell{i}];
     indiv_lme_reduced_se_All =[indiv_lme_reduced_se_All indiv_lme_reduced_se_Cell{i}];
     indiv_lme_reduced_coef_All =[indiv_lme_reduced_coef_All indiv_lme_reduced_coef_Cell{i}];
     r_squared_criteria_All=[r_squared_criteria_All; r_squared_criteria_Cell{i}];
     r_squared_adjusted_criteria_All = [r_squared_adjusted_criteria_All; r_squared_adjusted_criteria_Cell{i}];
     r_squared_complex_criteria_All = [r_squared_complex_criteria_All; r_squared_complex_criteria_Cell{i}];
     r_squared_adjusted_complex_criteria_All = [r_squared_adjusted_complex_criteria_All; r_squared_adjusted_complex_criteria_Cell{i}];
     LP_Durations_All = [LP_Durations_All; LP_Durations_Cell{i}]; 
     Logical_LP_All = [Logical_LP_All; Logcial_LP_Cell{i}];
     LP_Timestamps_All = [LP_Timestamps_All; LP_Timestamps_Cell{i}];
     Rewarded_Durations_All = [Rewarded_Durations_All; Rewarded_Durations_Cell{i}];
     Unrewarded_Durations_All = [Unrewarded_Durations_All; Unrewarded_Durations_Cell{i}];
     IPI_All =[IPI_All; IPI_Cell{i}];
     Shuffled_Durations_All = [Shuffled_Durations_All; Shuffled_LP_Duration_Cell{i}];
     Shufffled_Logical_Lever_Press_All =[Shuffled_Durations_All; Shufffled_Logical_Lever_Press_Cell{i}];
     HE_Indicator_All = [HE_Indicator_All; HE_Indicator_Cell{i}];
     HE_n_1_Indicator_All =[HE_n_1_Indicator_All; HE_n_1_Indicator_Cell{i}]; 
     HE_n_2_Indicator_All = [HE_n_2_Indicator_All; HE_n_2_Indicator_Cell{i}]; 
     HE_n_3_Indicator_All = [HE_n_3_Indicator_All; HE_n_3_Indicator_Cell{i}];
     criteria_indicator_All =[criteria_indicator_All;criteria_indicator_Cell{i}];
     criteria_percent_indicator_All =[criteria_percent_indicator_All;criteria_percent_indicator_Cell{i}];
     avg_duration_indicator_All = [avg_duration_indicator_All; avg_duration_indicator_Cell{i}];
     stim_indicator_All = [stim_indicator_All; stim_indicator_cell{i}];
     stim_indicator_n1_All = [stim_indicator_n1_All; stim_indicator_n1_cell{i}];
     stim_indicator_n2_All = [stim_indicator_n2_All; stim_indicator_n2_cell{i}];
     stim_indicator_n3_All = [stim_indicator_n3_All; stim_indicator_n3_cell{i}];
     reward_indicator_All = [reward_indicator_All; reward_indicator_Cell{i}];
     reward_indicator_n1_All =[reward_indicator_n1_All; reward_indicator_n1_Cell{i}];
     reward_indicator_n2_All =[reward_indicator_n2_All; reward_indicator_n2_Cell{i}];
     reward_indicator_n3_All =[reward_indicator_n3_All; reward_indicator_n3_Cell{i}];
     histogram_counts_All =[histogram_counts_All;  histogram_counts_Cell{i}];
     histogram_bins_All=[histogram_bins_All; histogram_bins_Cell{i}];
     histogam_indexs_All=[histogam_indexs_All; histogam_indexs_Cell{i}];      
end
  
%% the n-back variables are strucuted a bit differently, so loop through
%those on their own
press_n_back_Durations_Cell ={};
logical_n_back_cell ={};
ipi_n_back_Durations_Cell ={};
for i = 1:day_number
     ipi_n_back_Durations_Cell = [ipi_n_back_Durations_Cell Data.Day(i).Mouse(1:total_mice).n_back_IPIs];
     press_n_back_Durations_Cell = [press_n_back_Durations_Cell Data.Day(i).Mouse(1:total_mice).n_back_Lengths];
     logical_n_back_cell = [logical_n_back_cell Data.Day(i).Mouse(1:total_mice).Logical_n_back_Lengths];
end
% Create variables for n-back durations, outcomes, ipis, etc.
n_back_All = {};
logical_n_back_All ={};
ipi_n_back_All ={};
    for i = 1:nback_iter
    n_back_All = [n_back_All; press_n_back_Durations_Cell{i,1:end}];
    logical_n_back_All =[logical_n_back_All; logical_n_back_cell{i,1:end}];
    ipi_n_back_All =[ipi_n_back_All; ipi_n_back_Durations_Cell{i,1:end}];
    end
%end result is a matrix where rows are n back durations, columns are
%lever presses (lined up so n-1 and n-2 are from the same n lever press
%across all animals/days) 
   n_back_All_matrix = cell2mat( n_back_All);
   %check to make sure the cell array is the same length for all nbacks.
   %Need to remove it if animals made fewer than 10 lps, as this will cause
   %the cell to have unequal size and be unable to convert. 
   check_length = length(logical_n_back_All{1,1});
   for check_length_it = 1:length(logical_n_back_All)
       if length(logical_n_back_All{check_length_it}) > check_length_it
           logical_n_back_All{check_length_it} =  logical_n_back_All{check_length_it}(1:check_length);
       end
   end           
logical_n_back_All_matrix = cell2mat(logical_n_back_All);
ipi_n_back_matrix = cell2mat(ipi_n_back_All);
%add in press n duration to the n-back matrix
LP_n_back_DataRaster = [n_back_All_matrix; LP_Durations_All';];
%and press n reward logical
logical_LP_n_back_DataRaster = [logical_n_back_All_matrix; Logical_LP_All';];
%regressions require variables in a table, with each column being a
%variable, so flip
LP_n_back_DataRaster = LP_n_back_DataRaster';
logical_LP_n_back_DataRaster = logical_LP_n_back_DataRaster';
ipi_n_back_matrix = ipi_n_back_matrix';
%% Shuffled order Variables
%we need to do the same thing as above, but for shuffled versions of the
%variables
%this is slightly different, since we will have as many versions of the
%shuffled variables as we had shuffle iterations in MEDPC_Behavior_Extract_For_Regressions
%the cells will be number of mouse*days long, and each cell will have an
%individual animal's data, including all shuffle iterations
shuf_n_back_Lengths_distribution_Cell ={};
Shuffled_Durs_Distribution_Cell ={};
shuf_Logical_n_back_Lengths_distribution_Cell ={};
shuffled_ipi1_Cell = {};
shuffled_ipi2_Cell = {};
shuffled_ipi3_Cell = {};
shuffled_he1_Cell ={};
shuffled_he2_Cell ={};
shuffled_he3_Cell ={};
shuffled_he4_Cell ={};
shuffled_MA_Cell ={};
shuffled_ts_Cell ={};
shuffled_upstate_n1_Cell ={};
shuffled_total_reward_indicator_Cell={};
shuffled_MA_n7_Distribution_Cell={};
Shuffled_stim_ind_all_Cell = {};
Shuffled_stim_ind_all_n1_Cell={};
Shuffled_stim_ind_all_n2_Cell={};
Shuffled_stim_ind_all_n3_Cell={};
Shuffled_stim_ind_all_n4_Cell = {};
Shuffled_stim_ind_all_n5_Cell = {};
Shuffled_stim_ind_all_n6_Cell = {};
for i = 1:day_number
    shuf_n_back_Lengths_distribution_Cell = [shuf_n_back_Lengths_distribution_Cell Data.Day(i).Mouse(1:total_mice).shuf_n_back_Lengths_distribution];
    shuf_Logical_n_back_Lengths_distribution_Cell = [shuf_Logical_n_back_Lengths_distribution_Cell Data.Day(i).Mouse(1:total_mice).shuf_Logical_n_back_Lengths_distribution];
    Shuffled_Durs_Distribution_Cell = [Shuffled_Durs_Distribution_Cell Data.Day(i).Mouse(1:total_mice).Shuffled_Durs_Distribtuion];
    shuffled_ipi1_Cell = [shuffled_ipi1_Cell Data.Day(i).Mouse(1:total_mice).shuffled_ipi1_Distribution];
    shuffled_ipi2_Cell = [shuffled_ipi2_Cell  Data.Day(i).Mouse(1:total_mice).shuffled_ipi2_Distribution];
    shuffled_ipi3_Cell = [shuffled_ipi3_Cell  Data.Day(i).Mouse(1:total_mice).shuffled_ipi3_Distribution];
    shuffled_he1_Cell = [shuffled_he1_Cell  Data.Day(i).Mouse(1:total_mice).shuffled_he1_Distribution];
    shuffled_he2_Cell = [shuffled_he2_Cell  Data.Day(i).Mouse(1:total_mice).shuffled_he2_Distribution];
    shuffled_he3_Cell = [shuffled_he3_Cell  Data.Day(i).Mouse(1:total_mice).shuffled_he3_Distribution];
    shuffled_ts_Cell = [shuffled_ts_Cell  Data.Day(i).Mouse(1:total_mice).Shuffled_ts_Distribution];
    shuffled_total_reward_indicator_Cell=[shuffled_total_reward_indicator_Cell  Data.Day(i).Mouse(1:total_mice).shuffled_total_reward_indicator_dist];
    Shuffled_stim_ind_all_Cell = [Shuffled_stim_ind_all_Cell Data.Day(i).Mouse(1:total_mice).Shuffled_stim_ind_all_Distribution];
    Shuffled_stim_ind_all_n1_Cell = [Shuffled_stim_ind_all_n1_Cell Data.Day(i).Mouse(1:total_mice).Shuffled_stim_ind_all_n1_Distribution];
    Shuffled_stim_ind_all_n2_Cell = [Shuffled_stim_ind_all_n2_Cell Data.Day(i).Mouse(1:total_mice).Shuffled_stim_ind_all_n2_Distribution];
    Shuffled_stim_ind_all_n3_Cell = [Shuffled_stim_ind_all_n3_Cell Data.Day(i).Mouse(1:total_mice).Shuffled_stim_ind_all_n3_Distribution];
end
%now moving to the "all" versions of the cell, columns will be invidiual
%shuffle iterations, while all the mouse/day data will be lumped together
%in the rows, just as in the non-shuffled variables
shuffled_Durs_Dist_All = [];
shuffled_ipi1_All = [];
shuffled_ipi2_All = [];
shuffled_ipi3_All = [];
shuffled_he1_All =[];
shuffled_he2_All =[];
shuffled_he3_All =[];
shuffled_ts_All =[];
shuffled_total_reward_indicator_All=[];
Shuffled_stim_ind_all_All = [];
Shuffled_stim_ind_all_n1_All = [];
Shuffled_stim_ind_all_n2_All = [];
Shuffled_stim_ind_all_n3_All = [];
for i = 1:size(Shuffled_Durs_Distribution_Cell,2)
    shuffled_Durs_Dist_All = [shuffled_Durs_Dist_All; Shuffled_Durs_Distribution_Cell{1,i}]; 
    shuffled_ipi1_All = [shuffled_ipi1_All; shuffled_ipi1_Cell{1,i}];
    shuffled_ipi2_All = [shuffled_ipi2_All; shuffled_ipi2_Cell{1,i}];
    shuffled_ipi3_All = [shuffled_ipi3_All; shuffled_ipi3_Cell{1,i}];
    shuffled_he1_All = [shuffled_he1_All; shuffled_he1_Cell{1,i}];
    shuffled_he2_All = [shuffled_he2_All; shuffled_he2_Cell{1,i}];
    shuffled_he3_All = [shuffled_he3_All; shuffled_he3_Cell{1,i}];
    shuffled_ts_All = [shuffled_ts_All; shuffled_ts_Cell{1,i}];
    shuffled_total_reward_indicator_All=[shuffled_total_reward_indicator_All; shuffled_total_reward_indicator_Cell{1,i}];
    Shuffled_stim_ind_all_All = [Shuffled_stim_ind_all_All; Shuffled_stim_ind_all_Cell{1,i}];
    Shuffled_stim_ind_all_n1_All = [Shuffled_stim_ind_all_n1_All; Shuffled_stim_ind_all_n1_Cell{1,i}];
    Shuffled_stim_ind_all_n2_All = [Shuffled_stim_ind_all_n2_All; Shuffled_stim_ind_all_n2_Cell{1,i}];
    Shuffled_stim_ind_all_n3_All = [Shuffled_stim_ind_all_n3_All; Shuffled_stim_ind_all_n3_Cell{1,i}];
end
%% Now we need to create a table that includes all our variables,
%each variable to a column, rows are all the data together across mice/days
% logical for rewards
logical_LP_10_back_DataRaster = logical_LP_n_back_DataRaster(:,1:nback_iter);
logical_LP_10_back_DataRaster = [logical_LP_10_back_DataRaster logical_LP_n_back_DataRaster(:,nback_iter+1)];
%add in the durations as well
Logical_and_Durations  = [logical_LP_10_back_DataRaster LP_n_back_DataRaster(:,1:nback_iter)];
Logical_and_Durations = [Logical_and_Durations LP_n_back_DataRaster(:,nback_iter+1)];
%convert to a table, name all the variables. 
T10_Logical_and_Continuous = array2table(Logical_and_Durations,'VariableNames',{'n_minus_one_All',...
    'n_minus_two_All', 'n_minus_three_All','LP_All','n_minus_one_Durations_All',...
    'n_minus_two_Durations_All', 'n_minus_three_Durations_All', 'LP_Durations_All'});
%add in indicator variables
T10_Logical_and_Continuous.mouse_indicator =mouse_indicator; %each lp tagged with mouse number
T10_Logical_and_Continuous.day_indicator = day_indicator; %each lp tagged with session number
%make sure logical variables are coded as categorical
T10_Logical_and_Continuous.mouse_indicator = categorical(T10_Logical_and_Continuous.mouse_indicator);
T10_Logical_and_Continuous.day_indicator = categorical(T10_Logical_and_Continuous.day_indicator);
%indicator of mouse's avg. criteria% on a given day
T10_Logical_and_Continuous.criteria_percent_indicator = criteria_percent_indicator_All;
%Did mice make a HE on n-back presses?
T10_Logical_and_Continuous.HE_n_1_Indicator_All = HE_n_1_Indicator_All;
T10_Logical_and_Continuous.HE_n_1_Indicator_All = categorical(T10_Logical_and_Continuous.HE_n_1_Indicator_All);
T10_Logical_and_Continuous.HE_n_2_Indicator_All = HE_n_2_Indicator_All;
T10_Logical_and_Continuous.HE_n_2_Indicator_All = categorical(T10_Logical_and_Continuous.HE_n_2_Indicator_All);
T10_Logical_and_Continuous.HE_n_3_Indicator_All = HE_n_3_Indicator_All;
T10_Logical_and_Continuous.HE_n_3_Indicator_All = categorical(T10_Logical_and_Continuous.HE_n_3_Indicator_All);
%timestamps of when LPs happened 
T10_Logical_and_Continuous.LP_Timestamps_All = LP_Timestamps_All;
%add in IPI between n and various n-backs
T10_Logical_and_Continuous.ipi1 = ipi_n_back_matrix(:,1);
T10_Logical_and_Continuous.ipi2 = ipi_n_back_matrix(:,2);
T10_Logical_and_Continuous.ipi3 = ipi_n_back_matrix(:,3);
%avg duration on a day
T10_Logical_and_Continuous.avg_duration_indicator = avg_duration_indicator_All;
%LP_All is logical for reward, and then n-back logicals
T10_Logical_and_Continuous.LP_All =categorical(T10_Logical_and_Continuous.LP_All);
T10_Logical_and_Continuous.n_minus_one_All =categorical(T10_Logical_and_Continuous.n_minus_one_All);
T10_Logical_and_Continuous.n_minus_two_All =categorical(T10_Logical_and_Continuous.n_minus_two_All);
T10_Logical_and_Continuous.n_minus_three_All =categorical(T10_Logical_and_Continuous.n_minus_three_All);
%indicator for what the hold down criteria was (e.g., 800 or 1600)
T10_Logical_and_Continuous.criteria_indicator_All = criteria_indicator_All;
T10_Logical_and_Continuous.criteria_indicator_All = categorical(T10_Logical_and_Continuous.criteria_indicator_All);
% logical indicators for if stimulation occurred
T10_Logical_and_Continuous.stim_indicator_All = stim_indicator_All;
T10_Logical_and_Continuous.stim_indicator_n1_All = stim_indicator_n1_All;
T10_Logical_and_Continuous.stim_indicator_n2_All = stim_indicator_n2_All;
T10_Logical_and_Continuous.stim_indicator_n3_All = stim_indicator_n3_All;
%logical indicators of Reward - these only differ from the n_minus_x_All
%variables when reward is probabalistic.
T10_Logical_and_Continuous.reward_indicator_All = reward_indicator_All;
T10_Logical_and_Continuous.reward_indicator_n1_All =reward_indicator_n1_All;
T10_Logical_and_Continuous.reward_indicator_n2_All =reward_indicator_n2_All; 
T10_Logical_and_Continuous.reward_indicator_n3_All =reward_indicator_n3_All;
T10_Logical_and_Continuous.reward_indicator_n1_All = categorical(T10_Logical_and_Continuous.reward_indicator_n1_All);
T10_Logical_and_Continuous.reward_indicator_n2_All = categorical(T10_Logical_and_Continuous.reward_indicator_n2_All);
T10_Logical_and_Continuous.reward_indicator_n3_All = categorical(T10_Logical_and_Continuous.reward_indicator_n3_All);
 %% This section is primarily applicable for optogenetic Stimulation, and is used to find variables relative to stimulation
 %e.g., press duration during/after stimulation to look for direct effect
 %also some calculations of overall mouse mean durations and such
 %get lps and nbacks aligned with stim logical (1 = stim, 0 = no stim)
 no_stim_idx =  find(T10_Logical_and_Continuous.stim_indicator_All(:) ~=1);
 stim_idx =find(T10_Logical_and_Continuous.stim_indicator_All(:)==1);
 %press index after stim/nostim
 no_stim_idx_n1 =  find(T10_Logical_and_Continuous.stim_indicator_n1_All(:) ~=1);
 stim_idx_n1 =find(T10_Logical_and_Continuous.stim_indicator_n1_All(:)==1); 
 %press durations 
 lp_stim = T10_Logical_and_Continuous.LP_Durations_All(stim_idx);
 lp_no_stim = T10_Logical_and_Continuous.LP_Durations_All(no_stim_idx);
 %n + 1 press durations
 lp_stim_n1 = T10_Logical_and_Continuous.LP_Durations_All(stim_idx_n1);
 lp_no_stim_n1 = T10_Logical_and_Continuous.LP_Durations_All(no_stim_idx_n1);
 %ipi after stim
 ipi_stim_n1 = T10_Logical_and_Continuous.ipi1(stim_idx_n1);
 ipi_no_stim_n1 = T10_Logical_and_Continuous.ipi1(no_stim_idx_n1);
 %Get Proportion of HEs after stim or after nostim
 he_stim_n1 = T10_Logical_and_Continuous.HE_n_1_Indicator_All(stim_idx_n1);
 he_stim_n1=double(he_stim_n1);
 he_no_stim_n1 = T10_Logical_and_Continuous.HE_n_1_Indicator_All(no_stim_idx_n1);
 he_no_stim_n1=double(he_no_stim_n1);
%get means per animal
%briefly convert categorical mouse/day indicators to double to use indexing
T10_Logical_and_Continuous.day_indicator = double(T10_Logical_and_Continuous.day_indicator);
T10_Logical_and_Continuous.mouse_indicator =double(T10_Logical_and_Continuous.mouse_indicator);
%initialize variables
mouse_means =[];
mouse_std =[];
mouse_mean_over_std =[];
mouse_med=[];
mouse_iqr=[];
mouse_iqr_over_med=[];
mouse_no_stim_mean =[];
mouse_stim_mean =[];
mouse_stim_n1_mean =[];
mouse_no_stim_n1_mean =[];
mouse_stim_n1_ipi_mean =[];
mouse_no_stim_n1_ipi_mean =[];
mouse_stim_n1_std =[];
mouse_no_stim_n1_std =[];
mouse_stim_n_and_n1 =[];
mouse_no_stim_n_and_n1=[];
%create day and mouse iteration variables to loop through
it_d =length(unique(T10_Logical_and_Continuous.day_indicator));
it_m = length(unique(T10_Logical_and_Continuous.mouse_indicator));
for day_mean = 1:it_d
    %loop through days uisng da_index
    day_idx = find(T10_Logical_and_Continuous.day_indicator ==day_mean);
    for mouse_mean = 1:it_m
        %then within each day, loop through mice, using the index
        mouse_idx = find(T10_Logical_and_Continuous.mouse_indicator ==mouse_mean);
        %the intersection of mouse/day will give you data for individual mouse
        %sessions
        mouse_by_day_idx = intersect(day_idx,mouse_idx);
        %now calculate mean, std, med, interquartile range (IQR) overall per
        %mouse/day
        mouse_means = [mouse_means mean(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx))];
        mouse_std = [mouse_std std(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx))];
        mouse_mean_over_std =[mouse_mean_over_std mean(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx))/std(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx))];
        mouse_med=[mouse_med median(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx))];
        mouse_iqr=[mouse_iqr iqr(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx))];
        mouse_iqr_over_med=[mouse_iqr_over_med iqr(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx))/ median(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx))];
        %and now use the intersection of stim/nostim index to find lps that
        %were/not stimulated
        mouse_by_day_stim_idx = intersect(mouse_by_day_idx,stim_idx);
        mouse_by_day_no_stim_idx = intersect(mouse_by_day_idx,no_stim_idx);
        %get means, ipi, etc. for press with stimulation/w/o
        mouse_stim_mean = [mouse_stim_mean mean(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_stim_idx))];
        mouse_no_stim_mean = [mouse_no_stim_mean mean(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_no_stim_idx))];
        %now for n+1
        mouse_by_day_stim_idx_n1 = intersect(mouse_by_day_idx,stim_idx_n1);
        mouse_by_day_no_stim_idx_n1 = intersect(mouse_by_day_idx,no_stim_idx_n1);
        mouse_stim_n1_mean = [mouse_stim_n1_mean mean(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_stim_idx_n1))];
        mouse_no_stim_n1_mean = [mouse_no_stim_n1_mean mean(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_no_stim_idx_n1))];
        mouse_stim_n1_ipi_mean =[mouse_stim_n1_ipi_mean nanmean(T10_Logical_and_Continuous.ipi1(mouse_by_day_stim_idx_n1))];
        mouse_no_stim_n1_ipi_mean =[mouse_no_stim_n1_ipi_mean nanmean(T10_Logical_and_Continuous.ipi1(mouse_by_day_no_stim_idx_n1))];
        mouse_stim_n1_std = [mouse_stim_n1_std nanstd(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_stim_idx_n1))];
        mouse_no_stim_n1_std = [mouse_no_stim_n1_std nanstd(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_no_stim_idx_n1))];
        %we want to match the number of no stim to stim trials to calculate
        %correlations. E.g., how related are press n and n - 1 normally - does this
        %change with stimulation?
        %paste n and n+1 together into an array - will have to pad the n1
        %with nans appropriately to match length
        mouse_stim_lps =T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_stim_idx);
        mouse_stim_n1_lps = T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_stim_idx_n1);
        %use this to match length of n and n+1
        if length(mouse_stim_lps) > length(mouse_stim_n1_lps)
            mouse_stim_n1_lps = [mouse_stim_n1_lps;NaN];
        end
        stim_n_and_n1 = [mouse_stim_lps mouse_stim_n1_lps];
        %all mice lumped together
        mouse_stim_n_and_n1 = [mouse_stim_n_and_n1; stim_n_and_n1];
        %repeat for no stim.
        mouse_no_stim_lps =T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_no_stim_idx);
        mouse_no_stim_n1_lps = T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_no_stim_idx_n1);
        mouse_no_stim_n1_lps(1) =[];
        mouse_no_stim_n1_lps = [mouse_no_stim_n1_lps;NaN];
        if length(mouse_no_stim_lps) > length(mouse_no_stim_n1_lps)
            mouse_no_stim_n1_lps = [mouse_no_stim_n1_lps;NaN];
        end
        if length(mouse_no_stim_lps) < length(mouse_no_stim_n1_lps)
            mouse_no_stim_lps = [mouse_no_stim_lps;NaN];
        end
        no_stim_n_and_n1 = [mouse_no_stim_lps mouse_no_stim_n1_lps];
        mouse_no_stim_n_and_n1 = [mouse_no_stim_n_and_n1; no_stim_n_and_n1];
    end
end
%reshape to get mouse(row) x day(col), mostly for graphing
mouse_stim_mean = reshape(mouse_stim_mean,[total_mice,day_number]);
mouse_no_stim_mean = reshape(mouse_no_stim_mean,[total_mice,day_number]);
mouse_stim_n1_mean = reshape(mouse_stim_n1_mean,[total_mice,day_number]);
mouse_no_stim_n1_mean = reshape(mouse_no_stim_n1_mean,[total_mice,day_number]);
mouse_stim_n1_ipi_mean = reshape(mouse_stim_n1_ipi_mean,[total_mice,day_number]);
mouse_no_stim_n1_ipi_mean = reshape(mouse_no_stim_n1_ipi_mean,[total_mice,day_number]);
mouse_stim_n1_std = reshape(mouse_stim_n1_std,[total_mice,day_number]);
mouse_no_stim_n1_std = reshape(mouse_no_stim_n1_std,[total_mice,day_number]);
mouse_means = reshape(mouse_means,[total_mice,day_number]);
mouse_std = reshape(mouse_std,[total_mice,day_number]);
mouse_mean_over_std =reshape(mouse_mean_over_std,[total_mice,day_number]);
mouse_med=reshape(mouse_med,[total_mice,day_number]);
mouse_iqr=reshape(mouse_iqr,[total_mice,day_number]);
mouse_iqr_over_med=reshape(mouse_iqr_over_med,[total_mice,day_number]);
%% is there a correlation between individual model fit and task performance?
%export to R and use RMcorr on the table 
table_of_r2 =  array2table(r_squared_criteria_All,'VariableNames',{'R2','Criteria'});
table_of_r2_complex =  array2table(r_squared_complex_criteria_All,'VariableNames',{'R2','Criteria'});
%is there a correlation between n-1 coef value and r2?
table_of_r2.Coef = indiv_lme_10_indiv_coef_All(2,:)';
%or of model predictions
table_of_r2.Predictions = correctCIprop_criteria_All(:,1);
%add in mouse number
table_of_r2.Mouse = mouse_number_for_corrs;
%using complex LMEs
table_of_r2_complex.Pred = correctCIprop_criteria_complex_All(:,1);
%add in mouse number
table_of_r2_complex.Mouse = mouse_number_for_corrs;
table_of_r2_complex.Day = day_number_for_corrs;
%% return the indicator variables to categorical for use in LMEs
T10_Logical_and_Continuous.day_indicator = categorical(T10_Logical_and_Continuous.day_indicator);
T10_Logical_and_Continuous.mouse_indicator =categorical(T10_Logical_and_Continuous.mouse_indicator);
T10_Logical_and_Continuous.stim_indicator_All = categorical(T10_Logical_and_Continuous.stim_indicator_All);
T10_Logical_and_Continuous.stim_indicator_n1_All = categorical(T10_Logical_and_Continuous.stim_indicator_n1_All);
T10_Logical_and_Continuous.stim_indicator_n2_All = categorical(T10_Logical_and_Continuous.stim_indicator_n2_All);
T10_Logical_and_Continuous.stim_indicator_n3_All = categorical(T10_Logical_and_Continuous.stim_indicator_n3_All);
%% Complex LME Model
%Use BIC selected variables (selected from a full model that included all
%terms and their interactions with durations up to n -6)
model_spec_reduced   = 'LP_Durations_All ~ n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All + n_minus_one_All + HE_n_1_Indicator_All + n_minus_one_Durations_All:HE_n_1_Indicator_All + ipi1 + n_minus_one_Durations_All:ipi1 + LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
%create LME, and save the model
lme_reduced = fitlme(T10_Logical_and_Continuous,model_spec_reduced,'DummyVarCoding','reference');
lme_reduced_se = lme_reduced.Coefficients.SE;
lme_reduced_coef =  lme_reduced.Coefficients.Estimate;
lme_reduced_name =  lme_reduced.Coefficients.Name;
lme_reduced_pval =  lme_reduced.Coefficients.pValue;
lme_reduced_t =  lme_reduced.Coefficients.tStat;
lme_reduced_df =  lme_reduced.Coefficients.DF;
lme_reduced_upper =  lme_reduced.Coefficients.Upper;
lme_reduced_lower =  lme_reduced.Coefficients.Lower;
lme_reduced_AIC = lme_reduced.ModelCriterion.AIC;
lme_reduced_anova = anova(lme_reduced);
lme_reduced_fstat = lme_reduced_anova.FStat;
lme_reduced_fpval = lme_reduced_anova.pValue;
[b bnames bstats] =randomEffects(lme_reduced);
%% Simple LME Model
%create n10 back durations model continaing control vars for shuffling
model_spec_10_ma_totrew = 'LP_Durations_All~ n_minus_one_Durations_All + HE_n_1_Indicator_All + n_minus_one_All + ipi1 + LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
lme_10_ma_totrew = fitlme(T10_Logical_and_Continuous,model_spec_10_ma_totrew);
lme_10_ma_totrew_se = lme_10_ma_totrew.Coefficients.SE;
lme_10_ma_totrew_coef =  lme_10_ma_totrew.Coefficients.Estimate;
lme_10_ma_totrew_name =  lme_10_ma_totrew.Coefficients.Name;
lme_10_ma_totrew_pval =  lme_10_ma_totrew.Coefficients.pValue;
lme_10_ma_totrew_t =  lme_10_ma_totrew.Coefficients.tStat;
lme_10_ma_totrew_df =  lme_10_ma_totrew.Coefficients.DF;
lme_10_ma_totrew_upper =  lme_10_ma_totrew.Coefficients.Upper;
lme_10_ma_totrew_lower =  lme_10_ma_totrew.Coefficients.Lower;
lme_10_ma_totrew_AIC = lme_10_ma_totrew.ModelCriterion.AIC;
lme_10_ma_totrew_anova = anova(lme_10_ma_totrew);
lme_10_ma_totrew_fstat =lme_10_ma_totrew_anova.FStat;
lme_10_ma_totrew_fpval = lme_10_ma_totrew_anova.pValue;
%% shuffling
%For the simple LME, we will shuffle each term one at a time (e.g., only
%n-3) and compare the beta coef from the order shuffled datasets to the
%actual coeffs. 
%Also have for the complex model, but only to double check BIC results,
%comment out for speed of running
%initialize variables
%simple n1:n3 and control vars
%the size will be the number of coefficients in model by number of shuffles
shuf_n1_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n1_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_o1_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_o1_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_he1_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_he1_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_IPI1_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_IPI1_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n2_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n2_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n3_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n3_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_ts_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_ts_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_crit_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_crit_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
%variables from full beh model
shuf_n1_complex_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n1_complex_SE_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ma_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n1dur_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n2dur_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n3dur_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi_me_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi_n1_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi_avg_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi2_me_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi2_n2_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_rew_me_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_rew_n1_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_rew_avg_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ts_me_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ts_n1_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ts_avg_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_crit_me_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_crit_n1_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_crit_avg_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_he_me_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_he_n1_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_he_avg_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
%using AIC instead of Coefs.
shuf_full_n1dur_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n2dur_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n3dur_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n4dur_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n5dur_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n6dur_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ma_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi_me_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi_n1_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi_avg_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi2_me_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi2_n2_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_rew_me_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_rew_n1_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_rew_avg_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ts_me_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ts_n1_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ts_avg_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_crit_me_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_crit_n1_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_crit_avg_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_he_me_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_he_n1_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_he_avg_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi2_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi2_n2_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
%souter loop is for the number of shuffles made
for j =  1:size(shuf_n_back_Lengths_distribution_Cell,1)
    %need to erase the i_shuffled array every time, this will be used to create
    %the table for regression
    i_shuffled =[];
    %inner loop is for the number of mice x sessions to get n-backs
    for i = 1:size(shuf_n_back_Lengths_distribution_Cell,2)
        %first shuffle for all mice
        i_shuffled = [i_shuffled; shuf_n_back_Lengths_distribution_Cell{j,i}'];
    end
    %get logical reward n-backs inner loop
    shuf_logical_n1=[];
    shuf_logical_n2=[];
    shuf_logical_n3=[];
    for i = 1:size(shuf_n_back_Lengths_distribution_Cell,2)
        shuf_logical_n1 = [shuf_logical_n1; shuf_Logical_n_back_Lengths_distribution_Cell{j,i}(1,:)'];
        shuf_logical_n2 = [shuf_logical_n2; shuf_Logical_n_back_Lengths_distribution_Cell{j,i}(2,:)'];
        shuf_logical_n3 = [shuf_logical_n3; shuf_Logical_n_back_Lengths_distribution_Cell{j,i}(3,:)'];
    end
    %both nback durations and logicals for reward
    i_shuffled = [i_shuffled shuf_logical_n1];
    i_shuffled = [i_shuffled shuf_logical_n2];
    i_shuffled = [i_shuffled shuf_logical_n3];
    %add in the appropriate shuffled n lp (ie, the shuffled dist that built the
    %n-back array)
    i_shuffled = [i_shuffled shuffled_Durs_Dist_All(:,j)];
    %add in appropriate shuffled moving average and other vars
    i_shuffled = [i_shuffled shuffled_total_reward_indicator_All(:,j)];
    i_shuffled = [i_shuffled shuffled_ipi1_All(:,j)];
    i_shuffled = [i_shuffled shuffled_ipi2_All(:,j)];
    i_shuffled = [i_shuffled shuffled_ipi3_All(:,j)];
    i_shuffled = [i_shuffled shuffled_he1_All(:,j)];
    i_shuffled = [i_shuffled shuffled_he2_All(:,j)];
    i_shuffled = [i_shuffled shuffled_he3_All(:,j)];
    i_shuffled = [i_shuffled shuffled_ts_All(:,j)];
    %add in the shuffled indicator variables. As the order of these within day
    %cannot be shuffled, they are shuffled across mice/days
    shuf_criteria_percent_indicator_All = criteria_percent_indicator_All(randperm(length(criteria_percent_indicator_All)));
    shuf_criteria_indicator_All = criteria_indicator_All(randperm(length(criteria_indicator_All)));
    shuf_mouse_indicator = mouse_indicator(randperm(length(mouse_indicator)));
    shuf_day_indicator = day_indicator(randperm(length(day_indicator)));
    i_shuffled = [i_shuffled shuf_mouse_indicator];
    i_shuffled = [i_shuffled shuf_day_indicator];
    i_shuffled = [i_shuffled shuf_criteria_percent_indicator_All];
    i_shuffled = [i_shuffled shuf_criteria_indicator_All];
    %create a table from the shuffled array
    T_shuffled = array2table(i_shuffled,'VariableNames',{'shuf_n_minus_one_Durations_All',...
        'shuf_n_minus_two_Durations_All', 'shuf_n_minus_three_Durations_All', ...
        'shuf_n_minus_one_All','shuf_n_minus_two_All','shuf_n_minus_three_All',...
        'shuf_LP_Durations_All','shuf_total_reward_indicator_All',...
        'shuf_ipi1','shuf_ipi2','shuf_ipi3','shuf_HE_n_1_Indicator_All',...
        'shuf_HE_n_2_Indicator_All','shuf_HE_n_3_Indicator_All','shuf_LP_Timestamps_All',...
        'shuf_mouse_indicator','shuf_day_indicator','shuf_criteria_percent_indicator','shuf_criteria_indicator_All'});
    %create atable with both the actual and shuffled versions of all variables
    T_Both = [T_shuffled T10_Logical_and_Continuous];
    %add in stim variables if needed
    T_Both.shuf_stim_indicator_All = Shuffled_stim_ind_all_All(:,j);
    T_Both.shuf_stim_indicator_n1_All = Shuffled_stim_ind_all_n1_All(:,j);
    T_Both.shuf_stim_indicator_n2_All = Shuffled_stim_ind_all_n2_All(:,j);
    T_Both.shuf_stim_indicator_n3_All = Shuffled_stim_ind_all_n3_All(:,j);
    %% Simple Shuffles
    % shuffle the simple n10 ma tot rew model one at a time for permutation
    %tests
    %N-1 Duration
    model_spec_shuf_n1_simple= 'LP_Durations_All ~ shuf_n_minus_one_Durations_All + HE_n_1_Indicator_All + n_minus_one_All + ipi1 + LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator) + (1|day_indicator)';
    shuf_n1_simple = fitlme(T_Both,model_spec_shuf_n1_simple);
    shuf_n1_simple_coef_all(:,j) =  shuf_n1_simple.Coefficients.Estimate;
    shuf_n1_simple_SE_all(:,j) = shuf_n1_simple.Coefficients.SE;
	%O-1
    model_spec_shuf_o1_simple = 'LP_Durations_All ~ n_minus_one_Durations_All + HE_n_1_Indicator_All + shuf_n_minus_one_All + ipi1 + LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator) + (1|day_indicator)';
    shuf_o1_simple = fitlme(T_Both,model_spec_shuf_o1_simple);
    shuf_o1_simple_coef_all(:,j) =  shuf_o1_simple.Coefficients.Estimate;
    shuf_o1_simple_SE_all(:,j) = shuf_o1_simple.Coefficients.SE;
    %HE-1
    model_spec_shuf_he1_simple = 'LP_Durations_All ~ n_minus_one_Durations_All + shuf_HE_n_1_Indicator_All + n_minus_one_All + ipi1 + LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator) + (1|day_indicator)';
    shuf_he1_simple = fitlme(T_Both,model_spec_shuf_he1_simple);
    shuf_he1_simple_coef_all(:,j) =  shuf_he1_simple.Coefficients.Estimate;
    shuf_he1_simple_SE_all(:,j) = shuf_he1_simple.Coefficients.SE;
    %IPI-1
    model_spec_shuf_IPI1_simple = 'LP_Durations_All ~ n_minus_one_Durations_All + HE_n_1_Indicator_All + n_minus_one_All + shuf_ipi1 + LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator) + (1|day_indicator)';
    shuf_IPI1_simple = fitlme(T_Both,model_spec_shuf_IPI1_simple);
    shuf_IPI1_simple_coef_all(:,j) =  shuf_IPI1_simple.Coefficients.Estimate;
    shuf_IPI1_simple_SE_all(:,j) = shuf_IPI1_simple.Coefficients.SE;
    %Timestamps
    model_spec_shuf_ts_simple= 'LP_Durations_All~ n_minus_one_Durations_All + HE_n_1_Indicator_All + n_minus_one_All + ipi1 + shuf_LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator) + (1|day_indicator)';
    shuf_ts_simple = fitlme(T_Both,model_spec_shuf_ts_simple);
    shuf_ts_simple_coef_all(:,j) = shuf_ts_simple.Coefficients.Estimate;
    shuf_ts_simple_SE_all(:,j) = shuf_ts_simple.Coefficients.SE;
    %Criteria%
    model_spec_shuf_crit_simple= 'LP_Durations_All~ n_minus_one_Durations_All + HE_n_1_Indicator_All + n_minus_one_All + ipi1 + LP_Timestamps_All + shuf_criteria_percent_indicator + (1|mouse_indicator) + (1|day_indicator)';
    shuf_crit_simple = fitlme(T_Both,model_spec_shuf_crit_simple);
    shuf_crit_simple_coef_all(:,j) = shuf_crit_simple.Coefficients.Estimate;
    shuf_crit_simple_SE_all(:,j) = shuf_crit_simple.Coefficients.SE;
    %% Complex Shuffles
    % shuffle the Complex n10 ma tot rew model one at a time for permutation
    %tests
    %N-1 Duration
    model_spec_shuf_n1_complex   = 'LP_Durations_All ~ shuf_n_minus_one_Durations_All + shuf_n_minus_one_Durations_All:n_minus_one_All + n_minus_one_All + HE_n_1_Indicator_All + shuf_n_minus_one_Durations_All:HE_n_1_Indicator_All + ipi1 + shuf_n_minus_one_Durations_All:ipi1 + LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
    shuf_n1_complex = fitlme(T_Both, model_spec_shuf_n1_complex, 'DummyVarCoding','reference');
    shuf_n1_complex_coef_all(:,j) =  shuf_n1_complex.Coefficients.Estimate;
    shuf_n1_complex_SE_all(:,j) = shuf_n1_complex.Coefficients.SE;
    
    %% Complex Shuffles Stim, Uncomment For Optogenetics Experiments
    % shuffle the Complex n10 ma tot rew model one at a time for permutation
    %tests
%     %N-0 Stim
%     model_spec_shuf_n0_stim_complex   = 'LP_Durations_All ~ shuf_stim_indicator_All:n_minus_one_Durations_All + stim_indicator_n1_All:n_minus_one_Durations_All + shuf_stim_indicator_All + stim_indicator_n1_All + n_minus_one_Durations_All + n_minus_one_All + HE_n_1_Indicator_All + ipi1 + n_minus_one_Durations_All:n_minus_one_All + n_minus_one_Durations_All:HE_n_1_Indicator_All +  n_minus_one_Durations_All:ipi1 +  shuf_stim_indicator_All:ipi1 +  stim_indicator_n1_All:ipi1 + LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)'; 
%     shuf_n0_stim_complex = fitlme(T_Both, model_spec_shuf_n0_stim_complex, 'DummyVarCoding','reference');
%     shuf_n0_stim_complex_coef_all(:,j) =  shuf_n0_stim_complex.Coefficients.Estimate;
%     shuf_n0_stim_complex_SE_all(:,j) = shuf_n0_stim_complex.Coefficients.SE;
%     %N-1 Stim
%     model_spec_shuf_n1_stim_complex   = 'LP_Durations_All ~ stim_indicator_All:n_minus_one_Durations_All + shuf_stim_indicator_n1_All:n_minus_one_Durations_All + stim_indicator_All + shuf_stim_indicator_n1_All + n_minus_one_Durations_All + n_minus_one_All + HE_n_1_Indicator_All + ipi1 + n_minus_one_Durations_All:n_minus_one_All + n_minus_one_Durations_All:HE_n_1_Indicator_All +  n_minus_one_Durations_All:ipi1 +  stim_indicator_All:ipi1 +  shuf_stim_indicator_n1_All:ipi1 + LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)'; 
%     shuf_n1_stim_complex = fitlme(T_Both, model_spec_shuf_n1_stim_complex, 'DummyVarCoding','reference');
%     shuf_n1_stim_complex_coef_all(:,j) =  shuf_n1_stim_complex.Coefficients.Estimate;
%     shuf_n1_stim_complex_SE_all(:,j) = shuf_n1_stim_complex.Coefficients.SE;
%     %N-0 and N-1 Stim
%     model_spec_shuf_both_stim_complex   = 'LP_Durations_All ~ shuf_stim_indicator_All:n_minus_one_Durations_All + shuf_stim_indicator_n1_All:n_minus_one_Durations_All + shuf_stim_indicator_All + shuf_stim_indicator_n1_All + n_minus_one_Durations_All + n_minus_one_All + HE_n_1_Indicator_All + ipi1 + n_minus_one_Durations_All:n_minus_one_All + n_minus_one_Durations_All:HE_n_1_Indicator_All +  n_minus_one_Durations_All:ipi1 +  shuf_stim_indicator_All:ipi1 +  shuf_stim_indicator_n1_All:ipi1 + LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)'; 
%     shuf_both_stim_complex = fitlme(T_Both, model_spec_shuf_both_stim_complex, 'DummyVarCoding','reference');
%     shuf_both_stim_complex_coef_all(:,j) =  shuf_both_stim_complex.Coefficients.Estimate;
%     shuf_both_stim_complex_SE_all(:,j) = shuf_both_stim_complex.Coefficients.SE;
end
%% shuffle v acutal permutation tests for simple lme
%asking, how many times is there a value in the shuffled daataset that is
%at least as extreme as the actually obtained value? Set alpha to .05, so
%if, for instance, we have 1000 shuffled datasets, we say there can be no
%more than 50 times that a shuffled variable is larger than the actual.
n_p_o1_simple = sum(abs(shuf_o1_simple_coef_all(find(contains(shuf_o1_simple.Coefficients.Name,'shuf_n_minus_one_All')),:)) > abs(lme_10_ma_totrew_coef(find(contains(lme_10_ma_totrew.Coefficients.Name,'n_minus_one_All_1')))))/size(shuf_o1_simple_coef_all,2);
n_p_n1_simple = sum(abs(shuf_n1_simple_coef_all(find(contains(shuf_n1_simple.Coefficients.Name,'shuf_n_minus_one_Durations_All')),:)) > abs(lme_10_ma_totrew_coef(find(contains(lme_10_ma_totrew.Coefficients.Name,'n_minus_one_Durations_All')))))/size(shuf_n1_simple_coef_all,2);
n_p_crit_simple = sum(abs(shuf_crit_simple_coef_all(find(contains(shuf_crit_simple.Coefficients.Name,'shuf_criteria_percent_indicator')),:)) > abs(lme_10_ma_totrew_coef(find(contains(lme_10_ma_totrew.Coefficients.Name,'criteria_percent_indicator')))))/size(shuf_crit_simple_coef_all,2);
n_p_he1_simple = sum(abs(shuf_he1_simple_coef_all(find(contains(shuf_he1_simple.Coefficients.Name,'shuf_HE_n_1_Indicator_All')),:)) > abs(lme_10_ma_totrew_coef(find(contains(lme_10_ma_totrew.Coefficients.Name,'HE_n_1_Indicator_All_1')))))/size(shuf_he1_simple_coef_all,2);
n_p_ts_simple = sum(abs(shuf_ts_simple_coef_all(find(contains(shuf_ts_simple.Coefficients.Name,'shuf_LP_Timestamps_All')),:)) > abs(lme_10_ma_totrew_coef(find(contains(lme_10_ma_totrew.Coefficients.Name,'LP_Timestamps_All')))))/size(shuf_ts_simple_coef_all,2);
n_p_IPI1_simple = sum(abs(shuf_IPI1_simple_coef_all(find(contains(shuf_IPI1_simple.Coefficients.Name,'shuf_ipi1')),:)) > abs(lme_10_ma_totrew_coef(find(contains(lme_10_ma_totrew.Coefficients.Name,'ipi1')))))/size(shuf_IPI1_simple_coef_all,2);
%get means and sem for each nback
n1_shuf_mean_simple = mean(shuf_n1_simple_coef_all(find(contains(shuf_n1_simple.Coefficients.Name,'shuf_n_minus_one_Durations_All')),:)) ;
n1_shuf_sem_simple = std(shuf_n1_simple_coef_all(find(contains(shuf_n1_simple.Coefficients.Name,'shuf_n_minus_one_Durations_All')),:))/(sqrt(size(shuf_n1_simple_coef_all(find(contains(shuf_n1_simple.Coefficients.Name,'shuf_n_minus_one_Durations_All')),:),2)));
o1_shuf_mean_simple = mean(shuf_o1_simple_coef_all(find(contains(shuf_o1_simple.Coefficients.Name,'shuf_n_minus_one_All')),:)) ;
o1_shuf_sem_simple = std(shuf_o1_simple_coef_all(find(contains(shuf_o1_simple.Coefficients.Name,'shuf_n_minus_one_All')),:))/(sqrt(size(shuf_n1_simple_coef_all(find(contains(shuf_o1_simple.Coefficients.Name,'shuf_n_minus_one_All')),:),2)));
crit_shuf_mean_simple = mean(shuf_crit_simple_coef_all(find(contains(shuf_crit_simple.Coefficients.Name,'shuf_criteria_percent_indicator')),:)) ;
crit_shuf_sem_simple = std(shuf_crit_simple_coef_all(find(contains(shuf_crit_simple.Coefficients.Name,'shuf_criteria_percent_indicator')),:))/(sqrt(size(shuf_crit_simple_coef_all(find(contains(shuf_crit_simple.Coefficients.Name,'shuf_criteria_percent_indicator')),:),2)));
he1_shuf_mean_simple = mean(shuf_he1_simple_coef_all(find(contains(shuf_he1_simple.Coefficients.Name,'shuf_HE_n_1_Indicator_All')),:)) ;
he1_shuf_sem_simple = std(shuf_he1_simple_coef_all(find(contains(shuf_he1_simple.Coefficients.Name,'shuf_HE_n_1_Indicator_All')),:))/(sqrt(size(shuf_he1_simple_coef_all(find(contains(shuf_he1_simple.Coefficients.Name,'shuf_HE_n_1_Indicator_All')),:),2)));
ts_shuf_mean_simple = mean(shuf_ts_simple_coef_all(find(contains(shuf_ts_simple.Coefficients.Name,'shuf_LP_Timestamps_All')),:)) ;
ts_shuf_sem_simple = std(shuf_ts_simple_coef_all(find(contains(shuf_ts_simple.Coefficients.Name,'shuf_LP_Timestamps_All')),:))/(sqrt(size(shuf_ts_simple_coef_all(find(contains(shuf_ts_simple.Coefficients.Name,'shuf_LP_Timestamps_All')),:),2)));
IPI1_shuf_mean_simple = mean(shuf_IPI1_simple_coef_all(find(contains(shuf_IPI1_simple.Coefficients.Name,'shuf_ipi1')),:)) ;
IPI1_shuf_sem_simple = std(shuf_IPI1_simple_coef_all(find(contains(shuf_IPI1_simple.Coefficients.Name,'shuf_ipi1')),:))/(sqrt(size(shuf_IPI1_simple_coef_all(find(contains(shuf_IPI1_simple.Coefficients.Name,'shuf_ipi1')),:),2)));
%put the mean and sem together in coef order for graping
n1_to_n10_shuf_mean = [n1_shuf_mean_simple; o1_shuf_mean_simple; crit_shuf_mean_simple; he1_shuf_mean_simple; ts_shuf_mean_simple; IPI1_shuf_mean_simple;];
n1_to_n10_shuf_sem = [n1_shuf_sem_simple; o1_shuf_sem_simple; crit_shuf_sem_simple; he1_shuf_sem_simple; ts_shuf_sem_simple; IPI1_shuf_sem_simple;];
%put all the permutation pvals together in order of the coefs
n_p_1_to_10_totrew_name_simple = [{'n-1 Duration'}; {'n-1 Outcome'}; {'% Criteria'}; {'n-1 Head Entry'}; {'Lever Press Timestamp'}; {'n-1 Interpress Interval'}];
n_p_1_to_10_totrew_pvals_simple = [n_p_n1_simple; n_p_o1_simple; n_p_crit_simple; n_p_he1_simple; n_p_ts_simple; n_p_IPI1_simple];
n_p_1_to_10_totrew_mean_simple = n1_to_n10_shuf_mean;
n_p_1_to_10_totrew_sem_simple = n1_to_n10_shuf_sem;
%% shuffle v acutal permutation tests for complex lme
%asking, how many times is there a value in the shuffled daataset that is
%at least as extreme as the actually obtained value? Set alpha to .05, so
%if, for instance, we have 1000 shuffled datasets, we say there can be no
%more than 50 times that a shuffled variable is larger than the actual.
n_p_n1_complex = sum(abs(shuf_n1_complex_coef_all(2,:)) > abs(lme_reduced_coef(3)))/size(shuf_n1_complex_coef_all,2);
n_p_n1_outcome_complex = sum(abs(shuf_n1_complex_coef_all(find(contains(shuf_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:n_minus_one_All_1')),:)) > abs(lme_reduced_coef(find(contains(lme_reduced.Coefficients.Name,'n_minus_one_All_1:n_minus_one_Durations_All')))))/size(shuf_n1_complex_coef_all,2);
n_p_n1_he_complex =  sum(abs(shuf_n1_complex_coef_all(find(contains(shuf_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:HE_n_1_Indicator_All_1')),:)) > abs(lme_reduced_coef(find(contains(lme_reduced.Coefficients.Name,'n_minus_one_Durations_All:HE_n_1_Indicator_All_1')))))/size(shuf_n1_complex_coef_all,2);
n_p_n1_ipi_complex = sum(abs(shuf_n1_complex_coef_all(find(contains(shuf_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:ipi1')),:)) > abs(lme_reduced_coef(find(contains(lme_reduced.Coefficients.Name,'n_minus_one_Durations_All:ipi1')))))/size(shuf_n1_complex_coef_all,2);
%get means and sem for each nback
n1_shuf_mean_complex = mean(shuf_n1_complex_coef_all(2,:)) ;
n1_shuf_sem_complex  = std(shuf_n1_complex_coef_all(2,:))/(sqrt(size( shuf_n1_complex_coef_all(2,:),2)));
n1_outcome_shuf_mean_complex = mean(shuf_n1_complex_coef_all(find(contains(shuf_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:n_minus_one_All_1')),:)) ;
n1_outcome_shuf_sem_complex  = std(shuf_n1_complex_coef_all(find(contains(shuf_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:n_minus_one_All_1')),:))/(sqrt(size(shuf_n1_complex_coef_all(find(contains(shuf_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:n_minus_one_All_1')),:),2)));
n1_he_shuf_mean_complex = mean(shuf_n1_complex_coef_all(find(contains(shuf_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:HE_n_1_Indicator_All_1')),:)) ;
n1_he_shuf_sem_complex  = std(shuf_n1_complex_coef_all(find(contains(shuf_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:HE_n_1_Indicator_All_1')),:))/(sqrt(size(shuf_n1_complex_coef_all(find(contains(shuf_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:HE_n_1_Indicator_All_1')),:),2)));
n1_ipi_shuf_mean_complex = mean(shuf_n1_complex_coef_all(find(contains(shuf_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:ipi1')),:)) ;
n1_ipi_shuf_sem_complex  = std(shuf_n1_complex_coef_all(find(contains(shuf_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:ipi1')),:))/(sqrt(size(shuf_n1_complex_coef_all(find(contains(shuf_n1_complex.Coefficients.Name,'shuf_n_minus_one_Durations_All:ipi1')),:),2)));
%put the mean and sem together in coef order for graping
n1_to_n10_shuf_mean_complex = [n1_shuf_mean_complex; n1_outcome_shuf_mean_complex; n1_he_shuf_mean_complex; n1_ipi_shuf_mean_complex];
n1_to_n10_shuf_sem_complex = [n1_shuf_sem_complex; n1_outcome_shuf_sem_complex; n1_he_shuf_sem_complex; n1_ipi_shuf_sem_complex];
%rather than the SEM of the shuffle coefs, lets calculaute the average SEM of the shuffled models 
n1_to_n10_shuf_avg_model_se = [mean(shuf_o1_simple_SE_all(2,:)); mean(shuf_n1_simple_SE_all(2,:)); mean(shuf_crit_simple_SE_all(2,:)); mean(shuf_he1_simple_SE_all(2,:)); mean(shuf_ts_simple_SE_all(2,:)); mean(shuf_IPI1_simple_SE_all(2,:))];
n1_to_n10_shuf_avg_model_se_complex = [mean(shuf_n1_complex_SE_all(2,:))];
%put all the permutation pvals together in order of the coefs
n_p_1_to_10_totrew_name_complex = [{'n-1 Duration'}; {'n-1 Duration * n-1 Outcome'}; {'n-1 Duration * n-1 Head Entry'}; {'n-1 Duration * n-1 Interpress Interval'} ];
n_p_1_to_10_totrew_pvals_complex = [n_p_n1_complex; n_p_n1_outcome_complex; n_p_n1_he_complex; n_p_n1_ipi_complex];
n_p_1_to_10_totrew_mean_complex = [n1_shuf_mean_complex; n1_outcome_shuf_mean_complex; n1_he_shuf_mean_complex; n1_ipi_shuf_mean_complex];
n_p_1_to_10_totrew_sem_complex = [n1_shuf_sem_complex; n1_outcome_shuf_sem_complex; n1_he_shuf_sem_complex; n1_ipi_shuf_sem_complex];
%% Predictions 
%predict LP_Durations from the model (predicted response using model data by default)
 [yhat yhatCI yhatDF]= predict(lme_reduced,'Simultaneous',true);
%calculate how often the actual value lies between the CI for the
 %predicted
 correctish_pred =T10_Logical_and_Continuous.LP_Durations_All <=yhatCI(:,2) &  T10_Logical_and_Continuous.LP_Durations_All >=yhatCI(:,1);
 correctCI_prop = sum(correctish_pred)/length(T10_Logical_and_Continuous.LP_Durations_All);
%predict a smoothed dataset
%smooth the 10 presses to the left
 smootheddata = smoothdata(T10_Logical_and_Continuous.LP_Durations_All,'gaussian',[0 10]);
 correctish_pred_smooth =smootheddata<=yhatCI(:,2) & smootheddata >=yhatCI(:,1);
 correctCI_prop_Smooth = sum(correctish_pred_smooth)/length(smootheddata) ;
%What about predicting order shuffled Durations?
 correctish_shuf_pred =Shuffled_Durations_All <=yhatCI(:,2) & Shuffled_Durations_All >=yhatCI(:,1);
 correctCI_shuf_prop = sum(correctish_shuf_pred)/length(Shuffled_Durations_All);
 %How well does the model simply predict an increase a decrease relative to
%the preceding press?
diff_yhat = diff(yhat);
diff_x = diff(T10_Logical_and_Continuous.LP_Durations_All);
dixx_x_yhat = [diff_x diff_yhat];
%convert the numbers to either -1 or +1 to see how often it predicts a
%increase or decrease accurately
 dixx_x_yhat(dixx_x_yhat > 0) =1;
 dixx_x_yhat(dixx_x_yhat < 0) =-1;
inc_pred =0;
dec_pred = 0;
for i = 1:length(dixx_x_yhat)
    if dixx_x_yhat(i,1) == 1 && dixx_x_yhat(i,2) == 1
        inc_pred = inc_pred+1;
    elseif  dixx_x_yhat(i,1) == -1 && dixx_x_yhat(i,2) == -1
        dec_pred = dec_pred +1;
    end
    
end
actual_inc = sum(dixx_x_yhat(:,1)==1);
correct_inc_pred = inc_pred/actual_inc;
actual_dec=sum(dixx_x_yhat(:,1)==-1);
correct_dec_pred =dec_pred/actual_dec;
sum(dixx_x_yhat(:,1) == dixx_x_yhat(:,2))/length(dixx_x_yhat);
%How well does the model accurately predict rewarded/unrewarded presses?
%use the upper CI to see if the prediciton would be rewarded or not
%need numerical crteria% indicator for comparisons of duration to criterion
T10_Logical_and_Continuous.criteria_indicator_All_numerical = criteria_indicator_All;
pred_rew = yhatCI(:,2) >= T10_Logical_and_Continuous.criteria_indicator_All_numerical; 
pred_rew = double(pred_rew);
% ll here is the n-0 reward indicator var
rew_red_actual = [logical_LP_n_back_DataRaster(:,4) pred_rew];
success_pred =0;
fail_pred = 0;
for i = 1:length(rew_red_actual)
%     for i = (41267:414460 )
%if both the actual and predicted reward are 1
    if rew_red_actual(i,1) == 1 && rew_red_actual(i,2) == 1
        %then iterate the accurate success_pred variable
        success_pred = success_pred+1;
        %ditto for if it agrees on a failed press
    elseif  rew_red_actual(i,1) == 0 && rew_red_actual(i,2) == 0
        fail_pred = fail_pred +1;
    end
end
actual_rewards = sum(logical_LP_n_back_DataRaster(:,4));
correct_reward_pred = success_pred/actual_rewards;
fail_idxs=sum(logical_LP_n_back_DataRaster(:,4) == 0);
correct_fail_pred =fail_pred/fail_idxs;
sum(pred_rew == logical_LP_n_back_DataRaster(:,4))/length(pred_rew);
%% predict with the simple n10 and control variable model
 [yhat_simple yhatCI_simple yhatDF_simple] = predict(lme_10_ma_totrew,'Simultaneous',true);
 %calculate how often the actual value lies between the CI for the
 %predicted
 correctish_pred_simple =T10_Logical_and_Continuous.LP_Durations_All <=yhatCI_simple(:,2) &  T10_Logical_and_Continuous.LP_Durations_All >=yhatCI_simple(:,1);
 correctCI_prop_simple = sum(correctish_pred_simple)/length(T10_Logical_and_Continuous.LP_Durations_All);
 %% for probability add in main effect and interaction of success + rew vs success no reward
% remove n_minus_one_all effects, since it only cares whether duration was long enough, not whether there was reward
% (this only matters if reward is probbalistic)
% instead use the reward_indicator_n1_All variable
% where 0 = fail, no reward, 1 = success, no reward, and 2 = success, yes
% reward
T10_Logical_and_Continuous.reward_indicator_n1_All = categorical(T10_Logical_and_Continuous.reward_indicator_n1_All);
T10_Logical_and_Continuous.reward_indicator_All = categorical(T10_Logical_and_Continuous.reward_indicator_All);
%BIC selected model (Simple)
model_spec_moderate_remove_proability = 'LP_Durations_All ~ n_minus_one_Durations_All + reward_indicator_n1_All + HE_n_1_Indicator_All + ipi1 + LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
lme_moderate_remove_probability = fitlme(T10_Logical_and_Continuous,model_spec_moderate_remove_proability);
lme_moderate_remove_probability_se = lme_moderate_remove_probability.Coefficients.SE;
lme_moderate_remove_probability_coef =  lme_moderate_remove_probability.Coefficients.Estimate;
lme_moderate_remove_probability_name =  lme_moderate_remove_probability.Coefficients.Name;
lme_moderate_remove_probability_pval =  lme_moderate_remove_probability.Coefficients.pValue;
lme_moderate_remove_probability_t =  lme_moderate_remove_probability.Coefficients.tStat;
lme_moderate_remove_probability_df =  lme_moderate_remove_probability.Coefficients.DF;
lme_moderate_remove_probability_upper =  lme_moderate_remove_probability.Coefficients.Upper;
lme_moderate_remove_probability_lower =  lme_moderate_remove_probability.Coefficients.Lower;
lme_moderate_remove_probability_AIC = lme_moderate_remove_probability.ModelCriterion.AIC;
lme_moderate_remove_probability_anova = anova(lme_moderate_remove_probability);
lme_moderate_remove_probability_fstat =lme_moderate_remove_probability_anova.FStat;
lme_moderate_remove_probability_fpval = lme_moderate_remove_probability_anova.pValue;


% %% Stimulation during the lever press (Simple), Uncomment Below for
% Optogenetics Experiments
% 
% % %% stimulation after lmes
% model_spec_reduced_stim = 'LP_Durations_All ~ stim_indicator_All + stim_indicator_n1_All + n_minus_one_Durations_All + n_minus_one_All + HE_n_1_Indicator_All + ipi1 + LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% lme_reduced_stim = fitlme(T10_Logical_and_Continuous,model_spec_reduced_stim, 'DummyVarCoding','reference');
% lme_reduced_stim_se = lme_reduced_stim.Coefficients.SE;
% lme_reduced_stim_coef =  lme_reduced_stim.Coefficients.Estimate;
% lme_reduced_stim_name =  lme_reduced_stim.Coefficients.Name;
% lme_reduced_stim_pval =  lme_reduced_stim.Coefficients.pValue;
% lme_reduced_stim_t =  lme_reduced_stim.Coefficients.tStat;
% lme_reduced_stim_df =  lme_reduced_stim.Coefficients.DF;
% lme_reduced_stim_upper =  lme_reduced_stim.Coefficients.Upper;
% lme_reduced_stim_lower =  lme_reduced_stim.Coefficients.Lower;
% lme_reduced_stim_AIC = lme_reduced_stim.ModelCriterion.AIC;
% lme_reduced_stim_anova = anova(lme_reduced_stim);
% lme_reduced_stim_fstat =lme_reduced_stim_anova.FStat;
% lme_reduced_stim_fpval = lme_reduced_stim_anova.pValue;
% 
% %% Stimulation during the lever press (Complex)
% %model_spec_complex_stim = 'LP_Durations_All ~ stim_indicator_All*n_minus_one_Durations_All + stim_indicator_n1_All*n_minus_one_Durations_All + stim_indicator_All + stim_indicator_n1_All + n_minus_one_Durations_All + n_minus_one_All + HE_n_1_Indicator_All + ipi1 + n_minus_one_Durations_All:n_minus_one_All + n_minus_one_Durations_All:HE_n_1_Indicator_All +  n_minus_one_Durations_All:ipi1 + LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% model_spec_complex_stim = 'LP_Durations_All ~ stim_indicator_All:n_minus_one_Durations_All + stim_indicator_n1_All:n_minus_one_Durations_All + stim_indicator_All + stim_indicator_n1_All + n_minus_one_Durations_All + n_minus_one_All + HE_n_1_Indicator_All + ipi1 + n_minus_one_Durations_All:n_minus_one_All + n_minus_one_Durations_All:HE_n_1_Indicator_All +  n_minus_one_Durations_All:ipi1 +  stim_indicator_All:ipi1 +  stim_indicator_n1_All:ipi1 + LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% lme3_time_stim_n_int_ctl_ma_int_full = fitlme(T10_Logical_and_Continuous, model_spec_complex_stim, 'DummyVarCoding','reference');
% lme3_time_stim_n_int_ctl_ma_int_full_coef = lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.Estimate;
% lme3_time_stim_n_int_ctl_ma_int_full_se = lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.SE;
% lme3_time_stim_n_int_ctl_ma_int_full_pval = lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.pValue;
% lme3_time_stim_n_int_ctl_ma_int_full_name = lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.Name;
% lme3_time_stim_n_int_ctl_ma_int_full_df = lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.DF;
% lme3_time_stim_n_int_ctl_ma_int_full_t = lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.tStat;
% lme3_time_stim_n_int_ctl_ma_int_full_upper = lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.Upper;
% lme3_time_stim_n_int_ctl_ma_int_full_lower = lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.Lower;
% lme3_time_stim_n_int_ctl_ma_int_full_anova = anova(lme3_time_stim_n_int_ctl_ma_int_full);
% lme3_time_stim_n_int_ctl_ma_int_full_fstat = lme3_time_stim_n_int_ctl_ma_int_full_anova.FStat;
% lme3_time_stim_n_int_ctl_ma_int_full_fpval = lme3_time_stim_n_int_ctl_ma_int_full_anova.pValue;
% [b bnames bstats] = randomEffects(lme3_time_stim_n_int_ctl_ma_int_full);
% lme3_time_stim_n_int_ctl_ma_int_full_b = b;
% lme3_time_stim_n_int_ctl_ma_int_full_bnames = bnames;
% lme3_time_stim_n_int_ctl_ma_int_full_bstats = bstats;
% 
% %% shuffle v acutal permutation tests for complex lme
% %asking, how many times is there a value in the shuffled daataset that is
% %at least as extreme as the actually obtained value? Set alpha to .05, so
% %if, for instance, we have 1000 shuffled datasets, we say there can be no
% %more than 50 times that a shuffled variable is larger than the actual.
% 
% n0_p_n1_stim_complex = sum(abs(shuf_n0_stim_complex_coef_all(3,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_full_coef(3)))/size(shuf_n0_stim_complex_coef_all,2);
% n0_p_duration_n1_stim_complex = sum(abs(shuf_n0_stim_complex_coef_all(find(contains(shuf_n0_stim_complex.Coefficients.Name,'n_minus_one_Durations_All:shuf_stim_indicator_All')),:)) > abs(lme3_time_stim_n_int_ctl_ma_int_full_coef(find(contains(lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.Name,'n_minus_one_Durations_All:stim_indicator_All_1')))))/size(shuf_n0_stim_complex_coef_all,2);
% n0_p_ipi_n1_stim_complex =  sum(abs(shuf_n0_stim_complex_coef_all(find(contains(shuf_n0_stim_complex.Coefficients.Name,'ipi1:shuf_stim_indicator_All')),:)) > abs(lme3_time_stim_n_int_ctl_ma_int_full_coef(find(contains(lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.Name,'ipi1:stim_indicator_All_1')))))/size(shuf_n0_stim_complex_coef_all,2);
% 
% n1_p_n1_stim_complex = sum(abs(shuf_n1_stim_complex_coef_all(3,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_full_coef(3)))/size(shuf_n1_stim_complex_coef_all,2);
% n1_p_duration_n1_stim_complex = sum(abs(shuf_n1_stim_complex_coef_all(find(contains(shuf_n1_stim_complex.Coefficients.Name,'n_minus_one_Durations_All:shuf_stim_indicator_n1_All')),:)) > abs(lme3_time_stim_n_int_ctl_ma_int_full_coef(find(contains(lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.Name,'n_minus_one_Durations_All:stim_indicator_n1_All_1')))))/size(shuf_n1_stim_complex_coef_all,2);
% n1_p_ipi_n1_stim_complex =  sum(abs(shuf_n1_stim_complex_coef_all(find(contains(shuf_n1_stim_complex.Coefficients.Name,'ipi1:shuf_stim_indicator_n1_All')),:)) > abs(lme3_time_stim_n_int_ctl_ma_int_full_coef(find(contains(lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.Name,'ipi1:stim_indicator_n1_All_1')))))/size(shuf_n1_stim_complex_coef_all,2);
% 
% n1_p_n1_both_stim_complex = sum(abs(shuf_both_stim_complex_coef_all(3,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_full_coef(3)))/size(shuf_both_stim_complex_coef_all,2);
% n0_p_duration_n1_both_stim_complex = sum(abs(shuf_both_stim_complex_coef_all(find(contains(shuf_both_stim_complex.Coefficients.Name,'n_minus_one_Durations_All:shuf_stim_indicator_All')),:)) > abs(lme3_time_stim_n_int_ctl_ma_int_full_coef(find(contains(lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.Name,'n_minus_one_Durations_All:stim_indicator_n1_All_1')))))/size(shuf_both_stim_complex_coef_all,2);
% n0_p_ipi_n1_both_stim_complex =  sum(abs(shuf_both_stim_complex_coef_all(find(contains(shuf_both_stim_complex.Coefficients.Name,'ipi1:shuf_stim_indicator_All')),:)) > abs(lme3_time_stim_n_int_ctl_ma_int_full_coef(find(contains(lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.Name,'ipi1:stim_indicator_n1_All_1')))))/size(shuf_both_stim_complex_coef_all,2);
% n1_p_duration_n1_both_stim_complex = sum(abs(shuf_both_stim_complex_coef_all(find(contains(shuf_both_stim_complex.Coefficients.Name,'n_minus_one_Durations_All:shuf_stim_indicator_n1_All')),:)) > abs(lme3_time_stim_n_int_ctl_ma_int_full_coef(find(contains(lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.Name,'n_minus_one_Durations_All:stim_indicator_n1_All_1')))))/size(shuf_both_stim_complex_coef_all,2);
% n1_p_ipi_n1_both_stim_complex =  sum(abs(shuf_both_stim_complex_coef_all(find(contains(shuf_both_stim_complex.Coefficients.Name,'ipi1:shuf_stim_indicator_n1_All')),:)) > abs(lme3_time_stim_n_int_ctl_ma_int_full_coef(find(contains(lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.Name,'ipi1:stim_indicator_n1_All_1')))))/size(shuf_both_stim_complex_coef_all,2);
% 
% 
% %get means and sem for each nback
% n0_shuf_mean_n1_stim_complex = mean(shuf_n0_stim_complex_coef_all(3,:));
% n0_shuf_sem_n1_stim_complex = std(shuf_n0_stim_complex_coef_all(3,:))/(sqrt(size(shuf_n0_stim_complex_coef_all(3,:),2)));
% n0_shuf_mean_duration_n1_stim_complex = mean(shuf_n0_stim_complex_coef_all(find(contains(shuf_n0_stim_complex.Coefficients.Name,'n_minus_one_Durations_All:shuf_stim_indicator_All')),:));
% n0_shuf_sem_duration_n1_stim_complex = std(shuf_n0_stim_complex_coef_all(find(contains(shuf_n0_stim_complex.Coefficients.Name,'n_minus_one_Durations_All:shuf_stim_indicator_All')),:))/(sqrt(size(shuf_n0_stim_complex_coef_all(find(contains(shuf_n0_stim_complex.Coefficients.Name,'n_minus_one_Durations_All:shuf_stim_indicator_All')),:),2)));
% n0_shuf_mean_ipi_n1_stim_complex = mean(shuf_n0_stim_complex_coef_all(find(contains(shuf_n0_stim_complex.Coefficients.Name,'ipi1:shuf_stim_indicator_All')),:));
% n0_shuf_sem_ipi_n1_stim_complex = std(shuf_n0_stim_complex_coef_all(find(contains(shuf_n0_stim_complex.Coefficients.Name,'ipi1:shuf_stim_indicator_All')),:))/(sqrt(size(shuf_n0_stim_complex_coef_all(find(contains(shuf_n0_stim_complex.Coefficients.Name,'ipi1:shuf_stim_indicator_All')),:),2)));
% 
% n1_shuf_mean_n1_stim_complex = mean(shuf_n1_stim_complex_coef_all(3,:));
% n1_shuf_sem_n1_stim_complex = std(shuf_n1_stim_complex_coef_all(3,:))/(sqrt(size(shuf_n1_stim_complex_coef_all(3,:),2)));
% n1_shuf_mean_duration_n1_stim_complex = mean(shuf_n1_stim_complex_coef_all(find(contains(shuf_n1_stim_complex.Coefficients.Name,'n_minus_one_Durations_All:shuf_stim_indicator_n1_All')),:));
% n1_shuf_sem_duration_n1_stim_complex = std(shuf_n1_stim_complex_coef_all(find(contains(shuf_n1_stim_complex.Coefficients.Name,'n_minus_one_Durations_All:shuf_stim_indicator_n1_All')),:))/(sqrt(size(shuf_n1_stim_complex_coef_all(find(contains(shuf_n1_stim_complex.Coefficients.Name,'n_minus_one_Durations_All:shuf_stim_indicator_n1_All')),:),2)));
% n1_shuf_mean_ipi_n1_stim_complex = mean(shuf_n1_stim_complex_coef_all(find(contains(shuf_n1_stim_complex.Coefficients.Name,'ipi1:shuf_stim_indicator_n1_All')),:));
% n1_shuf_sem_ipi_n1_stim_complex = std(shuf_n1_stim_complex_coef_all(find(contains(shuf_n1_stim_complex.Coefficients.Name,'ipi1:shuf_stim_indicator_n1_All')),:))/(sqrt(size(shuf_n1_stim_complex_coef_all(find(contains(shuf_n1_stim_complex.Coefficients.Name,'ipi1:shuf_stim_indicator_n1_All')),:),2)));
% 
% both_shuf_mean_n1_stim_complex = mean(shuf_both_stim_complex_coef_all(3,:));
% both_shuf_sem_n1_stim_complex = std(shuf_both_stim_complex_coef_all(3,:))/(sqrt(size(shuf_both_stim_complex_coef_all(3,:),2)));
% both_shuf_mean_duration_n1_stim_complex = mean(shuf_both_stim_complex_coef_all(find(contains(shuf_both_stim_complex.Coefficients.Name,'n_minus_one_Durations_All:shuf_stim_indicator_n1_All')),:));
% both_shuf_sem_duration_n1_stim_complex = std(shuf_both_stim_complex_coef_all(find(contains(shuf_both_stim_complex.Coefficients.Name,'n_minus_one_Durations_All:shuf_stim_indicator_n1_All')),:))/(sqrt(size(shuf_both_stim_complex_coef_all(find(contains(shuf_both_stim_complex.Coefficients.Name,'n_minus_one_Durations_All:shuf_stim_indicator_n1_All')),:),2)));
% both_shuf_mean_ipi_n1_stim_complex = mean(shuf_both_stim_complex_coef_all(find(contains(shuf_both_stim_complex.Coefficients.Name,'ipi1:shuf_stim_indicator_n1_All')),:));
% both_shuf_sem_ipi_n1_stim_complex = std(shuf_both_stim_complex_coef_all(find(contains(shuf_both_stim_complex.Coefficients.Name,'ipi1:shuf_stim_indicator_n1_All')),:))/(sqrt(size(shuf_both_stim_complex_coef_all(find(contains(shuf_both_stim_complex.Coefficients.Name,'ipi1:shuf_stim_indicator_n1_All')),:),2)));
% both_shuf_mean_duration_n0_stim_complex = mean(shuf_both_stim_complex_coef_all(find(contains(shuf_both_stim_complex.Coefficients.Name,'n_minus_one_Durations_All:shuf_stim_indicator_All')),:));
% both_shuf_sem_duration_n0_stim_complex = std(shuf_both_stim_complex_coef_all(find(contains(shuf_both_stim_complex.Coefficients.Name,'n_minus_one_Durations_All:shuf_stim_indicator_All')),:))/(sqrt(size(shuf_both_stim_complex_coef_all(find(contains(shuf_both_stim_complex.Coefficients.Name,'n_minus_one_Durations_All:shuf_stim_indicator_All')),:),2)));
% both_shuf_mean_ipi_n0_stim_complex = mean(shuf_both_stim_complex_coef_all(find(contains(shuf_both_stim_complex.Coefficients.Name,'ipi1:shuf_stim_indicator_All')),:));
% both_shuf_sem_ipi_n0_stim_complex = std(shuf_both_stim_complex_coef_all(find(contains(shuf_both_stim_complex.Coefficients.Name,'ipi1:shuf_stim_indicator_All')),:))/(sqrt(size(shuf_both_stim_complex_coef_all(find(contains(shuf_both_stim_complex.Coefficients.Name,'ipi1:shuf_stim_indicator_All')),:),2)));
% 
% 
% %put the mean and sem together in coef order for graping
% n0_shuf_mean_stim_complex = [n0_shuf_mean_n1_stim_complex; n0_shuf_mean_duration_n1_stim_complex; n0_shuf_mean_ipi_n1_stim_complex];
% n0_shuf_sem_stim_complex = [n0_shuf_sem_n1_stim_complex; n0_shuf_sem_duration_n1_stim_complex; n0_shuf_sem_ipi_n1_stim_complex];
% 
% n1_shuf_mean_stim_complex = [n1_shuf_mean_n1_stim_complex; n1_shuf_mean_duration_n1_stim_complex; n1_shuf_mean_ipi_n1_stim_complex];
% n1_shuf_sem_stim_complex = [n1_shuf_sem_n1_stim_complex; n1_shuf_sem_duration_n1_stim_complex; n1_shuf_sem_ipi_n1_stim_complex];
% 
% both_shuf_mean_stim_complex = [both_shuf_mean_n1_stim_complex; both_shuf_mean_duration_n1_stim_complex; both_shuf_mean_ipi_n1_stim_complex; both_shuf_mean_duration_n0_stim_complex; both_shuf_mean_ipi_n0_stim_complex];
% both_shuf_sem_stim_complex = [both_shuf_sem_n1_stim_complex; both_shuf_sem_duration_n1_stim_complex; both_shuf_sem_ipi_n1_stim_complex; both_shuf_sem_duration_n0_stim_complex; both_shuf_sem_ipi_n0_stim_complex];
% 
% %put all the permutation pvals together in order of the coefs
% n_0_name_stim_complex = [{'n-1 Duration'}; {'n-1 Duration * n-0 Stim'}; {'n-1 Interpress Interval * n-0 Stim'}];
% n_0_pvals_stim_complex = [n0_p_n1_stim_complex; n0_p_duration_n1_stim_complex; n0_p_ipi_n1_stim_complex];
% n_0_mean_stim_complex = n0_shuf_mean_stim_complex;
% n_0_sem_stim_complex = n0_shuf_sem_stim_complex;
% 
% n_1_name_stim_complex = [{'n-1 Duration'}; {'n-1 Duration * n-1 Stim'}; {'n-1 Interpress Interval * n-1 Stim'}];
% n_1_pvals_stim_complex = [n1_p_n1_stim_complex; n1_p_duration_n1_stim_complex; n1_p_ipi_n1_stim_complex];
% n_1_mean_stim_complex = n1_shuf_mean_stim_complex;
% n_1_sem_stim_complex = n1_shuf_sem_stim_complex;
% 
% n_both_name_stim_complex = [{'n-1 Duration'}; {'n-1 Duration * n-1 Stim'}; {'n-1 Interpress Interval * n-1 Stim'}; {'n-1 Duration * n-0 Stim'}; {'n-1 Interpress Interval * n-0 Stim'}];
% n_both_pvals_stim_complex = [n1_p_n1_both_stim_complex; n0_p_duration_n1_both_stim_complex; n0_p_ipi_n1_both_stim_complex; n1_p_duration_n1_both_stim_complex; n1_p_ipi_n1_both_stim_complex];
% n_both_mean_stim_complex = both_shuf_mean_stim_complex;
% n_both_sem_stim_complex = both_shuf_sem_stim_complex;

% %% Treatment comparisons

T10_Logical_and_Continuous.ID_indicator = ID_indicator; %each lp tagged with session number
T_Both.ID_indicator = ID_indicator;

%Find the index of each unique mouse ID
%Sham IDs
Control_mice = {}; %Sham
% Lesion IDs
Exp_mice = {}; %Lesion
% 
%% PV IDS
% % YFP IDs
%Control_mice = {};
% % CH2 IDs
%Exp_mice = {};
%% CamKII IDS
% % % YFP IDs
% Control_mice = {}; %YFP
% % CH2 IDs
% % Exp_mice = {}; % CH2
%% Segment Data
% For each Control Mouse
for mouse = 1:length(Control_mice)

    Current_Mouse = Control_mice{mouse};
    Index = find(contains(T10_Logical_and_Continuous.ID_indicator,Current_Mouse));
    T10_Logical_and_Continuous.Treatment(Index) = 0;
   % T_Both.Treatment(Index) = 0;
end
% For each Experimental Mouse
for mouse = 1:length(Exp_mice)

    Current_Mouse = Exp_mice{mouse};
    Index = find(contains(T10_Logical_and_Continuous.ID_indicator,Current_Mouse));
    T10_Logical_and_Continuous.Treatment(Index) = 1;
    %T_Both.Treatment(Index) = 1;
end
%make sure logical variables are coded as categorical
T10_Logical_and_Continuous.ID_indicator = categorical(T10_Logical_and_Continuous.ID_indicator);
T10_Logical_and_Continuous.Treatment = categorical(T10_Logical_and_Continuous.Treatment);
%T_Both.mouse_indicator = categorical(T_Both.Treatment);


%% Simple LME Model
model_spec_10_ma_totrew_Treatment = 'LP_Durations_All~ n_minus_one_Durations_All + HE_n_1_Indicator_All + n_minus_one_All + ipi1 + LP_Timestamps_All + criteria_percent_indicator + Treatment + (1|mouse_indicator)+(1|day_indicator)';
lme_10_ma_totrew_Treatment = fitlme(T10_Logical_and_Continuous,model_spec_10_ma_totrew_Treatment);
lme_10_ma_totrew_Treatment_se = lme_10_ma_totrew_Treatment.Coefficients.SE;
lme_10_ma_totrew_Treatment_coef =  lme_10_ma_totrew_Treatment.Coefficients.Estimate;
lme_10_ma_totrew_Treatment_name =  lme_10_ma_totrew_Treatment.Coefficients.Name;
lme_10_ma_totrew_Treatment_pval =  lme_10_ma_totrew_Treatment.Coefficients.pValue;
lme_10_ma_totrew_Treatment_t =  lme_10_ma_totrew_Treatment.Coefficients.tStat;
lme_10_ma_totrew_Treatment_df =  lme_10_ma_totrew_Treatment.Coefficients.DF;
lme_10_ma_totrew_Treatment_upper =  lme_10_ma_totrew_Treatment.Coefficients.Upper;
lme_10_ma_totrew_Treatment_lower =  lme_10_ma_totrew_Treatment.Coefficients.Lower;
lme_10_ma_totrew_Treatment_AIC = lme_10_ma_totrew_Treatment.ModelCriterion.AIC;
lme_10_ma_totrew_Treatment_anova = anova(lme_10_ma_totrew_Treatment);
lme_10_ma_totrew_Treatment_fstat =lme_10_ma_totrew_Treatment_anova.FStat;
lme_10_ma_totrew_Treatment_fpval = lme_10_ma_totrew_Treatment_anova.pValue;

%% Complex LME Model
model_spec_reduced_Treatment   = 'LP_Durations_All ~ n_minus_one_Durations_All + n_minus_one_All + HE_n_1_Indicator_All + ipi1 + LP_Timestamps_All + criteria_percent_indicator + Treatment + Treatment:n_minus_one_Durations_All + Treatment:n_minus_one_All + Treatment:ipi1 + Treatment:HE_n_1_Indicator_All + (1|mouse_indicator)+(1|day_indicator)';
lme_reduced_Treatment = fitlme(T10_Logical_and_Continuous,model_spec_reduced_Treatment,'DummyVarCoding','reference');
lme_reduced_Treatment_se = lme_reduced_Treatment.Coefficients.SE;
lme_reduced_Treatment_coef =  lme_reduced_Treatment.Coefficients.Estimate;
lme_reduced_Treatment_name =  lme_reduced_Treatment.Coefficients.Name;
lme_reduced_Treatment_pval =  lme_reduced_Treatment.Coefficients.pValue;
lme_reduced_Treatment_t =  lme_reduced_Treatment.Coefficients.tStat;
lme_reduced_Treatment_df =  lme_reduced_Treatment.Coefficients.DF;
lme_reduced_Treatment_upper =  lme_reduced_Treatment.Coefficients.Upper;
lme_reduced_Treatment_lower =  lme_reduced_Treatment.Coefficients.Lower;
lme_reduced_Treatment_AIC = lme_reduced_Treatment.ModelCriterion.AIC;
lme_reduced_Treatment_anova = anova(lme_reduced_Treatment);
lme_reduced_Treatment_fstat = lme_reduced_Treatment_anova.FStat;
lme_reduced_Treatment_fpval = lme_reduced_Treatment_anova.pValue;
[b bnames bstats] =randomEffects(lme_reduced_Treatment);
lme_reduced_Treatment_b = b;
lme_reduced_Treatment_bnames = bnames;
lme_reduced_Treatment_bstats = bstats;

%% Shuffled stim
shuf_n1_simple_coef_all_treatment = nan(size(lme_10_ma_totrew_Treatment_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n1_simple_SE_all_treatment = nan(size(lme_10_ma_totrew_Treatment_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n1_complex_coef_all_treatment = nan(size(lme_reduced_Treatment_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n1_complex_SE_all_treatment = nan(size(lme_reduced_Treatment_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_o1_complex_coef_all_treatment = nan(size(lme_reduced_Treatment_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_o1_complex_SE_all_treatment = nan(size(lme_reduced_Treatment_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_he1_complex_coef_all_treatment = nan(size(lme_reduced_Treatment_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_he1_complex_SE_all_treatment = nan(size(lme_reduced_Treatment_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_IPI1_complex_coef_all_treatment = nan(size(lme_reduced_Treatment_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_IPI1_complex_SE_all_treatment = nan(size(lme_reduced_Treatment_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
%souter loop is for the number of shuffles made
for j =  1:size(shuf_n_back_Lengths_distribution_Cell,1)
    %need to erase the i_shuffled array every time, this will be used to create
    %the table for regression
    i_shuffled =[];
    %inner loop is for the number of mice x sessions to get n-backs
    for i = 1:size(shuf_n_back_Lengths_distribution_Cell,2)
        %first shuffle for all mice
        i_shuffled = [i_shuffled; shuf_n_back_Lengths_distribution_Cell{j,i}'];
    end
    %get logical reward n-backs inner loop
    shuf_logical_n1=[];
    shuf_logical_n2=[];
    shuf_logical_n3=[];
    for i = 1:size(shuf_n_back_Lengths_distribution_Cell,2)
        shuf_logical_n1 = [shuf_logical_n1; shuf_Logical_n_back_Lengths_distribution_Cell{j,i}(1,:)'];
        shuf_logical_n2 = [shuf_logical_n2; shuf_Logical_n_back_Lengths_distribution_Cell{j,i}(2,:)'];
        shuf_logical_n3 = [shuf_logical_n3; shuf_Logical_n_back_Lengths_distribution_Cell{j,i}(3,:)'];
    end
    %both nback durations and logicals for reward
    i_shuffled = [i_shuffled shuf_logical_n1];
    i_shuffled = [i_shuffled shuf_logical_n2];
    i_shuffled = [i_shuffled shuf_logical_n3];
    %add in the appropriate shuffled n lp (ie, the shuffled dist that built the
    %n-back array)
    i_shuffled = [i_shuffled shuffled_Durs_Dist_All(:,j)];
    %add in appropriate shuffled moving average and other vars
    i_shuffled = [i_shuffled shuffled_total_reward_indicator_All(:,j)];
    i_shuffled = [i_shuffled shuffled_ipi1_All(:,j)];
    i_shuffled = [i_shuffled shuffled_ipi2_All(:,j)];
    i_shuffled = [i_shuffled shuffled_ipi3_All(:,j)];
    i_shuffled = [i_shuffled shuffled_he1_All(:,j)];
    i_shuffled = [i_shuffled shuffled_he2_All(:,j)];
    i_shuffled = [i_shuffled shuffled_he3_All(:,j)];
    i_shuffled = [i_shuffled shuffled_ts_All(:,j)];
    
    %add in the shuffled indicator variables. As the order of these within day
    %cannot be shuffled, they are shuffled across mice/days
    shuf_criteria_percent_indicator_All = criteria_percent_indicator_All(randperm(length(criteria_percent_indicator_All)));
    shuf_criteria_indicator_All = criteria_indicator_All(randperm(length(criteria_indicator_All)));
    shuf_mouse_indicator = mouse_indicator(randperm(length(mouse_indicator)));
    shuf_day_indicator = day_indicator(randperm(length(day_indicator)));
    
    i_shuffled = [i_shuffled shuf_mouse_indicator];
    i_shuffled = [i_shuffled shuf_day_indicator];
    i_shuffled = [i_shuffled shuf_criteria_percent_indicator_All];
    i_shuffled = [i_shuffled shuf_criteria_indicator_All];

    %create a table from the shuffled array
    T_shuffled = array2table(i_shuffled,'VariableNames',{'shuf_n_minus_one_Durations_All',...
        'shuf_n_minus_two_Durations_All', 'shuf_n_minus_three_Durations_All', ...
        'shuf_n_minus_one_All','shuf_n_minus_two_All','shuf_n_minus_three_All',...
        'shuf_LP_Durations_All','shuf_total_reward_indicator_All',...
        'shuf_ipi1','shuf_ipi2','shuf_ipi3','shuf_HE_n_1_Indicator_All',...
        'shuf_HE_n_2_Indicator_All','shuf_HE_n_3_Indicator_All','shuf_LP_Timestamps_All',...
        'shuf_mouse_indicator','shuf_day_indicator','shuf_criteria_percent_indicator','shuf_criteria_indicator_All'});
    
    %create atable with both the actual and shuffled versions of all variables
    T_Both = [T_shuffled T10_Logical_and_Continuous];
    
    %% Simple Shuffles
    % shuffle the Complex n10 ma tot rew model one at a time for permutation
    %tests
    %N-1
    model_spec_shuf_n1_simple_treatment   = 'LP_Durations_All~ shuf_n_minus_one_Durations_All + HE_n_1_Indicator_All + n_minus_one_All + ipi1 + LP_Timestamps_All + criteria_percent_indicator + Treatment + (1|mouse_indicator)+(1|day_indicator)';
    shuf_n1_simple_treatment = fitlme(T_Both, model_spec_shuf_n1_simple_treatment);
    shuf_n1_simple_coef_all_treatment(:,j) =  shuf_n1_simple_treatment.Coefficients.Estimate;
    shuf_n1_simple_SE_all_treatment(:,j) = shuf_n1_simple_treatment.Coefficients.SE;
    
    %% Complex Shuffles
    % shuffle the Complex n10 ma tot rew model one at a time for permutation
    %tests
    %N-1
    model_spec_shuf_n1_complex_treatment   = 'LP_Durations_All ~ shuf_n_minus_one_Durations_All + n_minus_one_All + HE_n_1_Indicator_All + ipi1 + LP_Timestamps_All + criteria_percent_indicator + Treatment + Treatment:shuf_n_minus_one_Durations_All + Treatment:n_minus_one_All + Treatment:ipi1 + Treatment:HE_n_1_Indicator_All + (1|mouse_indicator)+(1|day_indicator)';
    shuf_n1_complex_treatment = fitlme(T_Both, model_spec_shuf_n1_complex_treatment, 'DummyVarCoding','reference');
    shuf_n1_complex_coef_all_treatment(:,j) =  shuf_n1_complex_treatment.Coefficients.Estimate;
    shuf_n1_complex_SE_all_treatment(:,j) = shuf_n1_complex_treatment.Coefficients.SE;
    
    %O-1
    model_spec_shuf_o1_complex_treatment   = 'LP_Durations_All ~ n_minus_one_Durations_All + shuf_n_minus_one_All + HE_n_1_Indicator_All + ipi1 + LP_Timestamps_All + criteria_percent_indicator + Treatment + Treatment:n_minus_one_Durations_All + Treatment:shuf_n_minus_one_All + Treatment:ipi1 + Treatment:HE_n_1_Indicator_All + (1|mouse_indicator)+(1|day_indicator)';
    shuf_o1_complex_treatment = fitlme(T_Both, model_spec_shuf_o1_complex_treatment, 'DummyVarCoding','reference');
    shuf_o1_complex_coef_all_treatment(:,j) =  shuf_o1_complex_treatment.Coefficients.Estimate;
    shuf_o1_complex_SE_all_treatment(:,j) = shuf_o1_complex_treatment.Coefficients.SE;
    
    %HE-1
    model_spec_shuf_he1_complex_treatment   = 'LP_Durations_All ~ n_minus_one_Durations_All + n_minus_one_All + shuf_HE_n_1_Indicator_All + ipi1 + LP_Timestamps_All + criteria_percent_indicator + Treatment + Treatment:n_minus_one_Durations_All + Treatment:n_minus_one_All + Treatment:ipi1 + Treatment:shuf_HE_n_1_Indicator_All + (1|mouse_indicator)+(1|day_indicator)';
    shuf_he1_complex_treatment = fitlme(T_Both, model_spec_shuf_he1_complex_treatment, 'DummyVarCoding','reference');
    shuf_he1_complex_coef_all_treatment(:,j) =  shuf_he1_complex_treatment.Coefficients.Estimate;
    shuf_he1_complex_SE_all_treatment(:,j) = shuf_he1_complex_treatment.Coefficients.SE;
    
    %IPI-1
    model_spec_shuf_IPI1_complex_treatment   = 'LP_Durations_All ~ n_minus_one_Durations_All + n_minus_one_All + HE_n_1_Indicator_All + shuf_ipi1 + LP_Timestamps_All + criteria_percent_indicator + Treatment + Treatment:n_minus_one_Durations_All + Treatment:n_minus_one_All + Treatment:shuf_ipi1 + Treatment:HE_n_1_Indicator_All + (1|mouse_indicator)+(1|day_indicator)';
    shuf_IPI1_complex_treatment = fitlme(T_Both, model_spec_shuf_IPI1_complex_treatment, 'DummyVarCoding','reference');
    shuf_IPI1_complex_coef_all_treatment(:,j) =  shuf_IPI1_complex_treatment.Coefficients.Estimate;
    shuf_IPI1_complex_SE_all_treatment(:,j) = shuf_IPI1_complex_treatment.Coefficients.SE;
end
%% shuffle v acutal permutation tests for simple lme
%asking, how many times is there a value in the shuffled daataset that is
%at least as extreme as the actually obtained value? Set alpha to .05, so
%if, for instance, we have 1000 shuffled datasets, we say there can be no
%more than 50 times that a shuffled variable is larger than the actual.
n_p_n1_treatment_simple = sum(abs(shuf_n1_simple_coef_all_treatment(find(contains(shuf_n1_simple_treatment.Coefficients.Name,'shuf_n_minus_one_Durations_All')),:)) > abs(lme_10_ma_totrew_Treatment_coef(find(contains(lme_10_ma_totrew_Treatment.Coefficients.Name,'n_minus_one_Durations_All')))))/size(shuf_n1_simple_coef_all,2);

%get means and sem for each nback
n1_shuf_mean_treatment_simple = mean(shuf_n1_simple_coef_all_treatment(2,:)) ;
n1_shuf_sem_treatment_simple  = std(shuf_n1_simple_coef_all_treatment(2,:))/(sqrt(size(shuf_n1_simple_coef_all_treatment(2,:),2)));

%put the mean and sem together in coef order for graping
n1_to_n10_shuf_mean_treatment_simple = [n1_shuf_mean_treatment_simple];
n1_to_n10_shuf_sem_treatment_simple = [n1_shuf_sem_treatment_simple];

%rather than the SEM of the shuffle coefs, lets calculaute the average SEM of the shuffled models 
n1_to_n10_shuf_avg_model_se_simple_treatment = [mean(shuf_n1_simple_SE_all_treatment(2,:))];

%put all the permutation pvals together in order of the coefs
n_p_1_to_10_totrew_name_treatment_simple = [{'n-1 Duration'}];
n_p_1_to_10_totrew_pval_treatment_simple = [n_p_n1_treatment_simple];
n_p_1_to_10_totrew_mean_treatment_simple = [n1_to_n10_shuf_mean_treatment_simple];
n_p_1_to_10_totrew_sem_treatment_simple = [n1_to_n10_shuf_sem_treatment_simple];

%% shuffle v acutal permutation tests for complex lme
%asking, how many times is there a value in the shuffled daataset that is
%at least as extreme as the actually obtained value? Set alpha to .05, so
%if, for instance, we have 1000 shuffled datasets, we say there can be no
%more than 50 times that a shuffled variable is larger than the actual.
n_p_n1_treatment_complex = sum(abs(shuf_n1_complex_coef_all_treatment(find(contains(shuf_n1_complex_treatment.Coefficients.Name,'shuf_n_minus_one_Durations_All:Treatment')),:)) > abs(lme_reduced_Treatment_coef(find(contains(lme_reduced_Treatment.Coefficients.Name,'n_minus_one_Durations_All:Treatment')))))/size(shuf_n1_complex_coef_all,2);
n_p_n1_outcome_treatment_complex = sum(abs(shuf_o1_complex_coef_all_treatment(find(contains(shuf_o1_complex_treatment.Coefficients.Name,'shuf_n_minus_one_All:Treatment')),:)) > abs(lme_reduced_Treatment_coef(find(contains(lme_reduced_Treatment.Coefficients.Name,'n_minus_one_All_1:Treatment')))))/size(shuf_n1_complex_coef_all,2);
n_p_n1_he_treatment_complex = sum(abs(shuf_he1_complex_coef_all_treatment(find(contains(shuf_he1_complex_treatment.Coefficients.Name,'shuf_HE_n_1_Indicator_All:Treatment')),:)) > abs(lme_reduced_Treatment_coef(find(contains(lme_reduced_Treatment.Coefficients.Name,'HE_n_1_Indicator_All_1:Treatment')))))/size(shuf_n1_complex_coef_all,2);
n_p_n1_ipi_treatment_complex = sum(abs(shuf_IPI1_complex_coef_all_treatment(find(contains(shuf_IPI1_complex_treatment.Coefficients.Name,'shuf_ipi1:Treatment')),:)) > abs(lme_reduced_Treatment_coef(find(contains(lme_reduced_Treatment.Coefficients.Name,'ipi1:Treatment')))))/size(shuf_n1_complex_coef_all,2);

%get means and sem for each nback
n1_shuf_mean_treatment_complex = mean(shuf_n1_complex_coef_all_treatment(2,:)) ;
n1_shuf_sem_treatment_complex  = std(shuf_n1_complex_coef_all_treatment(2,:))/(sqrt(size(shuf_n1_complex_coef_all_treatment(2,:),2)));
n1_outcome_shuf_mean_treatment_complex = mean(shuf_o1_complex_coef_all_treatment(find(contains(shuf_o1_complex_treatment.Coefficients.Name,'shuf_n_minus_one_All:Treatment')),:)) ;
n1_outcome_shuf_sem_treatment_complex  = std(shuf_o1_complex_coef_all_treatment(find(contains(shuf_o1_complex_treatment.Coefficients.Name,'shuf_n_minus_one_All:Treatment')),:))/(sqrt(size(shuf_o1_complex_coef_all_treatment(find(contains(shuf_o1_complex_treatment.Coefficients.Name,'shuf_n_minus_one_All:Treatment')),:),2)));
n1_he1_shuf_mean_treatment_complex = mean(shuf_he1_complex_coef_all_treatment(find(contains(shuf_he1_complex_treatment.Coefficients.Name,'shuf_HE_n_1_Indicator_All:Treatment')),:)) ;
n1_he1_shuf_sem_treatment_complex  = std(shuf_he1_complex_coef_all_treatment(find(contains(shuf_he1_complex_treatment.Coefficients.Name,'shuf_HE_n_1_Indicator_All:Treatment')),:))/(sqrt(size(shuf_he1_complex_coef_all_treatment(find(contains(shuf_he1_complex_treatment.Coefficients.Name,'shuf_HE_n_1_Indicator_All:Treatment')),:),2)));
n1_ipi_shuf_mean_treatment_complex = mean(shuf_IPI1_complex_coef_all_treatment(find(contains(shuf_IPI1_complex_treatment.Coefficients.Name,'shuf_ipi1:Treatment')),:)) ;
n1_ipi_shuf_sem_treatment_complex  = std(shuf_IPI1_complex_coef_all_treatment(find(contains(shuf_IPI1_complex_treatment.Coefficients.Name,'shuf_ipi1:Treatment')),:))/(sqrt(size(shuf_IPI1_complex_coef_all_treatment(find(contains(shuf_IPI1_complex_treatment.Coefficients.Name,'shuf_ipi1:Treatment')),:),2)));

%put the mean and sem together in coef order for graping
n1_to_n10_shuf_mean_treatment_complex = [n1_shuf_mean_treatment_complex; n1_outcome_shuf_mean_treatment_complex; n1_he1_shuf_mean_treatment_complex; n1_ipi_shuf_mean_treatment_complex];
n1_to_n10_shuf_sem_treatment_complex = [n1_shuf_sem_treatment_complex; n1_outcome_shuf_sem_treatment_complex; n1_he1_shuf_sem_treatment_complex; n1_ipi_shuf_sem_treatment_complex];

%rather than the SEM of the shuffle coefs, lets calculaute the average SEM of the shuffled models 
n1_to_n10_shuf_avg_model_se_complex_treatment = [mean(shuf_n1_complex_SE_all_treatment(2,:))];

%put all the permutation pvals together in order of the coefs
n_p_1_to_10_totrew_name_treatment_complex = [{'n-1 Duration * Treatment'}; {'n-1 Outcome * Treatment'}; {'n-1 Head Entry * Treatment'}; {'n-1 Interpress Interval * Treatment'} ];
n_p_1_to_10_totrew_pvals_treatment_complex = [n_p_n1_treatment_complex; n_p_n1_outcome_treatment_complex; n_p_n1_he_treatment_complex; n_p_n1_ipi_treatment_complex];
n_p_1_to_10_totrew_mean_treatment_complex = n1_to_n10_shuf_mean_treatment_complex;
n_p_1_to_10_totrew_sem_treatment_complex = n1_to_n10_shuf_sem_treatment_complex;
%% R2 split by treatment
mouse_list = unique(T10_Logical_and_Continuous.ID_indicator);
mouse_number = nan(length(mouse_list),1);
treatment_ID = nan(length(mouse_list),1);
table_of_r2_complex_treatment = [];
for mouse = 1:length(mouse_list)
    mouse_ID = mouse_list(mouse);
    Index = find(T10_Logical_and_Continuous.ID_indicator == mouse_ID);
    mouse_number(mouse) = unique(T10_Logical_and_Continuous.mouse_indicator(Index));
    treatment_ID(mouse) = unique(T10_Logical_and_Continuous.Treatment(Index));
    
    r2_index = find(table_of_r2_complex.Mouse == mouse_number(mouse));
    Treatment_table = table(ones(length(r2_index),1)*treatment_ID(mouse),'VariableNames',{'Treatment'});
    table_of_r2_complex_treatment = [table_of_r2_complex_treatment; [table_of_r2_complex(r2_index,:), Treatment_table]]; 
end
table_of_r2_complex_treatment = sortrows(table_of_r2_complex_treatment,{'Treatment', 'Mouse'});


toc   
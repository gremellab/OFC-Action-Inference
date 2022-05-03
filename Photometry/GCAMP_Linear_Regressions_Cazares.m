function [GCAMP] = GCAMP_Linear_Regressions_Cazares(GCAMP, base_time_start, base_time_end, time_end)

%want to create regressions to predict GCAMP activity given:
%current press duration
%n back press durations
%n back outcomes
%HEs
%and maybe other stuff
%start with relatively simple model
%GCAMP~current duration and n - 1 duration
Durations = GCAMP.HoldDown_times;
%calculate moving avg. from n-2 to n-116
ma_length = 60;
moving_average_lp_length = movmean(Durations,[ma_length,0]);
moving_average_lp_length = [NaN; NaN; moving_average_lp_length];
moving_average_lp_length = moving_average_lp_length(1:end-2);
%n7 and further back moving avg to exlude presses n1-n6
 moving_average_lp_length_n7andback = [NaN; NaN; NaN; NaN; NaN; moving_average_lp_length];
    moving_average_lp_length_n7andback = moving_average_lp_length_n7andback(1:end-5);
shuff_num = 1000;
%% total cumulative rewards across session
total_reward_indicator = [];
reward_indicator = 0;
for q = 1:length(Durations)
    total_reward_indicator = [total_reward_indicator; reward_indicator];
        if Durations(q) >= GCAMP.Criteria
    reward_indicator = reward_indicator+1;
        end
end
          
        %% n-x lever press lengths
       n_back_Lengths = {};
           n_lengths_array = [];
           Logical_n_lengths_array= [];
           Logical_n_back_Lengths ={};
          press_indices = 1:length(Durations);
           press_indices = press_indices';
                      %for n-backs 1 to 10
       for x = 1:3 
       for n_it = 1:length(Durations)
           %need to check that if there are at least x presses before the
           %current state press.
           if  press_indices(n_it) <= x
           %if there are not, then put a NaN there
           n_back = NaN;
           elseif  press_indices(n_it) > x
           n_back = Durations(n_it- x);               
           end
           n_lengths_array = [n_lengths_array n_back];
           
           %add in logical lengths 0/1 for failure/success
           Logical_n_lengths_array = n_lengths_array >= GCAMP.Criteria;
                       
       end
             
         n_back_Lengths = [n_back_Lengths; n_lengths_array]; 
         n_lengths_array =[];
         Logical_n_back_Lengths = [Logical_n_back_Lengths; Logical_n_lengths_array];  
         Logical_n_lengths_array =[];
       end
       
       %the logical tests turn NaNs into 0, but we need the NaNs to pad the
       %n-back array
       %loop through the array and convert the first x logical values to
       %Nans to make sure they line up with the absolute values
       for i = 1:length(Logical_n_back_Lengths)
           Logical_n_back_Lengths{i,1} =double(Logical_n_back_Lengths{i,1});
           Logical_n_back_Lengths{i,1}(1,1:i) = NaN;
       end
    
%% create data table with the behavioral info.
%can do within a single animal/day linear regression just to check stuff
%out
n_back_All = [];
for i = 1:length(n_back_Lengths)
    n_back_All = [n_back_All; n_back_Lengths{i}];
end
n_back_All = n_back_All';

  data_table = array2table(n_back_All(:,1:3),'VariableNames',{'n_minus_one_Durations_All',...
    'n_minus_two_Durations_All', 'n_minus_three_Durations_All'});

data_table.n_minus_zero_Durations_All = Durations;
%re-arrange order so n-0 is first
data_table = [data_table(:,end) data_table(:,1) data_table(:,2:end-1)];
%add in all other beh variables
%moving average, hen-1, criteria%, n-1 upstate, n-1 reward or not, hen-2,
%ipi1, ipi2
%this lines up ipi with n-1
ipi1 = [NaN; diff(GCAMP.LP_ON_timestamps)];
ipi2 = [NaN; ipi1(1:end-1)];
ipi3 = [NaN; ipi2(1:end-1)];
data_table.ipi1 = ipi1;
data_table.ipi2 = ipi2;
data_table.ipi3 = ipi3;

data_table.n_minus_one_reward_All =  data_table.n_minus_one_Durations_All >=GCAMP.Criteria;
data_table.n_minus_one_reward_All =categorical(data_table.n_minus_one_reward_All);

data_table.n_minus_two_reward_All =  data_table.n_minus_two_Durations_All >=GCAMP.Criteria;
data_table.n_minus_two_reward_All =categorical(data_table.n_minus_two_reward_All);

data_table.n_minus_three_reward_All =  data_table.n_minus_three_Durations_All >=GCAMP.Criteria;
data_table.n_minus_three_reward_All =categorical(data_table.n_minus_three_reward_All);

data_table.moving_mean = moving_average_lp_length;
data_table.moving_average_lp_length_n7andback = moving_average_lp_length_n7andback;

data_table.total_reward_indicator = total_reward_indicator;
logical_rewards = Durations >=GCAMP.Criteria;
criteria_percent = length(Durations)/sum(logical_rewards);
criteria_percent_indicator =  ones(length(Durations),1)*criteria_percent;
data_table.criteria_percent_indicator = criteria_percent_indicator;
%the he idx is tricky. Need to paste all the LP, HE and RE timestamps into a
%single array (keeping a variable to identify them) and then sort based on
%time
LP_start_times = GCAMP.LP_ON_timestamps;
LP_start_times =[ LP_start_times repmat(10,[length(LP_start_times) 1])  ];
 
HE_start_times = GCAMP.HE_ON_timestamps;
HE_start_times =[ HE_start_times repmat(12,[length(HE_start_times) 1])  ];

RE_start_times = GCAMP.RE_ON_timestamps;
RE_start_times =[ RE_start_times repmat(20,[length(RE_start_times) 1])  ];

Event_start_times = [LP_start_times; HE_start_times; RE_start_times];
[spot og_idx ] = sort(Event_start_times(:,1));
sorted_Events = Event_start_times(og_idx,:);
GCAMP.sorted_Events = sorted_Events;
%now loop through to find when a lp 
nanpad = [NaN NaN; NaN NaN];
sorted_Events = [sorted_Events; nanpad];
HE_Indicator =[];
HE_counter_rewarded = 0;
HE_counter_general =0;
HE_counter_unrewarded =0;
for he_it = 1:length(sorted_Events)
        
        
if sorted_Events(he_it,2) == 10 && sorted_Events(he_it+1,2) ==20  %if lp and rewarded
    
    if sorted_Events(he_it+2,2) == 12 %if theres a 12 after the 20 add it to the general he counter and reward specific
        HE_counter_rewarded = HE_counter_rewarded +1;
        HE_counter_general = HE_counter_general +1;
        
        HE_Indicator = [HE_Indicator; 1];
        
    elseif sorted_Events(he_it+2,2) ~= 12
             HE_Indicator = [HE_Indicator; 0];
    end
end

if sorted_Events(he_it,2) == 10 && sorted_Events(he_it+1,2) ~=20  %if LP was unrewarded
    if sorted_Events(he_it+1,2) == 12 %if theres a 12 after the 15 add it to the general he counter and unreward specific
        HE_counter_unrewarded = HE_counter_unrewarded +1;
        HE_counter_general = HE_counter_general +1;
        
        HE_Indicator = [HE_Indicator; 1];
        
    elseif sorted_Events(he_it+1,2) ~= 12
             HE_Indicator = [HE_Indicator; 0];
    end
             
end
end
HE_n1_Indicator = [NaN; HE_Indicator];
HE_n1_Indicator(end) =[];
HE_n2_Indicator = [NaN; NaN; HE_Indicator];
HE_n2_Indicator(end-1:end) =[];
HE_n3_Indicator = [NaN; NaN; NaN; HE_Indicator];
HE_n3_Indicator(end-2:end) =[];

data_table.HE_indicator = HE_Indicator;
data_table.HE_n1_Indicator = HE_n1_Indicator;
data_table.HE_n2_Indicator = HE_n2_Indicator;
data_table.HE_n3_Indicator = HE_n3_Indicator;
data_table.HE_indicator = categorical(data_table.HE_indicator);
data_table.HE_n1_Indicator = categorical(data_table.HE_n1_Indicator);
data_table.HE_n2_Indicator = HE_n2_Indicator;
data_table.HE_n2_Indicator = categorical(data_table.HE_n2_Indicator);
data_table.HE_n3_Indicator = categorical(data_table.HE_n3_Indicator);
%we have lp timestamps, but they are in bonsai time. Basically just need to
%make them in relation to the session, so timestamp 1 will become 0, and
%then timed up from there
lp_on_times  = GCAMP.LP_ON_timestamps;
lp_on_times = lp_on_times- lp_on_times(1);
data_table.lp_on_times = lp_on_times;

%up state n-1 index
 [Up_State_idx,ilower,uppersum,lowersum] = cusum(Durations,2,1,mean(Durations),std(Durations),'all');

        Up_State_Lengths = Durations(Up_State_idx); % Lever Press Lengths in upstate
        
        %need to get a column of 1 and 0 to indicate up-state or not
        idxxx=[];
        for state_it = 1:length(Up_State_idx)
          idxxx=[idxxx  press_indices==Up_State_idx(state_it)];
        end
        Up_State_idx_Logical=any(idxxx,2);
        
        %get the state identity of n-1, to be used to see if the state
        %identy of n-1 interacts with n:n-1
        Up_State_idx_n1= [NaN; Up_State_idx_Logical];
        Up_State_idx_n1=Up_State_idx_n1(1:end-1);

data_table.Up_State_idx_n1 = Up_State_idx_n1;
data_table.Up_State_idx_n1 =categorical(data_table.Up_State_idx_n1 );

%% Code to shuffle the order x number of times in order to build a distribution of shuffled data.
%also need to build the n-back arrays from these shuffled datasets
Shuffled_Durs_Distribtuion = [];
shuffled_ipi1_Distribution = [];
shuffled_he1_Distribution =[];
shuffled_he2_Distribution =[];
shuffled_MA_Distribution =[];
Shuffled_ts_Distribution =[];
shuffled_MA_n7_Distribution =[];
shuffled_total_reward_indicator_dist =[];
shuf_n_back_Lengths_distribution = {};
shuf_Logical_n_back_Lengths_distribution = {};
shuf_n_back_moving_averages ={};

for shuf_it = 1:shuff_num
    Shuffled_Durs = Durations;
    Shuffled_Durs = Shuffled_Durs(randperm(length(Shuffled_Durs)));
    Shuffled_Durs_Distribtuion = [Shuffled_Durs_Distribtuion Shuffled_Durs]; 
          
    Shuffled_ts = lp_on_times;
    Shuffled_ts = Shuffled_ts(randperm(length(Shuffled_ts)));
    Shuffled_ts_Distribution = [Shuffled_ts_Distribution Shuffled_ts];
       
    %calculate a moving average for each shuffle
    shuffled_moving_average_lp_length = movmean(Shuffled_Durs,[ma_length,0]);      
    shuffled_moving_average_lp_length = [NaN; NaN; shuffled_moving_average_lp_length];  
    shuffled_moving_average_lp_length = shuffled_moving_average_lp_length(1:end-2);
     shuffled_MA_Distribution = [shuffled_MA_Distribution shuffled_moving_average_lp_length];
     
        shuffled_moving_average_lp_length_n7andback = [NaN; NaN; NaN; NaN; NaN; shuffled_moving_average_lp_length];
    shuffled_moving_average_lp_length_n7andback = shuffled_moving_average_lp_length_n7andback(1:end-5);
     shuffled_MA_n7_Distribution = [shuffled_MA_n7_Distribution shuffled_moving_average_lp_length_n7andback];
  
     shuffled_total_reward_indicator = total_reward_indicator;
shuffled_total_reward_indicator = shuffled_total_reward_indicator(randperm(length(shuffled_total_reward_indicator)));
 shuffled_total_reward_indicator_dist = [shuffled_total_reward_indicator_dist shuffled_total_reward_indicator];    
     
    %independentyly shuffle ipi and he
    shuf_ipi1= ipi1;
    shuffled_ipi1 = shuf_ipi1(randperm(length(shuf_ipi1)));
    %removed the original nan, now add it back in
    shuffled_ipi1 = [NaN; shuffled_ipi1];
    shuffled_ipi1 = shuffled_ipi1(1:end-1);
    shuffled_ipi1_Distribution=[shuffled_ipi1_Distribution  shuffled_ipi1];  
     
    shuf_he1 = HE_n1_Indicator(2:end);
    shuf_he1 =  shuf_he1(randperm(length(shuf_he1)));
    shuf_he1 = [NaN; shuf_he1];
    shuffled_he1_Distribution =[shuffled_he1_Distribution shuf_he1];
    shuf_he2 = HE_n2_Indicator(3:end);
    shuf_he2 =  shuf_he2(randperm(length(shuf_he2)));
    shuf_he2 = [NaN; NaN; shuf_he2];
    shuffled_he2_Distribution =[shuffled_he2_Distribution shuf_he2];
    
end

GCAMP.shuffled_MA = shuffled_MA_Distribution;
GCAMP.shuffled_ipi1_Distribution = shuffled_ipi1_Distribution;
GCAMP.shuffled_he1_Distribution =shuffled_he1_Distribution;
GCAMP.shuffled_he2_Distribution =shuffled_he2_Distribution;
GCAMP.Shuffled_ts_Distribution =Shuffled_ts_Distribution;
GCAMP.shuffled_MA_n7_Distribution =shuffled_MA_n7_Distribution;
GCAMP.shuffled_total_reward_indicator_dist =shuffled_total_reward_indicator_dist;
 
%now shuffle the nback lengths  
for shuffle_length = 1:size(Shuffled_Durs_Distribtuion,2)
    shuf_n_back_Lengths = [];
    shuf_n_lengths_array = [];
    shuf_Logical_n_lengths_array= [];
    shuf_Logical_n_back_Lengths =[];
   
           current_shuffle_lps = Shuffled_Durs_Distribtuion(:,shuffle_length);
           shuf_press_indices = 1:length(current_shuffle_lps);
           shuf_press_indices = shuf_press_indices';
           
           %for n-backs 1 to 10
       for x = 1:10 
   
       for n_it = 1:length(current_shuffle_lps)
           %need to check that if there are at least x presses before the
           %current state press.
           if  press_indices(n_it) <= x
           %if there are not, then put a NaN there
           n_back = NaN;
           elseif  press_indices(n_it) > x
           n_back = current_shuffle_lps(n_it- x);               
           end
           shuf_n_lengths_array = [shuf_n_lengths_array n_back];
           
           %add in logical lengths 0/1 for failure/success
          
           shuf_Logical_n_lengths_array = shuf_n_lengths_array >= GCAMP.Criteria;
           
           end
       % shuf_n_lengths_array
      shuf_n_back_Lengths = [shuf_n_back_Lengths; shuf_n_lengths_array]; 
      shuf_n_lengths_array =[];
      
      %the logical test converts all NaNs to 0, need to put them back so we
      %don't inappropriately put lever press prior to the start of the
      %session
               shuf_Logical_n_lengths_array = double(shuf_Logical_n_lengths_array);
        shuf_Logical_n_lengths_array(1:x) = NaN;
         shuf_Logical_n_back_Lengths = [shuf_Logical_n_back_Lengths; shuf_Logical_n_lengths_array];  
         shuf_Logical_n_lengths_array =[];
   
       end
       shuf_n_back_Lengths_distribution = [shuf_n_back_Lengths_distribution; shuf_n_back_Lengths];
       shuf_Logical_n_back_Lengths_distribution = [shuf_Logical_n_back_Lengths_distribution; shuf_Logical_n_back_Lengths];    
       
end
% GCAMP.shuf_n1_reward_Distribution = shuf_Logical_n_back_Lengths_distribution
GCAMP.shuf_n_back_Lengths_distribution = shuf_n_back_Lengths_distribution;
GCAMP.Shuffled_Durs_Distribtuion = Shuffled_Durs_Distribtuion;
GCAMP.shuf_Logical_n_back_Lengths_distribution =shuf_Logical_n_back_Lengths_distribution;

%% get gcamp data for relevant events. Segment into timepoints of interest and mean and AUC
%see if baseline, baseline_raw (no running filter), or z-scored raw data
%affects the direction/magnitude of any effects.

%baseline z-score, either raw, or with the running filter
Onset_baselinezscore_Raw = GCAMP.LP_ON;
Offset_baselinezscore_Raw = GCAMP.LP_OFF;
Onset_baselinezscore =  GCAMP.baseline_norm_LP_ON;
Offset_baselinezscore =  GCAMP.baseline_norm_LP_OFF;
%lp offset for rewarded lp is the same is reward
RE_baselinezscore_raw = GCAMP.LP_OFF_Met;
RE_baselinezscore = GCAMP.baseline_norm_LP_OFF_Met;

%get the baselines themselves
baselines_on = GCAMP.baseline;
baselines_on_raw = GCAMP.baseline;

%based off baseline with running filter data
interpolated_base = GCAMP.baseline_norm_Duration;

%now we want to take several periods around press onset and get both AVG
%and AUC
 event_onset_idx = 1+abs(GCAMP.SR*base_time_end);
% Make Table for SVM
LP_ON_window = [find(base_time_end == GCAMP.plot_time) : find(0 == GCAMP.plot_time)];
LP_OFF_window = [find(0 == GCAMP.plot_time) : find(time_end == GCAMP.plot_time)];
%% lp onset
base_neg_1_to_onset = Onset_baselinezscore(:,LP_ON_window);
base_neg_1_to_onset_AUC_raw = trapz(base_neg_1_to_onset,2);
base_neg_1_to_onset_Mean_raw = mean(base_neg_1_to_onset,2);
base_neg_1_to_onset_Slope_raw= (mean(base_neg_1_to_onset(:,1:4),2)-mean(base_neg_1_to_onset(:,end-4:end),2))/ (17);
data_table.base_neg_1_to_onset_AUC_raw = base_neg_1_to_onset_AUC_raw ;
data_table.base_neg_1_to_onset_Mean_raw = base_neg_1_to_onset_Mean_raw ;
data_table.braw_neg_1_to_onset_Slope_raw= base_neg_1_to_onset_Slope_raw;
%% lp offset
base_neg_1_to_offset = Offset_baselinezscore(:,LP_OFF_window);
base_offset_to1_AUC_raw = trapz(base_neg_1_to_offset,2);
base_offset_to1_Mean_raw = mean(base_neg_1_to_offset,2);
base_neg_1_to_offset_Slope_raw= (mean(base_neg_1_to_offset(:,1:4),2)-mean(base_neg_1_to_offset(:,end-4:end),2))/ (17);
data_table.base_offset_to1_AUC_raw = base_offset_to1_AUC_raw;
data_table.base_offset_to1_Mean_raw = base_offset_to1_Mean_raw;
data_table.base_neg_1_to_offset_Slope_raw = base_neg_1_to_offset_Slope_raw;
%% reward on/offset
%%interpolated presses
%divide the interpolated presses into 4 bins 
interp_quarter = size(interpolated_base,2)/4;
interp_first_quart = interpolated_base(:,1:interp_quarter);
interp_second_quart = interpolated_base(:,interp_quarter+1:interp_quarter*2);
interp_third_quart = interpolated_base(:,1+interp_quarter*2:interp_quarter*3);
interp_fourth_quart = interpolated_base(:,1+interp_quarter*3:end);
%get mean and AUC
interp_all_AUC = trapz(interpolated_base,2);
interp_all_mean = mean(interpolated_base,2);
interp_firstquart_AUC = trapz(interp_first_quart,2);
interp_firstquart_mean =mean(interp_first_quart,2);
interp_secondquart_AUC = trapz(interp_second_quart,2);
interp_secondquart_mean =mean(interp_second_quart,2);
interp_thirdquart_AUC = trapz(interp_third_quart,2);
interp_thirdquart_mean =mean(interp_third_quart,2);
interp_fourthquart_AUC = trapz(interp_fourth_quart,2);
interp_fourthquart_mean =mean(interp_fourth_quart,2);
interp_all_Slope = (nanmean(interpolated_base(:,1:4),2)-nanmean(interpolated_base(:,end-4:end),2))/ (16);
data_table.interp_all_AUC = interp_all_AUC;
data_table.interp_all_mean =interp_all_mean;
data_table.interp_firstquart_AUC = interp_firstquart_AUC;
data_table.interp_firstquart_mean =interp_firstquart_mean;
data_table.interp_secondquart_AUC = interp_secondquart_AUC;
data_table.interp_secondquart_mean =interp_secondquart_mean;
data_table.interp_thirdquart_AUC = interp_thirdquart_AUC;
data_table.interp_thirdquart_mean =interp_thirdquart_mean;
data_table.interp_fourthquart_AUC = interp_fourthquart_AUC;
data_table.interp_fourthquart_mean =interp_fourthquart_mean;
data_table.interp_all_Slope = interp_all_Slope;
%% baselines
baselines_mean = mean(baselines_on,2);
baselines_AUC = trapz(baselines_on,2);
baselines_raw_mean = mean(baselines_on_raw,2);
baselines_raw_AUC = trapz(baselines_on_raw,2);
data_table.baselines_mean = baselines_mean;
data_table.baselines_AUC = baselines_AUC;
data_table.baselines_raw_mean = baselines_raw_mean;
data_table.baselines_raw_AUC = baselines_raw_AUC;
%% add in lagged variables
data_table.interp_all_AUC_n1 = [NaN; data_table.interp_all_AUC(1:end-1)]; 
data_table.interp_all_AUC_n2 = [NaN; data_table.interp_all_AUC_n1(1:end-1)]; 
data_table.interp_all_AUC_n3 = [NaN; data_table.interp_all_AUC_n2(1:end-1)]; 
data_table.interp_all_AUC_n4 = [NaN; data_table.interp_all_AUC_n3(1:end-1)]; 
data_table.base_neg_1_to_onset_AUC_raw_n1 = [NaN; data_table.base_neg_1_to_onset_AUC_raw(1:end-1)]; 
data_table.base_neg_1_to_onset_AUC_raw_n2 = [NaN; data_table.base_neg_1_to_onset_AUC_raw_n1(1:end-1)]; 
data_table.base_neg_1_to_onset_AUC_raw_n3 = [NaN; data_table.base_neg_1_to_onset_AUC_raw_n2(1:end-1)]; 
data_table.base_neg_1_to_onset_AUC_raw_n4 = [NaN; data_table.base_neg_1_to_onset_AUC_raw_n3(1:end-1)]; 
data_table.base_offset_to1_AUC_raw_n1 = [NaN; data_table.base_offset_to1_AUC_raw(1:end-1)];
data_table.base_offset_to1_AUC_raw_n2 = [NaN; data_table.base_offset_to1_AUC_raw_n1(1:end-1)];
data_table.base_offset_to1_AUC_raw_n3 = [NaN; data_table.base_offset_to1_AUC_raw_n2(1:end-1)];
data_table.base_offset_to1_AUC_raw_n4 = [NaN; data_table.base_offset_to1_AUC_raw_n3(1:end-1)];
%% Save
%now we need to save stuff into the GCAMP structure
%in particular need to save the table structures, and add in indicator
%variables for mouse and day of training
reps = size(data_table,1);
mouse_indicator = repmat(GCAMP.onlyID ,reps,1);
day_indicator = repmat(GCAMP.training_day,reps,1);
criteria_indicator = repmat(GCAMP.Criteria,reps,1);
data_table.mouse_indicator = mouse_indicator;
data_table.day_indicator = day_indicator;
data_table.criteria_indicator = criteria_indicator;
data_table_variables = data_table.Properties.VariableNames;
GCAMP.data_table_variables = data_table_variables;
%need to convert back to an array so that we can past all thetables
%together 
regression_cell= table2cell(data_table);
GCAMP.regression_cell = regression_cell;
end



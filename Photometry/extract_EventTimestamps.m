function [GCAMP] = extract_EventTimestamps(GCAMP)
%% Saves GCAMP and Behavior data into a data structure
%% Extract behavior timestamps from analog_inputs
if ischar(GCAMP.Criteria)
    GCAMP.Criteria = str2num(GCAMP.Criteria);
end
% Subtract logicals from each other in Lever Press column
switches_LP = diff(GCAMP.beh_data(:,2)); %2 %OFC paper
% Subtract logicals from each other in Reinforcement column
switches_RE = diff(GCAMP.beh_data(:,4)); %4 %OFC paper
% Subtract logicals from each other in HeadEntries column
switches_HE = diff(GCAMP.beh_data(:,3)); %3 %OFC paper

% When the difference is 1, that's event ON
LP_ON_index = find(switches_LP == 1)+1;
RE_ON_index = find(switches_RE == 1)+1;
HE_ON_index = find(switches_HE == 1)+1;

% When the difference is -1, that's event OFF
LP_OFF_index = find(switches_LP == -1)+1;
RE_OFF_index = find(switches_RE == -1)+1;
HE_OFF_index = find(switches_HE == -1)+1;

% Get the time stamps from those event indices
LP_ON_timestamps = GCAMP.beh_data(LP_ON_index, 1);
LP_OFF_timestamps = GCAMP.beh_data(LP_OFF_index, 1);

%code to exclude the first few lp if a lp happens too early during baseline 
%LP_ON_timestamps = LP_ON_timestamps(3:end); %OFC paper
%LP_OFF_timestamps = LP_OFF_timestamps(3:end); %OFC paper
RE_ON_timestamps = GCAMP.beh_data(RE_ON_index, 1);
RE_OFF_timestamps = GCAMP.beh_data(RE_OFF_index, 1);
HE_ON_timestamps = GCAMP.beh_data(HE_ON_index, 1);
HE_OFF_timestamps = GCAMP.beh_data(HE_OFF_index,1);

% Matrix of times where Column 1 is event ON, Column 2 is event OFF
LP_All = [LP_ON_timestamps LP_OFF_timestamps];
RE_All = [RE_ON_timestamps RE_OFF_timestamps];
HE_All = [HE_ON_timestamps HE_OFF_timestamps];

% Hold down time for each lever press
HD_time=(LP_All(:,2)-LP_All(:,1))-20;
HE_time =(HE_All(:,2)-HE_All(:,1))-20;

% Lever presses that met criteria and were reinforced
Criteria_met = find(HD_time >=  GCAMP.Criteria); %GCAMP.Criteria);
Criteria_fail = find(HD_time <   GCAMP.Criteria);
Total_Reinforcers = length(Criteria_met);

%HEs that were above or below the mean
HE_Above_median = find(HE_time >= median(HE_time));
HE_Below_median = find(HE_time < median(HE_time));
HE_idx = find(HE_time > 0);

% For peri-event data
SR = 20; % sampling rate

%Find the first HE made after a Reinforcer
First_HE_After_RE_idx =[];
for k = 1:length(RE_ON_timestamps) %for each RE_ON timestamp
    % Find index of this LP_ON timestamp in the vector of all of your
    % timestamps
    Closest_HE_After_RE_idx = nearestpoint(RE_ON_timestamps(k),HE_ON_timestamps(HE_idx),'next');
    First_HE_After_RE_idx = [First_HE_After_RE_idx; Closest_HE_After_RE_idx];
    % First_HE_After_RE_idx(isnan(First_HE_After_RE_idx)) = [];
end
% Remove NaNs (HE happening outside of last LP)
First_HE_After_RE_idx(isnan(First_HE_After_RE_idx)) = [];
% Remove instances in which the same head entry happens after two or more
% reinforcements to avoid sampling duplicate traces
First_HE_After_RE_idx = unique(First_HE_After_RE_idx,'rows');
% Remove rare instances in which nearestpoint captures a head entry before
% a lever press
remove_idx = find(HE_ON_timestamps(First_HE_After_RE_idx) <= min(LP_ON_timestamps));
First_HE_After_RE_idx(remove_idx) = [];

%Find the first HE made after a successful LP (useful for probability 
% experiments where some successful LPs are not rewarded)
First_HE_After_Success_idx  = [];
for k = 1:length(Criteria_met) %for each successful LP_ON timestamp
    % Find index of this LP_ON timestamp in the vector of all of your
    % timestamps
    Closest_HE_After_Success_idx = nearestpoint(LP_ON_timestamps(k),HE_ON_timestamps(HE_idx), 'next');
    First_HE_After_Success_idx = [First_HE_After_Success_idx; Closest_HE_After_Success_idx];
    % First_HE_After_RE_idx(isnan(First_HE_After_RE_idx)) = [];
end
% Remove NaNs (HE happening outside of last LP)
First_HE_After_Success_idx(isnan(First_HE_After_Success_idx)) = [];
% Remove instances in which the same head entry happens after two or more
% successful lever presses to avoid sampling duplicate traces
First_HE_After_Success_idx = unique(First_HE_After_Success_idx,'rows');
% Get index of head entries made after a successful lever press in which
% reward was not delivered (used for probability 
% experiments where some successful LPs are not rewarded)
First_HE_After_Success_noRE_idx = First_HE_After_Success_idx;
[a,b] = ismember(First_HE_After_RE_idx, First_HE_After_Success_noRE_idx);
 b(b==0) = [];
First_HE_After_Success_noRE_idx(b) = [];
% Remove rare instances in which nearestpoint captures a head entry before
% a lever press
remove_idx = find(HE_ON_timestamps(First_HE_After_Success_noRE_idx) <= min(LP_ON_timestamps));
First_HE_After_Success_noRE_idx(remove_idx) = [];


%Find the first HE made after a failed LP
First_HE_After_Failure_idx  = [];
for k = 1:length(Criteria_fail) %for each failed LP_ON timestamp
    % Find index of this LP_ON timestamp in the vector of all of your
    % timestamps
    Closest_HE_After_Failure_idx = nearestpoint(LP_ON_timestamps(k),HE_ON_timestamps(HE_idx), 'next');
    First_HE_After_Failure_idx = [First_HE_After_Failure_idx; Closest_HE_After_Failure_idx];
    % First_HE_After_RE_idx(isnan(First_HE_After_RE_idx)) = [];
end
% Remove NaNs (HE happening outside of last LP)
First_HE_After_Failure_idx(isnan(First_HE_After_Failure_idx)) = [];
% Remove instances in which the same head entry happens after two or more
% lever presses to avoid sampling duplicate traces.
First_HE_After_Failure_idx = unique(First_HE_After_Failure_idx,'rows');
% Get index of head entries made after a failed lever press in which
% reward was not delivered or a success was made beforehand (used for probability 
% experiments where some successful LPs are not rewarded). This essentially
% avoids including head entries made after a success or reward delivery in
% which failed lever presses happened in between.
First_HE_After_Failure_no_RE_idx = First_HE_After_Failure_idx;
[a,b] = ismember(First_HE_After_RE_idx, First_HE_After_Failure_no_RE_idx);
 b(b==0) = [];
First_HE_After_Failure_no_RE_idx(b) = [];
[a,b] = ismember(First_HE_After_Success_noRE_idx, First_HE_After_Failure_no_RE_idx);
 b(b==0) = [];
First_HE_After_Failure_no_RE_idx(b) = [];
% Remove rare instances in which nearestpoint captures a head entry before
% a lever press
remove_idx = find(HE_ON_timestamps(First_HE_After_Failure_no_RE_idx) <= min(LP_ON_timestamps));
First_HE_After_Failure_no_RE_idx(remove_idx) = [];


%% Save Data
% Create GCAMP data structure
GCAMP.HoldDown_times = HD_time;
GCAMP.Total_Reinforcers = Total_Reinforcers;
GCAMP.SR = SR;
GCAMP.Criteria_met = Criteria_met;
GCAMP.Criteria_fail = Criteria_fail;
GCAMP.LP_ON_timestamps = LP_ON_timestamps;
GCAMP.RE_ON_timestamps = RE_ON_timestamps;
GCAMP.LP_OFF_timestamps = LP_OFF_timestamps;
GCAMP.HE_ON_timestamps = HE_ON_timestamps;
GCAMP.HE_OFF_timestamps = HE_OFF_timestamps;
GCAMP.HE_Times = HE_time;
GCAMP.Short_HE = HE_Below_median;
GCAMP.Long_HE = HE_Above_median;
GCAMP.First_HE_After_RE = First_HE_After_RE_idx;
GCAMP.First_HE_After_Met_noRE = First_HE_After_Success_noRE_idx;
GCAMP.First_HE_After_Fail = First_HE_After_Failure_no_RE_idx;
GCAMP.HE_idx = HE_idx;
end


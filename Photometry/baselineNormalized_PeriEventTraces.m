function [GCAMP] = baselineNormalized_PeriEventTraces(GCAMP)
%% Express data as baseline change
%% Set up some basic variables
GCAMP.plot_time = GCAMP.base_time_end:1/GCAMP.SR:GCAMP.time_end; % Time series in second for per-event GCAMP traces
GCAMP.interp_length = 10; % Length of interpolated lever press duration GCAMP traces
% First Head Entry After Reward
% Get timestamps of the first head entry after true reward or no reward
First_HE_After_RE_timestamps = GCAMP.HE_ON_timestamps(GCAMP.First_HE_After_RE);
First_HE_After_Fail_timestamps = GCAMP.HE_ON_timestamps(GCAMP.First_HE_After_Fail);
% Remove head entries that occurred before lever presses
%% Pre-allocate data matrices
% Lever Press Onset, Head entries and Lever press Offset
delta_F_LP_ON = zeros(length(GCAMP.base_time_end:1/GCAMP.SR:GCAMP.time_end),...
     length(GCAMP.HoldDown_times))';
delta_F_LP_OFF = zeros(length(GCAMP.base_time_end:1/GCAMP.SR:GCAMP.time_end),...
     length(GCAMP.HoldDown_times))';
delta_F_Pre_LP_ON_Baselines = zeros(length(GCAMP.base_time_start:1/GCAMP.SR:GCAMP.base_time_end),...
     length(GCAMP.HoldDown_times))'; 
delta_F_Pre_LP_ON_Baselines_Mean = zeros(1,...
     length(GCAMP.HoldDown_times))';
delta_F_Pre_LP_ON_Baselines_STD = zeros(1,...
     length(GCAMP.HoldDown_times))';
 
delta_F_First_HE_After_RE = zeros(length(GCAMP.base_time_end:1/GCAMP.SR:GCAMP.time_end),...
     length(First_HE_After_RE_timestamps))'; 
delta_F_Pre_First_HE_After_RE_Baselines_Idx = zeros(1,...
     length(First_HE_After_RE_timestamps))'; 
 
delta_F_First_HE_After_Fail = zeros(length(GCAMP.base_time_end:1/GCAMP.SR:GCAMP.time_end),...
     length(First_HE_After_Fail_timestamps))';
delta_F_First_HE_After_Fail_Baselines_Idx = zeros(1,...
     length(First_HE_After_Fail_timestamps))'; 
 
%% For each Lever Press, extract peri-event (onset and offset) and pre-lever press onset baseline GCAMP traces
for lever_press = 1:length(GCAMP.LP_ON_timestamps)
    % Find indices of lever press onset and offset timestamps closest to GCAMP trace
    % timestamps
    Closest_LP_ON_idx = nearestpoint(GCAMP.LP_ON_timestamps(lever_press),GCAMP.gcampdata_timestamps);
    Closest_LP_OFF_idx = nearestpoint(GCAMP.LP_OFF_timestamps(lever_press),GCAMP.gcampdata_timestamps);
    % Extract peri-event traces for lever press onset and offset
    data_window_LP_ON = GCAMP.gcampdata(Closest_LP_ON_idx + GCAMP.base_time_end * GCAMP.SR : Closest_LP_ON_idx + GCAMP.time_end * GCAMP.SR)';
    data_window_LP_OFF = GCAMP.gcampdata(Closest_LP_OFF_idx + GCAMP.base_time_end * GCAMP.SR : Closest_LP_OFF_idx + GCAMP.time_end * GCAMP.SR)';
    delta_F_LP_ON(lever_press,:) = data_window_LP_ON;
    delta_F_LP_OFF(lever_press,:) = data_window_LP_OFF;
    % Extract pre-lever press onset baseline traces
    Pre_LP_ON_baseline = GCAMP.gcampdata(Closest_LP_ON_idx + GCAMP.base_time_start * GCAMP.SR : Closest_LP_ON_idx + GCAMP.base_time_end * GCAMP.SR);
    delta_F_Pre_LP_ON_Baselines(lever_press,:) = Pre_LP_ON_baseline;
    delta_F_Pre_LP_ON_Baselines_Mean(lever_press) = mean(Pre_LP_ON_baseline);
    delta_F_Pre_LP_ON_Baselines_STD(lever_press) = std(Pre_LP_ON_baseline);
end
% Normalize peri-event GCAMP traces to pre-lever press onset baselines
% Note: for the offset, we are actually going to use the baseline for lp onset.
% That way we will use the same baseline for onset/offset for a given lp
GCAMP.baseline_norm_LP_ON = (delta_F_LP_ON - delta_F_Pre_LP_ON_Baselines_Mean) ./ delta_F_Pre_LP_ON_Baselines_STD;
GCAMP.baseline_norm_LP_ON_Met = GCAMP.baseline_norm_LP_ON(GCAMP.Criteria_met, :);
GCAMP.baseline_norm_LP_ON_Fail = GCAMP.baseline_norm_LP_ON(GCAMP.Criteria_fail, :);
% Raw peri-event GCAMP traces
GCAMP.LP_ON = delta_F_LP_ON;
GCAMP.LP_ON_Met = GCAMP.LP_ON(GCAMP.Criteria_met, :);
GCAMP.LP_ON_Fail = GCAMP.LP_ON(GCAMP.Criteria_fail, :);
% Normalize peri-event GCAMP traces to pre-lever press onset baselines
% Note: for the offset, we are actually going to use the baseline for lp onset.
% That way we will use the same baseline for onset/offset for a given lp
GCAMP.baseline_norm_LP_OFF = (delta_F_LP_OFF - delta_F_Pre_LP_ON_Baselines_Mean)./ delta_F_Pre_LP_ON_Baselines_STD;
GCAMP.baseline_norm_LP_OFF_Met = GCAMP.baseline_norm_LP_OFF(GCAMP.Criteria_met,:);
GCAMP.baseline_norm_LP_OFF_Fail = GCAMP.baseline_norm_LP_OFF(GCAMP.Criteria_fail,:);
% Raw peri-event GCAMP traces
GCAMP.LP_OFF = delta_F_LP_OFF;
GCAMP.LP_OFF_Met = GCAMP.LP_OFF(GCAMP.Criteria_met, :);
GCAMP.LP_OFF_Fail = GCAMP.LP_OFF(GCAMP.Criteria_fail, :);
% Baselines
GCAMP.baseline = delta_F_Pre_LP_ON_Baselines;
    
%% For each Reward delivery, extract peri-event (first head entry after reward) GCAMP traces
    
for head_entry = 1:length(First_HE_After_RE_timestamps)
    % Find indices of first head entry after reward timestamps closes to GCAMP trace
    % timestamps
    Closest_HE_ON_idx = nearestpoint(First_HE_After_RE_timestamps(head_entry),GCAMP.gcampdata_timestamps);
    % Find index of the lever press that led to this head entry
    Closest_LP_before_HE_ON_idx = nearestpoint(First_HE_After_RE_timestamps(head_entry),GCAMP.LP_ON_timestamps,'previous');
    LP_before_HE_ON_timestamp = GCAMP.LP_ON_timestamps(Closest_LP_before_HE_ON_idx);
    LP_Idx_before_HE_ON_timestamp = nearestpoint(LP_before_HE_ON_timestamp,GCAMP.gcampdata_timestamps);
    % Extract peri-event traces for first head entry after reward
    data_window_HE_ON = GCAMP.gcampdata(Closest_HE_ON_idx + GCAMP.base_time_end * GCAMP.SR : Closest_HE_ON_idx + GCAMP.time_end * GCAMP.SR)';
    delta_F_First_HE_After_RE(head_entry,:) = data_window_HE_ON;
    % Extract pre-first head entry after reward baseline traces
    %Pre_First_HE_After_RE_baseline = GCAMP.gcampdata(LP_Idx_before_HE_ON_timestamp + GCAMP.base_time_start * GCAMP.SR : LP_Idx_before_HE_ON_timestamp + GCAMP.base_time_end * GCAMP.SR);
    delta_F_Pre_First_HE_After_RE_Baselines_Idx(head_entry,:) = Closest_LP_before_HE_ON_idx;

end
% Normalize peri-event GCAMP traces to pre-lever press onset baselines
GCAMP.baseline_norm_First_HE_After_RE = (delta_F_First_HE_After_RE - delta_F_Pre_LP_ON_Baselines_Mean(delta_F_Pre_First_HE_After_RE_Baselines_Idx)) ./ delta_F_Pre_LP_ON_Baselines_STD(delta_F_Pre_First_HE_After_RE_Baselines_Idx);


%% For each failed lever press, extract peri-event (first head entry after press) GCAMP traces
    
for head_entry = 1:length(First_HE_After_Fail_timestamps)
    % Find indices of first head entry after reward timestamps closes to GCAMP trace
    % timestamps
    Closest_HE_ON_idx = nearestpoint(First_HE_After_Fail_timestamps(head_entry),GCAMP.gcampdata_timestamps);
    % Find index of the lever press that led to this head entry
    Closest_LP_before_HE_ON_idx = nearestpoint(First_HE_After_Fail_timestamps(head_entry),GCAMP.LP_ON_timestamps,'previous');
    LP_before_HE_ON_timestamp = GCAMP.LP_ON_timestamps(Closest_LP_before_HE_ON_idx);
    LP_Idx_before_HE_ON_timestamp = nearestpoint(LP_before_HE_ON_timestamp,GCAMP.gcampdata_timestamps);
    % Extract peri-event traces for first head entry after reward
    data_window_HE_ON = GCAMP.gcampdata(Closest_HE_ON_idx + GCAMP.base_time_end * GCAMP.SR : Closest_HE_ON_idx + GCAMP.time_end * GCAMP.SR)';
    delta_F_First_HE_After_Fail(head_entry,:) = data_window_HE_ON;
    % Extract pre-first head entry after reward baseline traces
    %Pre_First_HE_After_RE_baseline = GCAMP.gcampdata(LP_Idx_before_HE_ON_timestamp + GCAMP.base_time_start * GCAMP.SR : LP_Idx_before_HE_ON_timestamp + GCAMP.base_time_end * GCAMP.SR);
    delta_F_First_HE_After_Fail_Baselines_Idx(head_entry,:) = Closest_LP_before_HE_ON_idx;

end
% Normalize peri-event GCAMP traces to pre-lever press onset baselines
GCAMP.baseline_norm_First_HE_After_Fail = (delta_F_First_HE_After_Fail - delta_F_Pre_LP_ON_Baselines_Mean(delta_F_First_HE_After_Fail_Baselines_Idx)) ./ delta_F_Pre_LP_ON_Baselines_STD(delta_F_First_HE_After_Fail_Baselines_Idx);



%% For each Lever Press, linearly interpolate the continuous activity during each ongoing lever press

    interp_data_norm = nan(length(GCAMP.HoldDown_times), GCAMP.interp_length); % placeholder matrix for norm interpolated data (lever press x length of interpolated trace)
    interp_data = nan(length(GCAMP.HoldDown_times), GCAMP.interp_length); % placeholder matrix for interpolated data (lever press x length of interpolated trace)
    zero_index = find(GCAMP.plot_time == 0); % find 0 second index in vector of peri-event lever press onset gcamp data
    max_length = zero_index + GCAMP.time_end * GCAMP.SR; % max vector length for interpolated 
    timeinsamples = ceil((GCAMP.HoldDown_times)/((1/GCAMP.SR)*1000)); % Length of each lever press duration rounded to next nearest sample

    for lever_press = 1:length(GCAMP.HoldDown_times) % for each lever press
        % Find indices of lever press onset  timestamps closest to GCAMP trace
        % timestamps
        Closest_LP_ON_idx = nearestpoint(GCAMP.LP_ON_timestamps(lever_press),GCAMP.gcampdata_timestamps);
        Closest_LP_OFF_idx = nearestpoint(GCAMP.LP_OFF_timestamps(lever_press),GCAMP.gcampdata_timestamps);
        % Extract GCAMP trace for activity between lever press onset and offset
        data_window_LP_Duration = GCAMP.gcampdata(Closest_LP_ON_idx : Closest_LP_OFF_idx)';
        delta_F_LP_Duration{lever_press} = data_window_LP_Duration;
        % Z-score normalize to pre-lever press onset window
        baseline_z_score_LP_Duration{lever_press} = (delta_F_LP_Duration{lever_press} - delta_F_Pre_LP_ON_Baselines_Mean(lever_press))./ delta_F_Pre_LP_ON_Baselines_STD(lever_press);
        % Get total samples in duration activity (corresponding to nearest time
        % point for lever press onset and offset detected by nearestpoing function above)
        total_samples_in_duration(lever_press,:) = length(baseline_z_score_LP_Duration{lever_press});
        if  total_samples_in_duration(lever_press,:) >= 2 % Interpolation requires at least two sample points of duration activity]
            duration_activity_norm = baseline_z_score_LP_Duration{lever_press}; % variable vector of norm F (determined by duration)
            duration_activity = delta_F_LP_Duration{lever_press}; % variable vector of F (determined by duration)
            t0 = linspace(1,total_samples_in_duration(lever_press,:),total_samples_in_duration(lever_press,:)); %original time vector
            t1 = linspace(1,total_samples_in_duration(lever_press,:),GCAMP.interp_length); % new time vector (specifying the time points at which you want to interpolate)
            interp_data_norm(lever_press,:) = interp1(t0,duration_activity_norm,t1,'makikma'); %norm data interpolated to fixed sample space
            interp_data(lever_press,:) = interp1(t0,duration_activity,t1,'makikma'); %data interpolated to fixed sample space
        end
    end

% Normalize peri-event GCAMP traces to pre-lever press onset baselines
GCAMP.baseline_norm_Duration = interp_data_norm;
GCAMP.baseline_norm_Duration_Met = GCAMP.baseline_norm_Duration(GCAMP.Criteria_met, :);
GCAMP.baseline_norm_Duration_Fail = GCAMP.baseline_norm_Duration(GCAMP.Criteria_fail, :);
% Raw peri-event GCAMP traces
GCAMP.Duration = interp_data;
GCAMP.Duration_Met = GCAMP.Duration(GCAMP.Criteria_met, :);
GCAMP.Duration_Fail = GCAMP.Duration(GCAMP.Criteria_fail, :);


end


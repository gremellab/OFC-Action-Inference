function [GroupedGCAMP] = nBackReward(GroupedGCAMP)
%% Create empty vectors for GCAMP traces
total_sessions = length(GroupedGCAMP.Mice);
Grouped_RE_n1_LP_ON = [];
Grouped_noRE_n1_LP_ON = [];
Grouped_RE_n1_LP_OFF = [];
Grouped_noRE_n1_LP_OFF = [];
Grouped_RE_n1_LPInterp = [];
Grouped_noRE_n1_LPInterp = [];


Grouped_n0_RE_RE_n1_LP_ON = [];
Grouped_n0_RE_noRE_n1_LP_ON = [];
Grouped_n0_RE_RE_n1_LP_OFF = [];
Grouped_n0_RE_noRE_n1_LP_OFF = [];
Grouped_n0_RE_RE_n1_LPInterp = [];
Grouped_n0_RE_noRE_n1_LPInterp = [];

Grouped_n0_noRE_RE_n1_LP_ON = [];
Grouped_n0_noRE_noRE_n1_LP_ON = [];
Grouped_n0_noRE_RE_n1_LP_OFF = [];
Grouped_n0_noRE_noRE_n1_LP_OFF = [];
Grouped_n0_noRE_RE_n1_LPInterp = [];
Grouped_n0_noRE_noRE_n1_LPInterp = [];


for session = 1:total_sessions
    %% Find lever press durations for n-0 and n-1 lever presses
    n0_durations = cell2mat(GroupedGCAMP.Mice{session}.regression_cell(:,1));
    n1_durations = cell2mat(GroupedGCAMP.Mice{session}.regression_cell(:,2));
    % Was n-1 rewarded?
    RE_n1_indices = find(n1_durations >= GroupedGCAMP.Mice{session}.Criteria);
    noRE_n1_indices = find(n1_durations < GroupedGCAMP.Mice{session}.Criteria);
    % Was n-0 rewarded?
    RE_n0_indices = find(n0_durations >= GroupedGCAMP.Mice{session}.Criteria);
    noRE_n0_indices = find(n0_durations < GroupedGCAMP.Mice{session}.Criteria);
    % If n-0 was rewaded, was n-1 rewarded?
    RE_n0_RE_n1_indices = find((n0_durations >= GroupedGCAMP.Mice{session}.Criteria) & (n1_durations >= GroupedGCAMP.Mice{session}.Criteria));
    RE_n0_noRE_n1_indices = find((n0_durations >= GroupedGCAMP.Mice{session}.Criteria) & (n1_durations < GroupedGCAMP.Mice{session}.Criteria));
    % If n-0 was not rewaded, was n-1 rewarded?
    noRE_n0_RE_n1_indices = find((n0_durations < GroupedGCAMP.Mice{session}.Criteria) & (n1_durations >= GroupedGCAMP.Mice{session}.Criteria));
    noRE_n0_noRE_n1_indices = find((n0_durations < GroupedGCAMP.Mice{session}.Criteria) & (n1_durations < GroupedGCAMP.Mice{session}.Criteria));
    
    %% Split GCAMP peri-event traces accroding to these categories:
    %% Was n-1 rewarded?
    % Lever Press Onset
    Grouped_RE_n1_LP_ON = [Grouped_RE_n1_LP_ON; GroupedGCAMP.Mice{session}.baseline_norm_LP_ON(RE_n1_indices,:)];
    Grouped_noRE_n1_LP_ON = [Grouped_noRE_n1_LP_ON; GroupedGCAMP.Mice{session}.baseline_norm_LP_ON(noRE_n1_indices,:)];
    % Lever Press Offset
    Grouped_RE_n1_LP_OFF = [Grouped_RE_n1_LP_OFF; GroupedGCAMP.Mice{session}.baseline_norm_LP_OFF(RE_n1_indices,:)];
    Grouped_noRE_n1_LP_OFF = [Grouped_noRE_n1_LP_OFF; GroupedGCAMP.Mice{session}.baseline_norm_LP_OFF(noRE_n1_indices,:)];
    % Lever Press Duration (interpolated activity)
    Grouped_RE_n1_LPInterp = [Grouped_RE_n1_LPInterp; GroupedGCAMP.Mice{session}.baseline_norm_Duration(RE_n1_indices,:)];
    Grouped_noRE_n1_LPInterp = [Grouped_noRE_n1_LPInterp;GroupedGCAMP.Mice{session}.baseline_norm_Duration(noRE_n1_indices,:)];
    
    %% If n-0 was rewaded, was n-1 rewarded?
    % Lever Press Onset
    Grouped_n0_RE_RE_n1_LP_ON = [Grouped_n0_RE_RE_n1_LP_ON; GroupedGCAMP.Mice{session}.baseline_norm_LP_ON(RE_n0_RE_n1_indices,:)];
    Grouped_n0_RE_noRE_n1_LP_ON = [Grouped_n0_RE_noRE_n1_LP_ON; GroupedGCAMP.Mice{session}.baseline_norm_LP_ON(RE_n0_noRE_n1_indices,:)];
    % Lever Press Offset
    Grouped_n0_RE_RE_n1_LP_OFF = [Grouped_n0_RE_RE_n1_LP_OFF; GroupedGCAMP.Mice{session}.baseline_norm_LP_OFF(RE_n0_RE_n1_indices,:)];
    Grouped_n0_RE_noRE_n1_LP_OFF = [Grouped_n0_RE_noRE_n1_LP_OFF; GroupedGCAMP.Mice{session}.baseline_norm_LP_OFF(RE_n0_noRE_n1_indices,:)];
    % Lever Press Duration (interpolated activity)
    Grouped_n0_RE_RE_n1_LPInterp = [Grouped_n0_RE_RE_n1_LPInterp; GroupedGCAMP.Mice{session}.baseline_norm_Duration(RE_n0_RE_n1_indices,:)];
    Grouped_n0_RE_noRE_n1_LPInterp = [Grouped_n0_RE_noRE_n1_LPInterp;GroupedGCAMP.Mice{session}.baseline_norm_Duration(RE_n0_noRE_n1_indices,:)];
    %% If n-0 was not rewaded, was n-1 rewarded?
    % Lever Press Onset
    Grouped_n0_noRE_RE_n1_LP_ON = [Grouped_n0_noRE_RE_n1_LP_ON; GroupedGCAMP.Mice{session}.baseline_norm_LP_ON(noRE_n0_RE_n1_indices,:)];
    Grouped_n0_noRE_noRE_n1_LP_ON = [Grouped_n0_noRE_noRE_n1_LP_ON; GroupedGCAMP.Mice{session}.baseline_norm_LP_ON(noRE_n0_noRE_n1_indices,:)];
    % Lever Press Offset
    Grouped_n0_noRE_RE_n1_LP_OFF = [Grouped_n0_noRE_RE_n1_LP_OFF; GroupedGCAMP.Mice{session}.baseline_norm_LP_OFF(noRE_n0_RE_n1_indices,:)];
    Grouped_n0_noRE_noRE_n1_LP_OFF = [Grouped_n0_noRE_noRE_n1_LP_OFF; GroupedGCAMP.Mice{session}.baseline_norm_LP_OFF(noRE_n0_noRE_n1_indices,:)];
    % Lever Press Duration (interpolated activity)
    Grouped_n0_noRE_RE_n1_LPInterp = [Grouped_n0_noRE_RE_n1_LPInterp; GroupedGCAMP.Mice{session}.baseline_norm_Duration(noRE_n0_RE_n1_indices,:)];
    Grouped_n0_noRE_noRE_n1_LPInterp = [Grouped_n0_noRE_noRE_n1_LPInterp;GroupedGCAMP.Mice{session}.baseline_norm_Duration(noRE_n0_noRE_n1_indices,:)];
    
end

%% Insert GCAMP traces into Grouped data structure
%% Was n-1 rewarded?
GroupedGCAMP.RE_n1_LP_ON = Grouped_RE_n1_LP_ON;
GroupedGCAMP.noRE_n1_LP_ON = Grouped_noRE_n1_LP_ON;
GroupedGCAMP.RE_n1_LP_OFF = Grouped_RE_n1_LP_OFF;
GroupedGCAMP.noRE_n1_LP_OFF = Grouped_noRE_n1_LP_OFF;
GroupedGCAMP.RE_n1_LPInterp = Grouped_RE_n1_LPInterp;
GroupedGCAMP.noRE_n1_LPInterp = Grouped_noRE_n1_LPInterp;
%% If n-0 was rewaded, was n-1 rewarded?
GroupedGCAMP.n0_RE_RE_n1_LP_ON = Grouped_n0_RE_RE_n1_LP_ON;
GroupedGCAMP.n0_RE_noRE_n1_LP_ON = Grouped_n0_RE_noRE_n1_LP_ON;
GroupedGCAMP.n0_RE_RE_n1_LP_OFF = Grouped_n0_RE_RE_n1_LP_OFF;
GroupedGCAMP.n0_RE_noRE_n1_LP_OFF = Grouped_n0_RE_noRE_n1_LP_OFF;
GroupedGCAMP.n0_RE_RE_n1_LPInterp = Grouped_n0_RE_RE_n1_LPInterp;
GroupedGCAMP.n0_RE_noRE_n1_LPInterp = Grouped_n0_RE_noRE_n1_LPInterp;
%% If n-0 was not rewaded, was n-1 rewarded?
GroupedGCAMP.n0_noRE_RE_n1_LP_ON = Grouped_n0_noRE_RE_n1_LP_ON;
GroupedGCAMP.n0_noRE_noRE_n1_LP_ON = Grouped_n0_noRE_noRE_n1_LP_ON;
GroupedGCAMP.n0_noRE_RE_n1_LP_OFF = Grouped_n0_noRE_RE_n1_LP_OFF;
GroupedGCAMP.n0_noRE_noRE_n1_LP_OFF = Grouped_n0_noRE_noRE_n1_LP_OFF;
GroupedGCAMP.n0_noRE_RE_n1_LPInterp = Grouped_n0_noRE_RE_n1_LPInterp;
GroupedGCAMP.n0_noRE_noRE_n1_LPInterp = Grouped_n0_noRE_noRE_n1_LPInterp;

end


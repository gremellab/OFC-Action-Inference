%% GCAMP Analysis Master Script
% Christian Cazares, Drew Schreiner, Christina Gremel
% University of California, San Diego
% Last updated by Christian Cazares on 5/02/22
%% Loop through each Excel File in Selected Folder to extract raw behavioral events and fluorescent signals
% Saves each session as a .mat containing that data
GCAMP_Save_Dir = '';
extractData(GCAMP_Save_Dir)
%% Loops through each .mat GCAMP session and peforms peri-event GCaMP trace analysis
% Function:
        % extract_PeriEventAnalysis(base_time_start, base_time_end, time_end)
    % Parameters:
    %   base_time_start: time (in seconds) before onset of lever press for
    %       start of baseline window
    %   base_time_end: time (in seconds) before onset of lever press for
    %       end of baseline window
    %   time_end: time (in seconds) after onset of event for
    %       end of analys window
    
    base_time_start = -10;
    base_time_end = -2;
    time_end = 5;
    extract_PeriEventAnalysis(base_time_start, base_time_end, time_end)
    
% 1. Extracts timestamps of each behavioral event (Onset, Offset, Head Entry, Reward)
    % Subfunction: 
            % extract_EventTimestamps(GCAMP)
        % Parameters:
            % GCAMP: Session including extract raw excel timestamps of
            % events and photometry signals
% 2. Baseline normalization of each session's peri-event GCAMP traces (Onset, Offset, Head Entry, Reward, Quantile, Interpolated)
    % Subfunction: 
            % baselineNormalized_PeriEventTraces(GCAMP)
        % Parameters:
            % GCAMP: Session including extract raw excel timestamps of
            % events and photometry signals
% 3. Performs Linear Mixed Effect Model (LME) Analysis on each session's peri-event GCAMP
%       traces (pre-Onset, post-Onset, post-Offset)

%% Plots example session peri-event GCAMP traces (onset, offset, first head entry after reward, duration)
plotExampleSession(GCAMP)
%% Group GCAMP data
[Grouped_GCAMP] = groupGCAMPData();
%% Regression
[Model] = GCAMP_grand_regression_indivshuffles_r(GroupedGCAMP_OFC);
save('Model','Model' , '-v7.3')
%% N-1
% If N-1 is Rewarded or not
GroupedGCAMP = nBackReward(GroupedGCAMP);
% What duration distribution quartile does N-1 belong to?
GroupedGCAMP = nBackQuartile(GroupedGCAMP);
%% Perform Statistical Comparisons in Grouped GCAMP data
[GroupedGCAMP] = groupGCAMPData_Stats(GroupedGCAMP);
%% Plots grouped session peri-event GCAMP traces (onset, offset, first head entry after reward, duration)
% Loops through each individual GCAMP session in folder and places
% peri-event data traces into a cohort-sized data structure
plotGroupedSessions(GroupedGCAMP);
function [GroupedGCAMP] = groupGCAMPData_Stats(GroupedGCAMP)


%% LP Onset Met vs Fail
[p_val, observeddifference] = permTest_array(GroupedGCAMP.LPON_Met, GroupedGCAMP.LPON_Fail, 1000);
threshold = .05;
sample_span = 5;
aboveThreshold = p_val <= threshold;

%aboveThreshold is a logical array, where 1 when above threshold, 0, below.
%we thus want to calculate the difference between rising and falling edges
aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
edges = diff(aboveThreshold);
rising = find(edges==1);     %rising/falling edges
falling = find(edges==-1);
spanWidth = falling - rising;  %width of span of 1's (above threshold)
wideEnough = spanWidth >= sample_span;
startPos = rising(wideEnough);    %start of each span
endPos = falling(wideEnough)-1;   %end of each span
%all points which are in the sample span (i.e. between startPos and endPos).
allInSpan = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));

%now use that to index into the p_val thing and turn the significant
%consecutive ones into 0s, and all others into Nans so that we can plot
%that into the original graph to show where sig. diff exist.
p_val_for_graph = p_val;
p_val_for_graph(allInSpan) =0;
non_sig_idx =p_val_for_graph > 0;
p_val_for_graph(non_sig_idx) = NaN;

GroupedGCAMP.Stats.LPON.p_val = p_val;
GroupedGCAMP.Stats.LPON.p_val_for_graph = p_val_for_graph;

%% LP Offset Met vs Fail
[p_val, observeddifference] = permTest_array(GroupedGCAMP.LPOFF_Met, GroupedGCAMP.LPOFF_Fail, 1000);
threshold = .05;
sample_span = 5;
aboveThreshold = p_val <= threshold;

%aboveThreshold is a logical array, where 1 when above threshold, 0, below.
%we thus want to calculate the difference between rising and falling edges
aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
edges = diff(aboveThreshold);
rising = find(edges==1);     %rising/falling edges
falling = find(edges==-1);
spanWidth = falling - rising;  %width of span of 1's (above threshold)
wideEnough = spanWidth >= sample_span;
startPos = rising(wideEnough);    %start of each span
endPos = falling(wideEnough)-1;   %end of each span
%all points which are in the sample span (i.e. between startPos and endPos).
allInSpan = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));

%now use that to index into the p_val thing and turn the significant
%consecutive ones into 0s, and all others into Nans so that we can plot
%that into the original graph to show where sig. diff exist.
p_val_for_graph = p_val;
p_val_for_graph(allInSpan) =0;
non_sig_idx =p_val_for_graph > 0;
p_val_for_graph(non_sig_idx) = NaN;

GroupedGCAMP.Stats.LPOFF.p_val = p_val;
GroupedGCAMP.Stats.LPOFF.p_val_for_graph = p_val_for_graph;

%% First Head Entry after Met vs Fail
%% Head Entry (Success Rewarded vs Failure)
[p_val, observeddifference] = permTest_array(GroupedGCAMP.FirstHE_RE, GroupedGCAMP.FirstHE_Fail, 1000);
threshold = .05;
sample_span = 5;
aboveThreshold = p_val <= threshold;

%aboveThreshold is a logical array, where 1 when above threshold, 0, below.
%we thus want to calculate the difference between rising and falling edges
aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
edges = diff(aboveThreshold);
rising = find(edges==1);     %rising/falling edges
falling = find(edges==-1);
spanWidth = falling - rising;  %width of span of 1's (above threshold)
wideEnough = spanWidth >= sample_span;
startPos = rising(wideEnough);    %start of each span
endPos = falling(wideEnough)-1;   %end of each span
%all points which are in the sample span (i.e. between startPos and endPos).
allInSpan = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));

%now use that to index into the p_val thing and turn the significant
%consecutive ones into 0s, and all others into Nans so that we can plot
%that into the original graph to show where sig. diff exist.
p_val_for_graph = p_val;
p_val_for_graph(allInSpan) =0;
non_sig_idx =p_val_for_graph > 0;
p_val_for_graph(non_sig_idx) = NaN;

GroupedGCAMP.Stats.HE_Fail.p_val = p_val;
GroupedGCAMP.Stats.HE_Fail.p_val_for_graph = p_val_for_graph;

%% LP Duration Interpolated Met vs Fail
% prior to comparison, remove any NaNs (which are LPs withouth at least 2
% samples)

Grouped_LPInterp_Met = rmmissing(GroupedGCAMP.LPInterp_Met);
Grouped_LPInterp_Fail = rmmissing(GroupedGCAMP.LPInterp_Fail);

[p_val, observeddifference] = permTest_array(Grouped_LPInterp_Met, Grouped_LPInterp_Fail, 1000); %permutations [, varargin]
threshold = .05;
sample_span = 3;
aboveThreshold = p_val <= threshold;

%aboveThreshold is a logical array, where 1 when above threshold, 0, below.
%we thus want to calculate the difference between rising and falling edges
aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
edges = diff(aboveThreshold);
rising = find(edges==1);     %rising/falling edges
falling = find(edges==-1);
spanWidth = falling - rising;  %width of span of 1's (above threshold)
wideEnough = spanWidth >= sample_span;
startPos = rising(wideEnough);    %start of each span
endPos = falling(wideEnough)-1;   %end of each span
%all points which are in the sample span (i.e. between startPos and endPos).
allInSpan = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));

%now use that to index into the p_val thing and turn the significant
%consecutive ones into 0s, and all others into Nans so that we can plot
%that into the original graph to show where sig. diff exist.
p_val_for_graph = p_val;
p_val_for_graph(allInSpan) =0;
non_sig_idx =p_val_for_graph > 0;
p_val_for_graph(non_sig_idx) = NaN;

GroupedGCAMP.Stats.LPInterp.p_val = p_val;
GroupedGCAMP.Stats.LPInterp.p_val_for_graph = p_val_for_graph;


%% LP Onset n-1 Rewarded vs not Rewarded
[p_val, observeddifference] = permTest_array(GroupedGCAMP.RE_n1_LP_ON, GroupedGCAMP.noRE_n1_LP_ON, 1000);
threshold = .05;
sample_span = 5;
aboveThreshold = p_val <= threshold;

%aboveThreshold is a logical array, where 1 when above threshold, 0, below.
%we thus want to calculate the difference between rising and falling edges
aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
edges = diff(aboveThreshold);
rising = find(edges==1);     %rising/falling edges
falling = find(edges==-1);
spanWidth = falling - rising;  %width of span of 1's (above threshold)
wideEnough = spanWidth >= sample_span;
startPos = rising(wideEnough);    %start of each span
endPos = falling(wideEnough)-1;   %end of each span
%all points which are in the sample span (i.e. between startPos and endPos).
allInSpan = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));

%now use that to index into the p_val thing and turn the significant
%consecutive ones into 0s, and all others into Nans so that we can plot
%that into the original graph to show where sig. diff exist.
p_val_for_graph = p_val;
p_val_for_graph(allInSpan) =0;
non_sig_idx =p_val_for_graph > 0;
p_val_for_graph(non_sig_idx) = NaN;

GroupedGCAMP.Stats.LPON_n1_reward.p_val = p_val;
GroupedGCAMP.Stats.LPON_n1_reward.p_val_for_graph = p_val_for_graph;

%% LP Offset n-1 Rewarded vs not Rewarded
[p_val, observeddifference] = permTest_array(GroupedGCAMP.RE_n1_LP_OFF, GroupedGCAMP.noRE_n1_LP_OFF, 1000);
threshold = .05;
sample_span = 5;
aboveThreshold = p_val <= threshold;

%aboveThreshold is a logical array, where 1 when above threshold, 0, below.
%we thus want to calculate the difference between rising and falling edges
aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
edges = diff(aboveThreshold);
rising = find(edges==1);     %rising/falling edges
falling = find(edges==-1);
spanWidth = falling - rising;  %width of span of 1's (above threshold)
wideEnough = spanWidth >= sample_span;
startPos = rising(wideEnough);    %start of each span
endPos = falling(wideEnough)-1;   %end of each span
%all points which are in the sample span (i.e. between startPos and endPos).
allInSpan = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));

%now use that to index into the p_val thing and turn the significant
%consecutive ones into 0s, and all others into Nans so that we can plot
%that into the original graph to show where sig. diff exist.
p_val_for_graph = p_val;
p_val_for_graph(allInSpan) =0;
non_sig_idx =p_val_for_graph > 0;
p_val_for_graph(non_sig_idx) = NaN;

GroupedGCAMP.Stats.LPOFF_n1_reward.p_val = p_val;
GroupedGCAMP.Stats.LPOFF_n1_reward.p_val_for_graph = p_val_for_graph;

%% LP Duration Interpolated n-1 Rewarded vs not Rewarded
% prior to comparison, remove any NaNs (which are LPs withouth at least 2
% samples)

Grouped_RE_n1_LPInterp = rmmissing(GroupedGCAMP.RE_n1_LPInterp);
Grouped_noRE_n1_LPInterp = rmmissing(GroupedGCAMP.noRE_n1_LPInterp);

[p_val, observeddifference] = permTest_array(Grouped_RE_n1_LPInterp, Grouped_noRE_n1_LPInterp, 1000); %permutations [, varargin]
threshold = .05;
sample_span = 3;
aboveThreshold = p_val <= threshold;

%aboveThreshold is a logical array, where 1 when above threshold, 0, below.
%we thus want to calculate the difference between rising and falling edges
aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
edges = diff(aboveThreshold);
rising = find(edges==1);     %rising/falling edges
falling = find(edges==-1);
spanWidth = falling - rising;  %width of span of 1's (above threshold)
wideEnough = spanWidth >= sample_span;
startPos = rising(wideEnough);    %start of each span
endPos = falling(wideEnough)-1;   %end of each span
%all points which are in the sample span (i.e. between startPos and endPos).
allInSpan = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));

%now use that to index into the p_val thing and turn the significant
%consecutive ones into 0s, and all others into Nans so that we can plot
%that into the original graph to show where sig. diff exist.
p_val_for_graph = p_val;
p_val_for_graph(allInSpan) =0;
non_sig_idx =p_val_for_graph > 0;
p_val_for_graph(non_sig_idx) = NaN;

GroupedGCAMP.Stats.LPInterp_n1_reward.p_val = p_val;
GroupedGCAMP.Stats.LPInterp_n1_reward.p_val_for_graph = p_val_for_graph;



%% LP Onset: If n-0 was rewarded, n-1 Rewarded vs not Rewarded
[p_val, observeddifference] = permTest_array(GroupedGCAMP.n0_RE_RE_n1_LP_ON, GroupedGCAMP.n0_RE_noRE_n1_LP_ON, 1000);
threshold = .05;
sample_span = 5;
aboveThreshold = p_val <= threshold;

%aboveThreshold is a logical array, where 1 when above threshold, 0, below.
%we thus want to calculate the difference between rising and falling edges
aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
edges = diff(aboveThreshold);
rising = find(edges==1);     %rising/falling edges
falling = find(edges==-1);
spanWidth = falling - rising;  %width of span of 1's (above threshold)
wideEnough = spanWidth >= sample_span;
startPos = rising(wideEnough);    %start of each span
endPos = falling(wideEnough)-1;   %end of each span
%all points which are in the sample span (i.e. between startPos and endPos).
allInSpan = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));

%now use that to index into the p_val thing and turn the significant
%consecutive ones into 0s, and all others into Nans so that we can plot
%that into the original graph to show where sig. diff exist.
p_val_for_graph = p_val;
p_val_for_graph(allInSpan) =0;
non_sig_idx =p_val_for_graph > 0;
p_val_for_graph(non_sig_idx) = NaN;

GroupedGCAMP.Stats.LPON_n0_reward_n1_reward.p_val = p_val;
GroupedGCAMP.Stats.LPON_n0_reward_n1_reward.p_val_for_graph = p_val_for_graph;

%% LP Offset: If n-0 was rewarded, n-1 Rewarded vs not Rewarded
[p_val, observeddifference] = permTest_array(GroupedGCAMP.n0_RE_RE_n1_LP_OFF, GroupedGCAMP.n0_RE_noRE_n1_LP_OFF, 1000);
threshold = .05;
sample_span = 5;
aboveThreshold = p_val <= threshold;

%aboveThreshold is a logical array, where 1 when above threshold, 0, below.
%we thus want to calculate the difference between rising and falling edges
aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
edges = diff(aboveThreshold);
rising = find(edges==1);     %rising/falling edges
falling = find(edges==-1);
spanWidth = falling - rising;  %width of span of 1's (above threshold)
wideEnough = spanWidth >= sample_span;
startPos = rising(wideEnough);    %start of each span
endPos = falling(wideEnough)-1;   %end of each span
%all points which are in the sample span (i.e. between startPos and endPos).
allInSpan = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));

%now use that to index into the p_val thing and turn the significant
%consecutive ones into 0s, and all others into Nans so that we can plot
%that into the original graph to show where sig. diff exist.
p_val_for_graph = p_val;
p_val_for_graph(allInSpan) =0;
non_sig_idx =p_val_for_graph > 0;
p_val_for_graph(non_sig_idx) = NaN;

GroupedGCAMP.Stats.LPOFF_n0_reward_n1_reward.p_val = p_val;
GroupedGCAMP.Stats.LPOFF_n0_reward_n1_reward.p_val_for_graph = p_val_for_graph;

%% LP Duration Interpolated n-1 Rewarded vs not Rewarded
% prior to comparison, remove any NaNs (which are LPs withouth at least 2
% samples)

Grouped_n0_RE_RE_n1_LPInterp = rmmissing(GroupedGCAMP.n0_RE_RE_n1_LPInterp);
Grouped_n0_RE_noRE_n1_LPInterp = rmmissing(GroupedGCAMP.n0_RE_noRE_n1_LPInterp);

[p_val, observeddifference] = permTest_array(Grouped_n0_RE_RE_n1_LPInterp, Grouped_n0_RE_noRE_n1_LPInterp, 1000); %permutations [, varargin]
threshold = .05;
sample_span = 3;
aboveThreshold = p_val <= threshold;

%aboveThreshold is a logical array, where 1 when above threshold, 0, below.
%we thus want to calculate the difference between rising and falling edges
aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
edges = diff(aboveThreshold);
rising = find(edges==1);     %rising/falling edges
falling = find(edges==-1);
spanWidth = falling - rising;  %width of span of 1's (above threshold)
wideEnough = spanWidth >= sample_span;
startPos = rising(wideEnough);    %start of each span
endPos = falling(wideEnough)-1;   %end of each span
%all points which are in the sample span (i.e. between startPos and endPos).
allInSpan = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));

%now use that to index into the p_val thing and turn the significant
%consecutive ones into 0s, and all others into Nans so that we can plot
%that into the original graph to show where sig. diff exist.
p_val_for_graph = p_val;
p_val_for_graph(allInSpan) =0;
non_sig_idx =p_val_for_graph > 0;
p_val_for_graph(non_sig_idx) = NaN;

GroupedGCAMP.Stats.LPInterp_n0_reward_n1_reward.p_val = p_val;
GroupedGCAMP.Stats.LPInterp_n0_reward_n1_reward.p_val_for_graph = p_val_for_graph;





%% LP Onset: If n-0 was not rewarded, n-1 Rewarded vs not Rewarded
[p_val, observeddifference] = permTest_array(GroupedGCAMP.n0_noRE_RE_n1_LP_ON, GroupedGCAMP.n0_noRE_noRE_n1_LP_ON, 1000);
threshold = .05;
sample_span = 5;
aboveThreshold = p_val <= threshold;

%aboveThreshold is a logical array, where 1 when above threshold, 0, below.
%we thus want to calculate the difference between rising and falling edges
aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
edges = diff(aboveThreshold);
rising = find(edges==1);     %rising/falling edges
falling = find(edges==-1);
spanWidth = falling - rising;  %width of span of 1's (above threshold)
wideEnough = spanWidth >= sample_span;
startPos = rising(wideEnough);    %start of each span
endPos = falling(wideEnough)-1;   %end of each span
%all points which are in the sample span (i.e. between startPos and endPos).
allInSpan = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));

%now use that to index into the p_val thing and turn the significant
%consecutive ones into 0s, and all others into Nans so that we can plot
%that into the original graph to show where sig. diff exist.
p_val_for_graph = p_val;
p_val_for_graph(allInSpan) =0;
non_sig_idx =p_val_for_graph > 0;
p_val_for_graph(non_sig_idx) = NaN;

GroupedGCAMP.Stats.LPON_n0_no_reward_n1_reward.p_val = p_val;
GroupedGCAMP.Stats.LPON_n0_no_reward_n1_reward.p_val_for_graph = p_val_for_graph;

%% LP Offset: If n-0 was not rewarded, n-1 Rewarded vs not Rewarded
[p_val, observeddifference] = permTest_array(GroupedGCAMP.n0_noRE_RE_n1_LP_OFF, GroupedGCAMP.n0_noRE_noRE_n1_LP_OFF, 1000);
threshold = .05;
sample_span = 5;
aboveThreshold = p_val <= threshold;

%aboveThreshold is a logical array, where 1 when above threshold, 0, below.
%we thus want to calculate the difference between rising and falling edges
aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
edges = diff(aboveThreshold);
rising = find(edges==1);     %rising/falling edges
falling = find(edges==-1);
spanWidth = falling - rising;  %width of span of 1's (above threshold)
wideEnough = spanWidth >= sample_span;
startPos = rising(wideEnough);    %start of each span
endPos = falling(wideEnough)-1;   %end of each span
%all points which are in the sample span (i.e. between startPos and endPos).
allInSpan = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));

%now use that to index into the p_val thing and turn the significant
%consecutive ones into 0s, and all others into Nans so that we can plot
%that into the original graph to show where sig. diff exist.
p_val_for_graph = p_val;
p_val_for_graph(allInSpan) =0;
non_sig_idx =p_val_for_graph > 0;
p_val_for_graph(non_sig_idx) = NaN;

GroupedGCAMP.Stats.LPOFF_n0_no_reward_n1_reward.p_val = p_val;
GroupedGCAMP.Stats.LPOFF_n0_no_reward_n1_reward.p_val_for_graph = p_val_for_graph;

%% LP Duration: If n-0 was not rewarded, n-1 Rewarded vs not Rewarded
% prior to comparison, remove any NaNs (which are LPs withouth at least 2
% samples)

Grouped_n0_noRE_RE_n1_LPInterp = rmmissing(GroupedGCAMP.n0_noRE_RE_n1_LPInterp);
Grouped_n0_noRE_noRE_n1_LPInterp = rmmissing(GroupedGCAMP.n0_noRE_noRE_n1_LPInterp);

[p_val, observeddifference] = permTest_array(Grouped_n0_noRE_RE_n1_LPInterp, Grouped_n0_noRE_noRE_n1_LPInterp, 1000); %permutations [, varargin]
threshold = .05;
sample_span = 3;
aboveThreshold = p_val <= threshold;

%aboveThreshold is a logical array, where 1 when above threshold, 0, below.
%we thus want to calculate the difference between rising and falling edges
aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
edges = diff(aboveThreshold);
rising = find(edges==1);     %rising/falling edges
falling = find(edges==-1);
spanWidth = falling - rising;  %width of span of 1's (above threshold)
wideEnough = spanWidth >= sample_span;
startPos = rising(wideEnough);    %start of each span
endPos = falling(wideEnough)-1;   %end of each span
%all points which are in the sample span (i.e. between startPos and endPos).
allInSpan = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));

%now use that to index into the p_val thing and turn the significant
%consecutive ones into 0s, and all others into Nans so that we can plot
%that into the original graph to show where sig. diff exist.
p_val_for_graph = p_val;
p_val_for_graph(allInSpan) =0;
non_sig_idx =p_val_for_graph > 0;
p_val_for_graph(non_sig_idx) = NaN;

GroupedGCAMP.Stats.LPInterp_n0_no_reward_n1_reward.p_val = p_val;
GroupedGCAMP.Stats.LPInterp_n0_no_reward_n1_reward.p_val_for_graph = p_val_for_graph;


end


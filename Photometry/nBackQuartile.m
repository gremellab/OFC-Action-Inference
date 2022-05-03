function [GroupedGCAMP] = nBackQuartile(GroupedGCAMP)

total_sessions = length(GroupedGCAMP.Mice);


Grouped_Q1_n_LP_ON = [];
Grouped_Q2_n_LP_ON = [];
Grouped_Q3_n_LP_ON = [];
Grouped_Q4_n_LP_ON = [];

Grouped_Q1_n_LP_OFF = [];
Grouped_Q2_n_LP_OFF = [];
Grouped_Q3_n_LP_OFF = [];
Grouped_Q4_n_LP_OFF = [];

Grouped_Q1_n_LPInterp = [];
Grouped_Q2_n_LPInterp = [];
Grouped_Q3_n_LPInterp = [];
Grouped_Q4_n_LPInterp = [];


Grouped_Q1_n1_LP_ON = [];
Grouped_Q2_n1_LP_ON = [];
Grouped_Q3_n1_LP_ON = [];
Grouped_Q4_n1_LP_ON = [];

Grouped_Q1_n1_LP_OFF = [];
Grouped_Q2_n1_LP_OFF = [];
Grouped_Q3_n1_LP_OFF = [];
Grouped_Q4_n1_LP_OFF = [];

Grouped_Q1_n1_LPInterp = [];
Grouped_Q2_n1_LPInterp = [];
Grouped_Q3_n1_LPInterp = [];
Grouped_Q4_n1_LPInterp = [];

quartiles_n = [];
quartiles_n1 = [];
%% Quartiles
histogram_bin_indexs_all =[];
for session = 1:total_sessions
    Lengths =  cell2mat(GroupedGCAMP.Mice{session}.regression_cell(:,1));
    edge_start = [0 0.25 0.50 0.75];
    edge_end = [0.25 0.50 0.75 1];
    Counts = [];
    True_edges = zeros(4,2);
    for quart = 1:4
        edges = quantile(Lengths,[edge_start(quart) edge_end(quart)]);
        if quart == 1
            quart_lenghts = Lengths(Lengths >= edges(1) & Lengths <= edges(2));
        else
            quart_lenghts = Lengths(Lengths > edges(1) & Lengths <= edges(2));
        end
        True_edges(quart,:) = edges;
        [N,~,bin]= histcounts(quart_lenghts, edges);
        Counts = [Counts N];
    end
    
    hist_edges =[ True_edges(:,1);True_edges(4,end)];
    [histogram_counts,histogram_edges,histogram_bin_indexs] = histcounts( Lengths,[hist_edges]);
%     histogram_bin_indexs_all = [histogram_bin_indexs_all; histogram_bin_indexs];
%     quartiles_n = [quartiles_n; histogram_bin_indexs_all];
%     quartiles_n1_pad = [NaN ; quartiles_n(1:end-1)];
%     quartiles_n1 = [quartiles_n1; quartiles_n1_pad];
    
    quartiles_n = histogram_bin_indexs;
    quartiles_n1 = [NaN ; histogram_bin_indexs(1:end-1)];
    
    GroupedGCAMP.Mice{session}.regression_cell(:,end+1) = num2cell(histogram_bin_indexs, length(histogram_bin_indexs));
    GroupedGCAMP.Mice{session}.regression_cell(:,end+1) = num2cell(quartiles_n1, length(quartiles_n1));
    GroupedGCAMP.Mice{session}.data_table_variables(:,end+1) = {'n_quartile'};
    GroupedGCAMP.Mice{session}.data_table_variables(:,end+1) = {'n1_quartile'};
    
    n_Q1_idx = find(quartiles_n == 1);
    n_Q2_idx = find(quartiles_n == 2);
    n_Q3_idx = find(quartiles_n == 3);
    n_Q4_idx = find(quartiles_n == 4);
    
    n1_Q1_idx = find(quartiles_n1 == 1);
    n1_Q2_idx = find(quartiles_n1 == 2);
    n1_Q3_idx = find(quartiles_n1 == 3);
    n1_Q4_idx = find(quartiles_n1 == 4);
    
    % LP ON
    Grouped_Q1_n_LP_ON = [Grouped_Q1_n_LP_ON; GroupedGCAMP.Mice{session}.baseline_norm_LP_ON(n_Q1_idx,:)];
    Grouped_Q2_n_LP_ON = [Grouped_Q2_n_LP_ON; GroupedGCAMP.Mice{session}.baseline_norm_LP_ON(n_Q2_idx,:)];
    Grouped_Q3_n_LP_ON = [Grouped_Q3_n_LP_ON; GroupedGCAMP.Mice{session}.baseline_norm_LP_ON(n_Q3_idx,:)];
    Grouped_Q4_n_LP_ON = [Grouped_Q4_n_LP_ON; GroupedGCAMP.Mice{session}.baseline_norm_LP_ON(n_Q4_idx,:)];
    
    Grouped_Q1_n1_LP_ON = [Grouped_Q1_n1_LP_ON; GroupedGCAMP.Mice{session}.baseline_norm_LP_ON(n1_Q1_idx,:)];
    Grouped_Q2_n1_LP_ON = [Grouped_Q2_n1_LP_ON; GroupedGCAMP.Mice{session}.baseline_norm_LP_ON(n1_Q2_idx,:)];
    Grouped_Q3_n1_LP_ON = [Grouped_Q3_n1_LP_ON; GroupedGCAMP.Mice{session}.baseline_norm_LP_ON(n1_Q3_idx,:)];
    Grouped_Q4_n1_LP_ON = [Grouped_Q4_n1_LP_ON; GroupedGCAMP.Mice{session}.baseline_norm_LP_ON(n1_Q4_idx,:)];
    
    % LP OFF
    Grouped_Q1_n_LP_OFF = [Grouped_Q1_n_LP_OFF; GroupedGCAMP.Mice{session}.baseline_norm_LP_OFF(n_Q1_idx,:)];
    Grouped_Q2_n_LP_OFF = [Grouped_Q2_n_LP_OFF; GroupedGCAMP.Mice{session}.baseline_norm_LP_OFF(n_Q2_idx,:)];
    Grouped_Q3_n_LP_OFF = [Grouped_Q3_n_LP_OFF; GroupedGCAMP.Mice{session}.baseline_norm_LP_OFF(n_Q3_idx,:)];
    Grouped_Q4_n_LP_OFF = [Grouped_Q4_n_LP_OFF; GroupedGCAMP.Mice{session}.baseline_norm_LP_OFF(n_Q4_idx,:)];
    
    Grouped_Q1_n1_LP_OFF = [Grouped_Q1_n1_LP_OFF; GroupedGCAMP.Mice{session}.baseline_norm_LP_OFF(n1_Q1_idx,:)];
    Grouped_Q2_n1_LP_OFF = [Grouped_Q2_n1_LP_OFF; GroupedGCAMP.Mice{session}.baseline_norm_LP_OFF(n1_Q2_idx,:)];
    Grouped_Q3_n1_LP_OFF = [Grouped_Q3_n1_LP_OFF; GroupedGCAMP.Mice{session}.baseline_norm_LP_OFF(n1_Q3_idx,:)];
    Grouped_Q4_n1_LP_OFF = [Grouped_Q4_n1_LP_OFF; GroupedGCAMP.Mice{session}.baseline_norm_LP_OFF(n1_Q4_idx,:)];
    
    % Interpolated Duration
    Grouped_Q1_n_LPInterp = [Grouped_Q1_n_LPInterp; GroupedGCAMP.Mice{session}.baseline_norm_Duration(n_Q1_idx,:)];
    Grouped_Q2_n_LPInterp = [Grouped_Q2_n_LPInterp; GroupedGCAMP.Mice{session}.baseline_norm_Duration(n_Q2_idx,:)];
    Grouped_Q3_n_LPInterp = [Grouped_Q3_n_LPInterp; GroupedGCAMP.Mice{session}.baseline_norm_Duration(n_Q3_idx,:)];
    Grouped_Q4_n_LPInterp = [Grouped_Q4_n_LPInterp; GroupedGCAMP.Mice{session}.baseline_norm_Duration(n_Q4_idx,:)];
    
    Grouped_Q1_n1_LPInterp = [Grouped_Q1_n1_LPInterp; GroupedGCAMP.Mice{session}.baseline_norm_Duration(n1_Q1_idx,:)];
    Grouped_Q2_n1_LPInterp = [Grouped_Q2_n1_LPInterp; GroupedGCAMP.Mice{session}.baseline_norm_Duration(n1_Q2_idx,:)];
    Grouped_Q3_n1_LPInterp = [Grouped_Q3_n1_LPInterp; GroupedGCAMP.Mice{session}.baseline_norm_Duration(n1_Q3_idx,:)];
    Grouped_Q4_n1_LPInterp = [Grouped_Q4_n1_LPInterp; GroupedGCAMP.Mice{session}.baseline_norm_Duration(n1_Q4_idx,:)];

end


    % LP ON
    GroupedGCAMP.Q1_n_LP_ON = Grouped_Q1_n_LP_ON;
    GroupedGCAMP.Q2_n_LP_ON = Grouped_Q2_n_LP_ON;
    GroupedGCAMP.Q3_n_LP_ON = Grouped_Q3_n_LP_ON;
    GroupedGCAMP.Q4_n_LP_ON = Grouped_Q4_n_LP_ON;
    
    GroupedGCAMP.Q1_n1_LP_ON = Grouped_Q1_n1_LP_ON;
    GroupedGCAMP.Q2_n1_LP_ON = Grouped_Q2_n1_LP_ON;
    GroupedGCAMP.Q3_n1_LP_ON = Grouped_Q3_n1_LP_ON;
    GroupedGCAMP.Q4_n1_LP_ON = Grouped_Q4_n1_LP_ON;
    
    % LP OFF
    GroupedGCAMP.Q1_n_LP_OFF = Grouped_Q1_n_LP_OFF;
    GroupedGCAMP.Q2_n_LP_OFF = Grouped_Q2_n_LP_OFF;
    GroupedGCAMP.Q3_n_LP_OFF = Grouped_Q3_n_LP_OFF;
    GroupedGCAMP.Q4_n_LP_OFF = Grouped_Q4_n_LP_OFF;
    
    GroupedGCAMP.Q1_n1_LP_OFF = Grouped_Q1_n1_LP_OFF;
    GroupedGCAMP.Q2_n1_LP_OFF = Grouped_Q2_n1_LP_OFF;
    GroupedGCAMP.Q3_n1_LP_OFF = Grouped_Q3_n1_LP_OFF;
    GroupedGCAMP.Q4_n1_LP_OFF = Grouped_Q4_n1_LP_OFF;
    
    % Interpolated Duration
    GroupedGCAMP.Q1_n_LPInterp = Grouped_Q1_n_LPInterp;
    GroupedGCAMP.Q2_n_LPInterp = Grouped_Q2_n_LPInterp;
    GroupedGCAMP.Q3_n_LPInterp = Grouped_Q3_n_LPInterp;
    GroupedGCAMP.Q4_n_LPInterp = Grouped_Q4_n_LPInterp;
    
    GroupedGCAMP.Q1_n1_LPInterp = Grouped_Q1_n1_LPInterp;
    GroupedGCAMP.Q2_n1_LPInterp = Grouped_Q2_n1_LPInterp;
    GroupedGCAMP.Q3_n1_LPInterp = Grouped_Q3_n1_LPInterp;
    GroupedGCAMP.Q4_n1_LPInterp = Grouped_Q4_n1_LPInterp;

end



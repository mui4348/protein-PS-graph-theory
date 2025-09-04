clear

for ii = 270:10:340

temp = int2str(ii);

fileName1 = strcat('wt_1C_protNoPBC_1k_',temp,'K.pdb');
fileName2 = strcat('wt_2C_protNoPBC_1k_',temp,'K.pdb');
fileName3 = strcat('wt_3C_protNoPBC_1k_',temp,'K.pdb');

threshold_1 = 8;
chainNum = 60;
boxsize = 20;

[ChainPair1,sig_node_count1] = ChainPairDFF(fileName1,threshold_1);
[ChainPair2,sig_node_count2] = ChainPairDFF(fileName2,threshold_1);
[ChainPair3,sig_node_count3] = ChainPairDFF(fileName3,threshold_1);
%%
ChainPair = [ChainPair1,ChainPair2,ChainPair3];
%%
Chain_Adjacency_matrices = {};
Chain_Adjacency_matrices_Weighted = cell(1,numel(ChainPair));
Chain_degree = {};
Efficiency = {};
for i=1:numel(ChainPair)
    logical_indices{i} = ChainPair{i} < threshold_1 & ChainPair{i} ~= 0;
    Chain_Adjacency_matrices{i} = double(logical_indices{i});
    Chain_Adjacency_matrices_Weighted{i} = zeros(chainNum,chainNum);
    for j=1:chainNum
        for l=1:chainNum
            if ChainPair{i}(j,l) < threshold_1
                Chain_Adjacency_matrices_Weighted{i}(j,l) = 1/ChainPair{i}(j,l);
                Chain_Adjacency_matrices_distance{i}(j,l) = ChainPair{i}(j,l);
                Chain_Adjacency_matrices_Weighted{i}(j,j) = 0;
            end
        end
    end   

% Node degree
for j=1:chainNum
    Chain_degree{i}(j) = sum(Chain_Adjacency_matrices{i}(j,:));
end
Graph_average_degree(i) = mean(Chain_degree{i});
% Node weight
for j=1:chainNum
    Chain_weight{i}(j) = sum(Chain_Adjacency_matrices_Weighted{i}(j,:));
end

% Degree Centrality
for j=1:chainNum
    Chain_centrality{i}(j) = Chain_degree{i}(j)/chainNum-1;
end

% Weight Centrality
for j=1:chainNum
    Weight_centrality{i}(j) = Chain_weight{i}(j)/sum(Chain_weight{i});
end

% Closeness Centrality
for j=1:chainNum
    Closeness_centrality{i}(j) = Chain_centrality{i}(j)*Chain_degree{i}(j)/sum(Chain_Adjacency_matrices_distance{i}(j,:));
    if any(isnan(Closeness_centrality{i}(j)))
        Closeness_centrality{i}(j) = 0;
    end
end

% EigenVector Centrality
[eigen_vectors{i},eigen_values{i}] = eig(Chain_Adjacency_matrices{i});
lambda = max(max(eigen_values{i}));
for j=1:chainNum
    for l=1:chainNum
        Eigenvector_centrality{i}(j) = (1/lambda)*sum(Chain_Adjacency_matrices{i}(j,l)*eigen_vectors{i}(:,l));
        Eigenvector_centrality_weighted{i}(j) = (1/lambda)*sum(Chain_Adjacency_matrices_Weighted{i}(j,l)*eigen_vectors{i}(:,l));
    end
end
Chain_centrality_normalized{i} = Chain_centrality{i}/sum(Chain_centrality{i});
Weight_centrality_normalized{i} = Weight_centrality{i}/sum(Weight_centrality{i});
Closeness_centrality_normalized{i} = Closeness_centrality{i}/sum(Closeness_centrality{i});
Eigenvector_centrality_normalized{i} = Eigenvector_centrality{i}/sum(Eigenvector_centrality{i});
Eigenvector_centrality_weighted_normalized{i} = Eigenvector_centrality_weighted{i}/sum(Eigenvector_centrality_weighted{i});
for j=1:chainNum
if any(isnan(Eigenvector_centrality_normalized{i}(j)))
Eigenvector_centrality_normalized{i}(j)=0;
end
if any(isnan(Eigenvector_centrality_weighted_normalized{i}(j)))
Eigenvector_centrality_weighted_normalized{i}(j)=0;
end
end

% Average Clustering Coefficient
clustering_coef{i} = clustering_coef_bu(Chain_Adjacency_matrices{i});
clustering_coef_arry(i) = clustering_coef{i};
% Global Efficiency of the graph
Efficiency{i} = global_efficiency(Chain_Adjacency_matrices{i});
global_Efficiency(i) = Efficiency{i};
% graph k core number
k_core(i) = find_k_core(Chain_Adjacency_matrices{i});

% degree_assortativity_coefficient
[x, y] = compute_degrees(Chain_Adjacency_matrices{i});
exy = compute_edge_fraction_matrix(Chain_Adjacency_matrices{i}, x, y);
[ax, by] = compute_fraction_sets(exy);
[sigma_a, sigma_b] = compute_strand_deviations(ax, by);
network_degree_assortativity_coefficient(i) = degree_assortativity_coefficient(exy, ax, by, sigma_a, sigma_b);
% degree_assortativity_coefficient

% Network Desity
network_density(i) = calculate_network_density(Chain_Adjacency_matrices{i});

Chain_Adjacency_matrices_arr(:,:,i) = Chain_Adjacency_matrices{i};

end
network_degree_assortativity_coefficient_abs = abs(network_degree_assortativity_coefficient);
network_degree_assortativity_coefficient_norm = network_degree_assortativity_coefficient/max(network_degree_assortativity_coefficient_abs);

degree_average = mean(Graph_average_degree);
clustering_coef_average = mean(clustering_coef_arry);
network_density_average = mean(network_density);
global_Efficiency_average = mean(global_Efficiency);
k_core_average = mean(k_core);
network_degree_assortativity_coefficient_norm_average = mean(network_degree_assortativity_coefficient_norm);

for i = 1:numel(ChainPair)
    adjacency_array(:,:,i) = Chain_Adjacency_matrices{i};
end
    
boxx = int2str(boxsize);
chn = int2str(chainNum);
edgfile = strcat('wt_edgeList_box',boxx,'_chain',chn,'_3k_',temp,'K.xlsx');
degreefile = strcat('wt_degList_box',boxx,'_chain',chn,'_3k_',temp,'K.xlsx');
degreeAdjfile = strcat('wt_degAdjacency_box',boxx,'_chain',chn,'_3k_',temp,'K.xlsx');

sum_adjacency_array = sum(adjacency_array , 3);
csvwrite(degreeAdjfile,sum_adjacency_array);

threshold_2 = 0;
cumulated_adjacency_matrix = sum_adjacency_array > threshold_2;
adjacency_matrix_weighted = sum_adjacency_array.*cumulated_adjacency_matrix;

network=graph(adjacency_matrix_weighted);
% plot(network);
all_edges = network.Edges;
writetable(all_edges,edgfile);
%plot(network);
degree_list(:,2) = sum(cumulated_adjacency_matrix)';
degree_list(:,1) = 1:chainNum;
writematrix(degree_list,degreefile);
%%
for i=1:numel(ChainPair)
    
    Chain_centrality_avg (i,1) = i;
    Chain_centrality_avg (i,2) = mean(Chain_centrality{i});
    csvwrite(strcat('ChainCentrality_avg','_box',boxx,'_chain',chn,'_wt_3k_',temp,'K.csv'),Chain_centrality_avg);
    
    Chain_degree_avg(i,1) = i;
    Chain_degree_avg(i,2) = mean(Chain_degree{i});
    csvwrite(strcat('ChainDeg_avg','_box',boxx,'_chain',chn,'_wt_3k_',temp,'K.csv'),Chain_degree_avg);
    
    Chain_weight_avg (i,1) = i;
    Chain_weight_avg (i,2) = mean(Chain_weight{i});
    csvwrite(strcat('ChainWeight_avg','_box',boxx,'_chain',chn,'_wt_3k_',temp,'K.csv'),Chain_weight_avg);
    
    Closeness_centrality_avg(i,1) = i;
    Closeness_centrality_avg(i,2) = mean(Closeness_centrality{i});
    csvwrite(strcat('ClosenessCentrality_avg','_box',boxx,'_chain',chn,'_wt_3k_',temp,'K.csv'),Closeness_centrality_avg);
    
    clustering_coef_vector(i,1) = i;
    clustering_coef_vector(i,2) = clustering_coef{i};
    csvwrite(strcat('clustCoef_vector','_box',boxx,'_chain',chn,'_wt_3k_',temp,'K.csv'),clustering_coef_vector);
    
    Eigenvector_centrality_avg(i,1) = i;
    if any(isnan(Eigenvector_centrality{i}))
        Eigenvector_centrality{i} = 0;
    end
        Eigenvector_centrality_avg(i,2) = mean(Eigenvector_centrality{i});
        csvwrite(strcat('EigenvectorCentrality_avg','_box',boxx,'_chain',chn,'_wt_3k_',temp,'K.csv'),Eigenvector_centrality_avg);
    
    global_Efficiency_vector(i,1) = i;
    global_Efficiency_vector(i,2) = global_Efficiency(1,i);
    csvwrite(strcat('globalEfficiency_vector','_box',boxx,'_chain',chn,'_wt_3k_',temp,'K.csv'),global_Efficiency_vector);
    

    network_degree_assortativity_coefficient_vector(i,1) = i;
    if any(isnan(network_degree_assortativity_coefficient(1,i)))
        network_degree_assortativity_coefficient(1,i) = 0;
    end
        network_degree_assortativity_coefficient_vector(i,2) = network_degree_assortativity_coefficient(1,i);  
    csvwrite(strcat('networkDeg_assortativityCoeff_vector','_box',boxx,'_chain',chn,'_wt_3k_',temp,'K.csv'),network_degree_assortativity_coefficient_vector);

    network_density_vector(i,1) = i;
    network_density_vector(i,2) = network_density(1,i);
    csvwrite(strcat('networkDensity_vector','_box',boxx,'_chain',chn,'_wt_3k_',temp,'K.csv'),network_density_vector);
    
    Weight_centrality_vector(i,1) = i;
    if any(isnan(Weight_centrality{i}))
       Weight_centrality{i} = 0;
    end
    Weight_centrality_vector(i,2) = mean(Weight_centrality{i});
    csvwrite(strcat('WeightCentrality_vector','_box',boxx,'_chain',chn,'_wt_3k_',temp,'K.csv'),Weight_centrality_vector);
    
end
means_arr(1,1) = mean(Chain_centrality_avg(:,2));
means_arr(2,1) = mean(Chain_degree_avg(:,2));
means_arr(3,1) = mean(Chain_weight_avg(:,2));
means_arr(4,1) = mean(Closeness_centrality_avg(:,2));
means_arr(5,1) = mean(clustering_coef_vector(:,2));
means_arr(6,1) = mean(Eigenvector_centrality_avg(:,2));
means_arr(7,1) = mean(global_Efficiency_vector(:,2));
means_arr(8,1) = mean(network_degree_assortativity_coefficient_vector(:,2));
means_arr(9,1) = mean(network_density_vector(:,2));
means_arr(10,1) = mean(Weight_centrality_vector(:,2));

labels = {'Chain_centrality_avg'; 'Chain_degree_avg'; 'Chain_weight_avg';'Closeness_centrality_avg';'clustering_coef_vector';'Eigenvector_centrality_avg';'global_Efficiency_vector';'network_degree_assortativity_coefficient_vector';'network_density_vector';'Weight_centrality_vector'};
fileName_output = strcat('avg_ntwrk_properties','_box',boxx,'_chain',chn,'_wt_3k_',temp,'K.csv');
fileID = fopen(fileName_output, 'w');
for i = 1:numel(labels)
    fprintf(fileID, '%s,%d\n', labels{i}, means_arr(i));
end
fclose(fileID);

end
%%%%%%%%

% Write files
% Variable_table = readcell('Variable_list.xlsx');
% for i = 1:numel(Variable_table)
%     variable_cell = eval(Variable_table{i});
%     if iscell(variable_cell)
%     variable = cell2table(variable_cell);
%     writetable(variable,strcat(Variable_table{i},'_wt_unpack_alpha0.925_3k_',temp,'K.csv'));
% elseif strcmp(class(variable_cell), 'double')
%     variable = variable_cell;
%     csvwrite(strcat(Variable_table{i},'_wt_unpack_alpha0.925_3k_',temp,'K.csv'),variable);
% end
% end

%% 
%for i = 1:numel(Chain_Adjacency_matrices)
 %   longest_path_arr(i,1) = i;
 %   longest_path_arr(i,2) = longest_path_dag(Chain_Adjacency_matrices{i});
%end
%longest_path_ = max(longest_path_arr(:,2));
%longest_path_arr_normalized = longest_path_arr(:,2)/longest_path_;
%writematrix(longest_path_arr,'wt_unpack_alpha0.925_3k_linked_chains.xlsx');
%%
% Create a figure
%figure('Position', [100, 100, 1000, 300]);
%
%% Create a tiled layout with minimal gaps
%t = tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
%
%% Create the first tile for the scatter plot
%ax1 = nexttile;
%scatter(longest_path_arr(:, 1), longest_path_arr(:, 2), [], [0.62, 0.0, 1.0], 'MarkerFaceColor','none', 'SizeData', 9); % Scatter plot
%title('condensates formation vs time (ns)');
%xlabel('time (ns)');
%ylabel('condensates');
%
%% Set xlim and ylim for the scatter plot
%xlim(ax1, [min(longest_path_arr(:, 1)), max(longest_path_arr(:, 1))]);
%ylim(ax1, [min(longest_path_arr(:, 2))-0.5, max(longest_path_arr(:, 2))+0.5]);
%
%% Create the second tile for the histogram
%ax2 = nexttile;
%bin_width = 0.5; % Specify the bin width
%h = histogram(longest_path_arr(:, 2), 'FaceColor', 'b', 'NumBins', 9, 'BinWidth', bin_width); % Histogram with specified bin width
%title('condensates distribution');
%xlabel('condensates');
%ylabel('Frequency');
%view(ax2, 90, -90); % Rotate the histogram 90 degrees to the right
%
%% Set ylim for the histogram plot to match the scatter plot
%ylim(ax1, [min(longest_path_arr(:, 2))-0.5, max(longest_path_arr(:, 2))+0.5]);
%
%% Export the figure as a TIFF file
%exportgraphics(gcf, 'wt_unpack_alpha0.925_3k_chain_connection.tif', 'Resolution', 300);
%%
%conncomps = {};
% for i = 1:1000
%     A = Chain_Adjacency_matrices{i};
%     g = graph(A);
%     conncomps{i} = conncomp(g);
%     conncomps_m(:,i) = conncomps{i};
%     conncomps_max = max(conncomps_m);
% end
% [conncomps_min,loc] = min(conncomps_max);
% ggg = Chain_Adjacency_matrices{loc-1};
% gg = graph(ggg);
% plot(gg);

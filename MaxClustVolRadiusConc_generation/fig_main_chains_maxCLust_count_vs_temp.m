clear;

% Read system names from the Excel file
Variable_table = readcell('names_temperature.xlsx');

num_of_syst = numel(Variable_table);

% Initialize a cell array to store the data
wt_pack = cell(num_of_syst, 3);
wt_unpack = cell(num_of_syst, 3);

mt_pack = cell(num_of_syst, 3);
mt_unpack = cell(num_of_syst, 3);

% Loop through the trajectories and system names to read data files
for i = 1:3
    trajNum = int2str(i);
    for j = 1:num_of_syst
        wt_pack{j, i} = readtable(strcat('wt_', trajNum, 'C_60chains_clust_size_Num_1k_', Variable_table{j}, '.xlsx'));
        disp(['Loaded data1 for ', Variable_table{j}, ' in trajectory ', trajNum]);

        wt_unpack{j, i} = readtable(strcat('wt_unpack_', trajNum, 'C_60chains_clust_size_Num_1k_', Variable_table{j}, '.xlsx'));
        disp(['Loaded data2 for ', Variable_table{j}, ' in trajectory ', trajNum]);

        mt_pack{j, i} = readtable(strcat('mt_alpha0.925_', trajNum, 'C_60chains_clust_size_Num_1k_', Variable_table{j}, '.xlsx'));
        disp(['Loaded data1 for ', Variable_table{j}, ' in trajectory ', trajNum]);

        mt_unpack{j, i} = readtable(strcat('mt_unpack_', trajNum, 'C_60chains_clust_size_Num_1k_', Variable_table{j}, '.xlsx'));
        disp(['Loaded data2 for ', Variable_table{j}, ' in trajectory ', trajNum]);

    end
end

wt_pack_clustsize = cell(num_of_syst, 3);
wt_unpack_clustsize = cell(num_of_syst, 3);
wt_pack_totalclust = cell(num_of_syst, 3);
wt_unpack_totalclust = cell(num_of_syst, 3);

mt_pack_clustsize = cell(num_of_syst, 3);
mt_unpack_clustsize = cell(num_of_syst, 3);
mt_pack_totalclust = cell(num_of_syst, 3);
mt_unpack_totalclust = cell(num_of_syst, 3);

for j = 1:num_of_syst
    for i = 1:3
        if ~isempty(wt_pack{j, i})
            wt_pack_clustsize{j, i} = table2array(wt_pack{j, i}(:, 3));
            wt_pack_totalclust{j, i} = table2array(wt_pack{j, i}(:, 2));
            mt_pack_clustsize{j, i} = table2array(mt_pack{j, i}(:, 3));
            mt_pack_totalclust{j, i} = table2array(mt_pack{j, i}(:, 2));
        end
        if ~isempty(wt_unpack{j, i})
            wt_unpack_clustsize{j, i} = table2array(wt_unpack{j, i}(:, 3));
            wt_unpack_totalclust{j, i} = table2array(wt_unpack{j, i}(:, 2));
            mt_unpack_clustsize{j, i} = table2array(mt_unpack{j, i}(:, 3));
            mt_unpack_totalclust{j, i} = table2array(mt_unpack{j, i}(:, 2));
        end
    end
end

mean_last500_wt_pack_clustsize = zeros(num_of_syst, 1);
std_last500_wt_pack_clustsize = zeros(num_of_syst, 1);
mean_last500_wt_unpack_clustsize = zeros(num_of_syst, 1);
std_last500_wt_unpack_clustsize = zeros(num_of_syst, 1);

mean_last500_wt_pack_totalclust = zeros(num_of_syst, 1);
std_last500_wt_pack_totalclust = zeros(num_of_syst, 1);
mean_last500_wt_unpack_totalclust = zeros(num_of_syst, 1);
std_last500_wt_unpack_totalclust = zeros(num_of_syst, 1);

mean_last500_mt_pack_clustsize = zeros(num_of_syst, 1);
std_last500_mt_pack_clustsize = zeros(num_of_syst, 1);
mean_last500_mt_unpack_clustsize = zeros(num_of_syst, 1);
std_last500_mt_unpack_clustsize = zeros(num_of_syst, 1);

mean_last500_mt_pack_totalclust = zeros(num_of_syst, 1);
std_last500_mt_pack_totalclust = zeros(num_of_syst, 1);
mean_last500_mt_unpack_totalclust = zeros(num_of_syst, 1);
std_last500_mt_unpack_totalclust = zeros(num_of_syst, 1);

for j = 1:num_of_syst
    wt_pack_clust_means = zeros(3, 1);
    wt_unpack_clust_means = zeros(3, 1);
    wt_pack_total_means = zeros(3, 1);
    wt_unpack_total_means = zeros(3, 1);

    mt_pack_clust_means = zeros(3, 1);
    mt_unpack_clust_means = zeros(3, 1);
    mt_pack_total_means = zeros(3, 1);
    mt_unpack_total_means = zeros(3, 1);

    for i = 1:3
        if ~isempty(wt_pack_clustsize{j, i})
            wt_pack_clust_means(i) = mean(wt_pack_clustsize{j, i}(end-499:end));
            mt_pack_clust_means(i) = mean(mt_pack_clustsize{j, i}(end-499:end));
        end
        if ~isempty(wt_unpack_clustsize{j, i})
            wt_unpack_clust_means(i) = mean(wt_unpack_clustsize{j, i}(end-499:end));
            mt_unpack_clust_means(i) = mean(mt_unpack_clustsize{j, i}(end-499:end));
        end
        if ~isempty(wt_pack_totalclust{j, i})
            wt_pack_total_means(i) = mean(wt_pack_totalclust{j, i}(end-499:end));
            mt_pack_total_means(i) = mean(mt_pack_totalclust{j, i}(end-499:end));
        end
        if ~isempty(wt_unpack_totalclust{j, i})
            wt_unpack_total_means(i) = mean(wt_unpack_totalclust{j, i}(end-499:end));
            mt_unpack_total_means(i) = mean(mt_unpack_totalclust{j, i}(end-499:end));
        end
    end

    mean_last500_wt_pack_clustsize(j) = mean(wt_pack_clust_means);
    std_last500_wt_pack_clustsize(j) = std(wt_pack_clust_means);
    mean_last500_wt_unpack_clustsize(j) = mean(wt_unpack_clust_means);
    std_last500_wt_unpack_clustsize(j) = std(wt_unpack_clust_means);

    mean_last500_wt_pack_totalclust(j) = mean(wt_pack_total_means);
    std_last500_wt_pack_totalclust(j) = std(wt_pack_total_means);
    mean_last500_wt_unpack_totalclust(j) = mean(wt_unpack_total_means);
    std_last500_wt_unpack_totalclust(j) = std(wt_unpack_total_means);


    mean_last500_mt_pack_clustsize(j) = mean(mt_pack_clust_means);
    std_last500_mt_pack_clustsize(j) = std(mt_pack_clust_means);
    mean_last500_mt_unpack_clustsize(j) = mean(mt_unpack_clust_means);
    std_last500_mt_unpack_clustsize(j) = std(mt_unpack_clust_means);
    
    mean_last500_mt_pack_totalclust(j) = mean(mt_pack_total_means);
    std_last500_mt_pack_totalclust(j) = std(mt_pack_total_means);
    mean_last500_mt_unpack_totalclust(j) = mean(mt_unpack_total_means);
    std_last500_mt_unpack_totalclust(j) = std(mt_unpack_total_means);

end
%%
figure('Position', [100, 100, 650, 250]);  % Adjust figure size
tiledlayout(1, 2, 'Padding', 'none', 'TileSpacing', 'tight');

temp_labels = 270:10:340;  % Actual temperature values

% Max Cluster Size Plot
nexttile;
errorbar(temp_labels, mean_last500_wt_pack_clustsize, std_last500_wt_pack_clustsize, '-o', 'DisplayName', 'compact to PS (WT)', 'Color', [0 0.45 0.74], 'LineWidth', 1); hold on;
errorbar(temp_labels, mean_last500_wt_unpack_clustsize, std_last500_wt_unpack_clustsize, '-o', 'DisplayName', 'dispersed to PS (WT)', 'Color',    [0.85 0.33 0.1], 'LineWidth', 1); hold on;

errorbar(temp_labels, mean_last500_mt_pack_clustsize, std_last500_mt_pack_clustsize, '-^', 'DisplayName', 'compact to PS (MT)', 'Color', [0.47 0.67 0.19], 'LineWidth', 1); hold on;
errorbar(temp_labels, mean_last500_mt_unpack_clustsize, std_last500_mt_unpack_clustsize, '-^', 'DisplayName', 'dispersed to PS (MT)', 'Color',    [0.64 0.08 0.18], 'LineWidth', 1);

xlabel('Temperature (K)', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
ylabel('Max Cluster size', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex');
xlim([265 345]);  % Now it works correctly
legend('box', 'off', 'Location', 'best', 'FontSize', 6, 'Interpreter', 'latex');
box off;

% Total Cluster Count Plot
nexttile;
errorbar(temp_labels, mean_last500_wt_pack_totalclust, std_last500_wt_pack_totalclust, '-o', 'DisplayName', 'compact to PS (WT)', 'Color', [0 0.45 0.74] , 'LineWidth', 1)   ; hold on;
errorbar(temp_labels, mean_last500_wt_unpack_totalclust, std_last500_wt_unpack_totalclust, '-o', 'DisplayName', 'dispersed to PS (WT)', 'Color', [0.85 0.33 0.1] , 'LineWidth', 1); hold on;

errorbar(temp_labels, mean_last500_mt_pack_totalclust, std_last500_mt_pack_totalclust, '-^', 'DisplayName', 'compact to PS (MT)', 'Color', [0.47 0.67 0.19], 'LineWidth', 1)   ; hold on;
errorbar(temp_labels, mean_last500_mt_unpack_totalclust, std_last500_mt_unpack_totalclust, '-^', 'DisplayName', 'dispersed to PS (MT)', 'Color',    [0.64 0.08 0.18] , 'LineWidth', 1) ;

xlabel('Temperature (K)', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
ylabel('Total Clusters', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex');
xlim([265 345]);  % Still valid
box off;

 print(gcf, 'wt_vs_mt_chains_maxClust_count_means_vs_temp.tif', '-dtiff', '-r300');

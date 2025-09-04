clear;
fileName = 'mute_';

% Read system names from the Excel file
Variable_table = readcell('names_temperature.xlsx');

num_of_syst = numel(Variable_table);
numPlotsPerSystem = 4; % Four plots per system (ClustSize First50, ClustSize Last50, TotalClust First50, TotalClust Last50)
totalSubplots = num_of_syst * numPlotsPerSystem;

% Initialize a cell array to store the data
data_pack = cell(num_of_syst, 3);
data_unpack = cell(num_of_syst, 3);

% Loop through the trajectories and system names to read data files
for i = 1:3
    trajNum = int2str(i);
    for j = 1:num_of_syst
        data_pack{j, i} = readtable(strcat(fileName, 'alpha0.925_', trajNum, 'C_60chains_clust_size_Num_1k_', Variable_table{j}, '.xlsx'));
        disp(['Loaded data1 for ', Variable_table{j}, ' in trajectory ', trajNum]);

        data_unpack{j, i} = readtable(strcat(fileName, 'unpack_alpha0.925_', trajNum, 'C_60chains_clust_size_Num_1k_', Variable_table{j}, '.xlsx'));
        disp(['Loaded data2 for ', Variable_table{j}, ' in trajectory ', trajNum]);
    end
end

pack_framenum = cell(size(data_pack));
pack_clustsize = cell(size(data_pack));
unpack_framenum = cell(size(data_unpack));
unpack_clustsize = cell(size(data_unpack));
pack_totalclust =  cell(size(data_pack));
unpack_totalclust =  cell(size(data_unpack));

for j = 1:num_of_syst
    for i = 1:3
        if ~isempty(data_pack{j, i})
            pack_framenum{j, i} = table2array(data_pack{j, i}(1:1000, 1));
            pack_clustsize{j, i} = table2array(data_pack{j, i}(1:1000, 3));
            pack_totalclust{j, i} = table2array(data_pack{j, i}(1:1000, 2));
        end
        if ~isempty(data_unpack{j, i})
            unpack_framenum{j, i} = table2array(data_unpack{j, i}(1:1000, 1));
            unpack_clustsize{j, i} = table2array(data_unpack{j, i}(1:1000, 3));
            unpack_totalclust{j, i} = table2array(data_unpack{j, i}(1:1000, 2));
        end
    end
end

meanpack_clustsize = cell(num_of_syst, 1);
stderrpack_clustsize = cell(num_of_syst, 1);
meanunpack_clustsize = cell(num_of_syst, 1);
stderrunpack_clustsize = cell(num_of_syst, 1);

meanpack_totalclust = cell(num_of_syst, 1);
stderrpack_totalclust = cell(num_of_syst, 1);
meanunpack_totalclust = cell(num_of_syst, 1);
stderrunpack_totalclust = cell(num_of_syst, 1);

for j = 1:num_of_syst
    combinedpack_clustsize = [pack_clustsize{j, :}];
    combinedunpack_clustsize = [unpack_clustsize{j, :}];
    combinedpack_totalclust = [pack_totalclust{j, :}];
    combinedunpack_totalclust = [unpack_totalclust{j, :}];

    meanpack_clustsize{j} = mean(combinedpack_clustsize, 2);
    stderrpack_clustsize{j} = std(combinedpack_clustsize, 0, 2) / sqrt(3);
    meanunpack_clustsize{j} = mean(combinedunpack_clustsize, 2);
    stderrunpack_clustsize{j} = std(combinedunpack_clustsize, 0, 2) / sqrt(3);

    meanpack_totalclust{j} = mean(combinedpack_totalclust, 2);
    stderrpack_totalclust{j} = std(combinedpack_totalclust, 0, 2) / sqrt(3);
    meanunpack_totalclust{j} = mean(combinedunpack_totalclust, 2);
    stderrunpack_totalclust{j} = std(combinedunpack_totalclust, 0, 2) / sqrt(3);
end
%%
% function plot_shaded_error(x, y, err, color)
%     fill([x; flipud(x)], [y - err; flipud(y + err)], color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     plot(x, y, color, 'LineWidth', 1.5);
% end

function plot_shaded_error(x, y, err, color)
    fill([x; flipud(x)], ...
         [y - err; flipud(y + err)], ...
         color, ...
         'FaceAlpha', 0.3, 'EdgeColor', 'none');  % <-- critical
    plot(x, y, '-', 'Color', color, 'LineWidth', 0.3);
end


% Determine y-limits for cluster size plots
all_y_clust = [];
for j = 1:num_of_syst
    all_y_clust = [all_y_clust; meanpack_clustsize{j}; meanunpack_clustsize{j}];
end
y_lim_clust = [min(all_y_clust) - 0.1 * range(all_y_clust), max(all_y_clust) + 0.1 * range(all_y_clust)];

% Determine y-limits for total cluster number plots
all_y_total = [];
for j = 1:num_of_syst
    all_y_total = [all_y_total; meanpack_totalclust{j}; meanunpack_totalclust{j}];
end
y_lim_total = [min(all_y_total) - 0.1 * range(all_y_total), max(all_y_total) + 0.1 * range(all_y_total)];
% %
figure('Position', [100, 100, 500, 300]);  % 4 rows × 2 columns
tiledlayout(3, 3, 'Padding', 'tight', 'TileSpacing', 'tight');

% Define custom colors
col_compact_max   = [0    0.45 0.74];   % Blue
col_dispersed_max = [0.85 0.33 0.10];   % Orange
col_compact_total = [0.47 0.67 0.19];   % Green
col_dispersed_total = [0.64 0.08 0.18]; % Dark Red

for j = 1:num_of_syst
    nexttile(j);

    x1 = pack_framenum{j, 1};
    y1_max = meanpack_clustsize{j};
    err1_max = stderrpack_clustsize{j};
    y1_total = meanpack_totalclust{j};
    err1_total = stderrpack_totalclust{j};

    x2 = unpack_framenum{j, 1};
    y2_max = meanunpack_clustsize{j};
    err2_max = stderrunpack_clustsize{j};
    y2_total = meanunpack_totalclust{j};
    err2_total = stderrunpack_totalclust{j};

    % LEFT Y-axis: Max Cluster
    yyaxis left
    plot_shaded_error(x1, y1_max, err1_max, col_compact_max); hold on;
    plot_shaded_error(x2, y2_max, err2_max, col_dispersed_max);
    set(gca, 'YColor', 'k');  % Black axis
    ylim(y_lim_clust)

    if mod(j, 3) == 1  % Only for column 1
        ylabel('Max Clust', 'Interpreter', 'latex', 'FontSize', 8);
    else 
        set(gca, 'YTickLabel', []);
    end

    % RIGHT Y-axis: Total Cluster
    yyaxis right
    plot_shaded_error(x1, y1_total, err1_total, col_compact_total); hold on;
    plot_shaded_error(x2, y2_total, err2_total, col_dispersed_total);
    set(gca, 'YColor', 'k');
    ylim(y_lim_total)

    if mod(j, 3) == 0  % Only for column 3
        ylabel('Total Clust', 'Interpreter', 'latex', 'FontSize', 8);
    
    elseif j == 8
    set(gca, 'YTickLabel', {'0','5','10','15'}, 'TickLabelInterpreter', 'latex');
    else
        set(gca, 'YTickLabel', []);
    end

    xlim([min(x1), max(x1)]);
    if j > 6  % bottom row
        xlabel('Time ($\mu$s)', 'Interpreter', 'latex', 'FontSize', 8);
        set(gca, 'XTickLabel', {'2','4','6','8','10'}, 'TickLabelInterpreter', 'latex');
        elseif j == 6
    set(gca, 'XTickLabel', {'2','4','6','8','10'}, 'TickLabelInterpreter', 'latex');
    else
        set(gca, 'XTickLabel', []);
    end

    % Add system label
    text(0.05, 0.9, Variable_table{j}, 'Units', 'normalized', ...
         'FontSize', 8, 'Interpreter','latex', 'Color', [0.1 0.1 0.1]);

    % Legend in the first plot
    % if j == 8
    %     legend({'Compact Max', '', 'Dispersed Max', '', ...
    %             'Compact Total', '', 'Dispersed Total'}, ...
    %            'Location', 'bestoutside', 'Box', 'off', ...
    %            'FontSize', 8, 'Interpreter', 'latex');
    % end

    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 8);
    box off;
end

% After all tiles are plotted
lgd = legend({'compact to PS Max clust', '', 'dispersed to PS Max clust','', ...
    'compact to PS total clust', '', 'dispersed to PS total clust'}, ...
    'Interpreter', 'latex', 'FontSize', 7, ...
    'Box', 'off', ...
    'Orientation', 'vertical');


% Position the legend outside the tiled layout (top right)
%lgd.Layout.Tile = 'southeast';  % Moves it above the full figure

set(lgd, 'Units', 'normalized', 'Position', [0.71 0.12 0.1 0.2]);


figName = strcat(fileName, 'chains_maxClust_count_vs_time_temp.tif');
print(gcf, figName, '-dtiff', '-r300');
clear;

% File name prefix
fileName = 'avg_networkProperties_chain60_';

% Read system names from the Excel file
temp_table = readcell('names_temperature.xlsx');

% Read property names from the Excel file
properties_table = readcell('names_properties.xlsx');

% Initialize arrays
num_of_temp = numel(temp_table);
num_of_properties = 9;

% Initialize data arrays
data_array1 = nan(num_of_temp, num_of_properties + 1);
data_array2 = nan(num_of_temp, num_of_properties + 1);
data_array3 = nan(num_of_temp, num_of_properties + 1);
data_array4 = nan(num_of_temp, num_of_properties + 1);

for j = 1:num_of_temp
    temperature = regexp(temp_table{j}, '\d+', 'match');
    
    if ~isempty(temperature)
        temperature_value = round(str2double(temperature{1}));
        data_array1(j, 1) = temperature_value;
        data_array2(j, 1) = temperature_value;
        data_array3(j, 1) = temperature_value;
        data_array4(j, 1) = temperature_value;
    end

    fullFileName1 = strcat(fileName, 'wt_3k_end5mics_', temp_table{j}, '.csv');
    fullFileName2 = strcat(fileName, 'wt_unpack_alpha0.925_3k_end5mics_', temp_table{j}, '.csv');
    fullFileName3 = strcat(fileName, 'mt_3k_end5mics_', temp_table{j}, '.csv');
    fullFileName4 = strcat(fileName, 'mt_unpack_3k_end5mics_', temp_table{j}, '.csv');

    try
        table_data = {
            readtable(fullFileName1);
            readtable(fullFileName2);
            readtable(fullFileName3);
            readtable(fullFileName4);
        };

        for i = 1:num_of_properties
            property_name = properties_table{i, 1};
            for k = 1:4
                property_index = find(strcmp(table_data{k}{:, 1}, property_name));
                if ~isempty(property_index)
                    val = table_data{k}{property_index, 2};
                    if iscell(val), val = cell2mat(val); end
                else
                    val = NaN;
                end
                eval(sprintf('data_array%d(j, i + 1) = val;', k));
            end
        end
    catch
        warning('Missing one or more files at temperature %s', temp_table{j});
        data_array1(j, 2:end) = NaN;
        data_array2(j, 2:end) = NaN;
        data_array3(j, 2:end) = NaN;
        data_array4(j, 2:end) = NaN;
    end
end
%%
% Table headers
column_headers = ['temperature', properties_table(1:num_of_properties, 1)'];
data_table1 = array2table(data_array1, 'VariableNames', column_headers);
data_table2 = array2table(data_array2, 'VariableNames', column_headers);
data_table3 = array2table(data_array3, 'VariableNames', column_headers);
data_table4 = array2table(data_array4, 'VariableNames', column_headers);

subplot_labels = 'A':'Z';  % Predefine subplot labels (sufficient for up to 26 subplots)

% Plotting
figure('Position', [100, 100, 800, 800]);
custom_y_labels = {'avg chain centrality', 'avg chain degree', 'avg chain weight', 'avg closeness centrality', 'clustering coefficient', 'avg eigenvector centrality', 'global efficiency', 'degree assortativity coeff', 'density'};

plot_height = 0.7 / 3;
plot_width = 0.7 / 3;
vertical_gap = -0.01;
horizontal_gap = 0.07;

for i = 1:num_of_properties
    [row, col] = ind2sub([3, 3], i);
    left_margin = 0.1;
    bottom_margin = 0.9 - row * plot_height - (row - 1) * vertical_gap;
    width = plot_width;
    height = plot_height - 0.05;

    subplot('Position', [left_margin + (col - 1) * (width + horizontal_gap), bottom_margin, width, height]);

    plot(data_table1.temperature, data_table1{:, i + 1}, '-o', 'DisplayName', 'compact to PS (WT)', 'Color', [0 0.45 0.74]); hold on;
    plot(data_table2.temperature, data_table2{:, i + 1}, '-o', 'DisplayName', 'dispersed to PS', 'Color', [0.85 0.33 0.1]);
    plot(data_table3.temperature, data_table3{:, i + 1}, '-^', 'DisplayName', 'compact to PS (MT)', 'Color', [0.47 0.67 0.19]);
    plot(data_table4.temperature, data_table4{:, i + 1}, '-^', 'DisplayName', 'dispersed to PS', 'Color', [0.64 0.08 0.18]);

    xlim([265 345]);
    xticks(270:10:340);
    if ismember(i, [3, 6, 9])
        xlabel('temperature (K)', 'FontSize', 12, 'Interpreter','latex');
    else
        set(gca, 'XTickLabel', []);
    end

    if i == 6
        yticklabels({'-1', '-0.05', '0', '0.05', '1'});
        %text(0.11, 1.08, '×10^{-2}', 'Units', 'normalized', 'FontSize', 9);
    elseif i == 8
       text(0.0, 1.08, '$\times 10^{5}$', 'Units', 'normalized', 'FontSize', 9, 'Interpreter','latex');

    end
   set(gca, 'TickLabelInterpreter', 'latex');
    ylabel(custom_y_labels{i}, 'Interpreter', 'latex', 'FontSize', 12);

    % Add subplot label (e.g., A, B, C, ...)
    text(-0.15, 1.1, subplot_labels(i), 'Units', 'normalized', ...
         'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'latex');
    
    if i == 1
        legend('box', 'off', 'Location', 'best', 'FontSize', 6, 'Interpreter', 'latex');
    end
end

% Save figure
 figName1 = strcat(fileName, '_vs_temp_allSystems.tif');
  exportgraphics(gcf, figName1, 'Resolution', 300);

 %% Correlation heatmaps among systems for each property
systems_data = cat(3, data_array1(:, 2:end), data_array2(:, 2:end), ...
                      data_array3(:, 2:end), data_array4(:, 2:end));

system_labels = {'WT-Compact', 'WT-Dispersed', 'MT-Compact', 'MT-Dispersed'};

% Define tiled layout
figure('Position', [100 100 800 800]);
t = tiledlayout(3, 3, 'Padding', 'tight', 'TileSpacing', 'compact');

% Loop through properties
for i = 1:num_of_properties
    nexttile(t, i);

    % Compute correlation
    % (Assuming `corr_data` is a matrix of size [num_systems x num_properties])
    R = corrcoef([data_array1(:, i+1), data_array2(:, i+1), ...
                  data_array3(:, i+1), data_array4(:, i+1)], ...
                  'Rows', 'pairwise');

    % Create heatmap without axis handle (each heatmap manages its own axes)
    h = heatmap(system_labels, system_labels, R, ...
        'Colormap', jet, ...
        'ColorLimits', [-1 1], ...
        'CellLabelColor','none', ...
        'MissingDataColor', [1 1 1], ...
        'MissingDataLabel', '');

    % Hide colorbars except one
    if i ~= 9
        h.ColorbarVisible = 'off';
    end

    % Tick labels only on left column and bottom row
    if mod(i-1, 3) ~= 0
        h.YDisplayLabels = repmat({''}, size(system_labels));
    end
    if i <= 6
        h.XDisplayLabels = repmat({''}, size(system_labels));
    end

    % Title as property label
    h.Title = custom_y_labels{i};
end

% % Create dummy image for shared colorbar
% cb_ax = axes('Position', [0.93, 0.2, 0.015, 0.6]); % [left, bottom, width, height]
% imagesc(cb_ax, linspace(-1, 1, 100)');
% colormap(cb_ax, parula);
% set(cb_ax, 'YTick', [-1, 0, 1], 'XTick', [], 'YAxisLocation', 'right', 'Box', 'off');
% ylabel(cb_ax, 'Correlation', 'FontSize', 10, 'Interpreter', 'latex');


% Save figure as TIFF
corr_fig_name = strcat(fileName, '_correlation_across_systems.tif');
exportgraphics(gcf, corr_fig_name, 'Resolution', 300);

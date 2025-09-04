clear;

% File name prefix and temperature table
fileName = 'chain_residue_connections_';
temp_table = readcell('names_temperature.xlsx');
num_of_temp = numel(temp_table);
residue = 33;

% Preallocate arrays
WT_compact = nan(residue, num_of_temp);
WT_dispersed = nan(residue, num_of_temp);
MT_compact = nan(residue, num_of_temp);
MT_dispersed = nan(residue, num_of_temp);
names_array_WT = cell(residue, num_of_temp);
names_array_MT = cell(residue, num_of_temp);
temperature_array = zeros(1, num_of_temp);

% Load all 4 data sets
for j = 1:num_of_temp
    temperature = regexp(temp_table{j}, '\d+', 'match');
    temperature_array(j) = str2double(temperature{1});

    f1 = strcat(fileName, 'wt_', temp_table{j}, '.csv');
    f2 = strcat(fileName, 'wt_unpack_', temp_table{j}, '.csv');
    f3 = strcat(fileName, 'mt_', temp_table{j}, '.csv');
    f4 = strcat(fileName, 'mt_unpack_', temp_table{j}, '.csv');

    t1 = readtable(f1);
    t2 = readtable(f2);
    t3 = readtable(f3);
    t4 = readtable(f4);

    WT_compact(:, j) = t1{:, 2};
    WT_dispersed(:, j) = t2{:, 2};
    MT_compact(:, j) = t3{:, 2};
    MT_dispersed(:, j) = t4{:, 2};

    names_array_WT(:, j) = t1{:, 1};
    names_array_MT(:, j) = t3{:, 1};
end

%% %
unique_residue_labels_WT = strcat(names_array_WT(:, 1), "\_", string(1:residue)');
unique_residue_labels_MT = strcat(names_array_MT(:, 1), "\_", string(1:residue)');
%
% === Plotting ===
figure('Position', [100 100 900 500]);
tiledlayout(1, 4, 'Padding', 'compact', 'TileSpacing', 'tight');

% WT Compact
nexttile;
h1 = heatmap(string(temperature_array), unique_residue_labels_WT, WT_compact, ...
       'Interpreter', 'latex', 'Colormap', jet, 'ColorbarVisible', 'off');
xlabel('temperature (K)'); ylabel('Residue');


% Add LaTeX title manually as annotation
annotation('textbox', [0.02 0.95 0.3 0.05], 'String', ...
    'compact to PS', ...
    'Interpreter', 'latex', 'FontSize', 12, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');

% Add subplot label 'A'
annotation('textbox', [0.07 0.91 0.05 0.05], 'String', 'A', ...
    'Interpreter', 'latex', 'FontSize', 16, 'EdgeColor', 'none');


% WT Dispersed
nexttile;
h2 = heatmap(string(temperature_array), unique_residue_labels_WT, WT_dispersed, ...
      'Interpreter', 'latex', 'Colormap', jet, 'ColorbarVisible', 'off');
xlabel('temperature (K)');
h2.YDisplayLabels = repmat({''}, size(unique_residue_labels_WT));

% Add LaTeX title manually as annotation
annotation('textbox', [0.25 0.95 0.3 0.05], 'String', ...
    'dispersed to PS', ...
    'Interpreter', 'latex', 'FontSize', 12, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');


% MT Compact
nexttile;
h3 = heatmap(string(temperature_array), unique_residue_labels_MT, MT_compact, ...
         'Interpreter', 'latex', 'Colormap', jet, 'ColorbarVisible', 'off');
xlabel('temperature (K)');

% Add LaTeX title manually as annotation
annotation('textbox', [0.47 0.95 0.3 0.05], 'String', ...
    'compact to PS', ...
    'Interpreter', 'latex', 'FontSize', 12, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');

% Add subplot label 'B'
annotation('textbox', [0.52 0.91 0.05 0.05], 'String', 'B', ...
    'Interpreter', 'latex', 'FontSize', 16, 'EdgeColor', 'none');


% MT Dispersed
nexttile;
h4 = heatmap(string(temperature_array), unique_residue_labels_MT, MT_dispersed, ...
         'Interpreter', 'latex', 'Colormap', jet, 'ColorbarVisible', 'on');
xlabel('temperature (K)');
h4.YDisplayLabels = repmat({''}, size(unique_residue_labels_MT));

% Add LaTeX title manually as annotation
annotation('textbox', [0.7 0.95 0.3 0.05], 'String', ...
    'dispersed to PS', ...
    'Interpreter', 'latex', 'FontSize', 12, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');


% Sync color limits across all
clims = [min([h1.ColorLimits(1), h2.ColorLimits(1), h3.ColorLimits(1), h4.ColorLimits(1)]), ...
         max([h1.ColorLimits(2), h2.ColorLimits(2), h3.ColorLimits(2), h4.ColorLimits(2)])];

h1.ColorLimits = clims;
h2.ColorLimits = clims;
h3.ColorLimits = clims;
h4.ColorLimits = clims;

% Reduce the width of the heatmap axes (which also squeezes the colorbar)
% Get the underlying axes that contains the heatmap
ax2 = struct(h2).Axes;  % h2 is your heatmap object

% Now safely change the width to reduce the colorbar size
pos = ax2.Position;
pos(3) = pos(3);  % shrink width by 0.05 (adjust as needed)
ax2.Position = pos;


% Save figure
figName = strcat(fileName, 'vs_temp_WT_MT.tif');
 %saveas(gcf, figName); % Uncomment to save
 print(gcf, figName, '-dtiff', '-r300');

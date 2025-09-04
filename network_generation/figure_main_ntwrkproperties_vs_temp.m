clear;

% File name prefix
fileName = 'avg_networkProperties_';

% Read system names from the Excel file
temp_table = readcell('names_temperature.xlsx');

% Read property names from the Excel file
properties_table = readcell('names_properties.xlsx');

% Initialize the temperature_array
temperature_array = [];

% Get the number of elements in the temp_table
num_of_temp = numel(temp_table);
num_of_properties = 9; % Assuming there are 9 properties

% Initialize the data arrays to store temperatures and properties
data_array1 = nan(num_of_temp, num_of_properties + 1);
data_array2 = nan(num_of_temp, num_of_properties + 1);

% Loop through each element in the temp_table
for j = 1:num_of_temp
    % Extract numerical values from the current variable name
    temperature = regexp(temp_table{j}, '\d+', 'match');
    
    % Convert the extracted numeric strings to integer values and store them in temperature_array
    if ~isempty(temperature)
        temperature_value = round(str2double(temperature{1}));
        data_array1(j, 1) = temperature_value; % Store temperature in the first column
        data_array2(j, 1) = temperature_value; % Store temperature in the first column
    else
        data_array1(j, 1) = NaN; % Handle cases where no numeric value is found
        data_array2(j, 1) = NaN; % Handle cases where no numeric value is found
    end
    
    % Construct the full file names
  
    fullFileName1 = strcat(fileName, 'wt_3k_end5mics_', temp_table{j}, '.csv');
    fullFileName2 = strcat(fileName, 'wt_unpack_3k_end5mics_', temp_table{j}, '.csv');
    
    % Read the table from the CSV files
    try
        table_data1 = readtable(fullFileName1);
        table_data2 = readtable(fullFileName2);
        
        % Loop through each property
        for i = 1:num_of_properties
            % Find the row in the table_data where the property name matches
            property_name = properties_table{i, 1};
            property_index1 = find(strcmp(table_data1{:, 1}, property_name));
            property_index2 = find(strcmp(table_data2{:, 1}, property_name));
            
            if ~isempty(property_index1)
                % Extract the corresponding property value
                property_value1 = table_data1{property_index1, 2};
                % Convert the property value to numeric if it is a cell
                if iscell(property_value1)
                    property_value1 = cell2mat(property_value1);
                end
                
                % Store the property value in the corresponding column
                data_array1(j, i + 1) = property_value1;
            else
                data_array1(j, i + 1) = NaN; % Handle cases where the property name is not found
            end
            
            if ~isempty(property_index2)
                % Extract the corresponding property value
                property_value2 = table_data2{property_index2, 2};
                % Convert the property value to numeric if it is a cell
                if iscell(property_value2)
                    property_value2 = cell2mat(property_value2);
                end
                
                % Store the property value in the corresponding column
                data_array2(j, i + 1) = property_value2;
            else
                data_array2(j, i + 1) = NaN; % Handle cases where the property name is not found
            end
        end
    catch
        warning('Could not read file: %s or %s', fullFileName1, fullFileName2);
        data_array1(j, 2:num_of_properties + 1) = NaN; % Handle cases where the file can't be read
        data_array2(j, 2:num_of_properties + 1) = NaN; % Handle cases where the file can't be read
    end
end

% Define column headers
column_headers = ['temperature', properties_table(1:num_of_properties, 1)'];

% Create tables with the data arrays and column headers
data_table1 = array2table(data_array1, 'VariableNames', column_headers);
data_table2 = array2table(data_array2, 'VariableNames', column_headers);

% Display the final data tables
disp('Data Table 1:');
disp(data_table1);

disp('Data Table 2:');
disp(data_table2);

%%
% Plot the data as separate subplots with gaps
figure('Position', [100, 100, 800, 800]); % Adjust the figure size

% Define custom y-axis labels and titles
custom_y_labels = {'avg chain centrality', 'avg chain degree', 'avg chain weight', 'avg closeness centrality', 'clustering coefficient', 'avg eigenvector centrality', 'global efficiency', 'degree assortativity coeff', 'density'};

num_plots = num_of_properties; % Number of subplots
plot_height = 0.7 / 3; % Height of each subplot
plot_width = 0.7 / 3; % Width of each subplot
vertical_gap = -0.01; % Vertical gap between subplots
horizontal_gap = 0.07; % Horizontal gap between subplots



for i = 1:num_plots
    [row, col] = ind2sub([3, 3], i); % Calculate row and column for subplot
    
    % Calculate position with gaps
    left_margin = 0.1;
    bottom_margin = 0.9 - row * plot_height - (row - 1) * vertical_gap;
    width = plot_width;
    height = plot_height - 0.05;
    
    subplot('Position', [left_margin + (col - 1) * (width + horizontal_gap), bottom_margin, width, height]);
    % Plot data from both data_table1 and data_table2
    plot(data_table1.temperature, data_table1{:, i + 1}, '-o', 'DisplayName', 'Compact to PS');
    hold on;
    plot(data_table2.temperature, data_table2{:, i + 1}, '-s', 'DisplayName', 'Dispersed to PS');
    hold on;
    % Define custom x-ticks
 xlim([265 345]);
xticks([270 280 290 300 310 320 330 340]);
xtick_labels = {'270', '280', '290', '300', '310', '320', '330', '340'};
    if ismember(i, [3, 6, 9]) % Only add labels for these subplots
        set(gca, 'XTick', xticks, 'XTickLabel', xtick_labels, 'FontSize', 10);
        xlabel('Temperature (K)', 'FontSize', 12, 'Interpreter','latex');
    else
        set(gca, 'XTickLabel', []); % Remove x-tick labels for other subplots
    end
    % Add exponent annotation at the top of y-axis for subplots 6 and 8
    if i == 6
       %ytickformat('%.e'); % Format y-ticks in scientific notation for the second subplot
       yticklabels({'-1', '-0.05', '0', '0.05', '1'}); % Hide y-tick labels
       %text(0.11, 1, '×10^{-2}', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
    elseif i == 8
        text(0.1, 1, '×10^{5}', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
    end
    ylabel(custom_y_labels{i}, 'Interpreter', 'latex', 'FontSize', 12);
   % Show legend only in the first subplot
    if i == 4
        % Add legend
    legend('compact to PS', 'dispersed to PS', 'box', 'off', 'Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
    end
    
end

figName1 = strcat(fileName, '_vs_temp_wt.tif');
% Save the figure
exportgraphics(gcf, figName1, 'Resolution', 600);
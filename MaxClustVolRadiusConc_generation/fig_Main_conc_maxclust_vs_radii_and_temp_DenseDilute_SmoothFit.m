clear;

temp_table = readcell('names_temperature.xlsx');
%temp_table = {'300K', '310K'};
num_of_temp = numel(temp_table);
num_of_traj = 3; % Number of trajectories

% Storage for mean radii, concentration, temperatures, and packed/unpacked labels
mean_radii = [];
mean_conc = [];
temp_values_expanded = [];
packed_labels = [];

% Predefine arrays to store overall means and standard deviations for each temperature and file type (40 columns retained)
mean_radii_final = cell(num_of_temp, 2);
mean_conc_final = cell(num_of_temp, 2);
std_radii_final = cell(num_of_temp, 2);
std_conc_final = cell(num_of_temp, 2);

mean_maxclust_num_alltemp = cell(num_of_temp, 2);
std_maxclust_num_alltemp = cell(num_of_temp, 2);
mean_maxclust_radii_alltemp = cell(num_of_temp, 2);
std_maxclust_radii_alltemp = cell(num_of_temp, 2);

temp_values = zeros(num_of_temp, 1); % Store temperature values

% Initialize cell arrays to store trajectory-specific data
mean_conc_temp_traj = cell(num_of_temp, 2, num_of_traj);

for i = 1:num_of_temp
    % Extract numerical temperature value
    temperature = str2double(regexp(temp_table{i}, '\d+', 'match', 'once'));
    temp_values(i) = temperature;
    
    % Define filenames for compact and unpacked datasets
    radii_files = {['wt_radii_', temp_table{i}, '.xlsx'], ...
                   ['wt_unpack_radii_', temp_table{i}, '.xlsx']};
    conc_files = {['wt_concentrations_', temp_table{i}, '.xlsx'], ...
                  ['wt_unpack_concentrations_', temp_table{i}, '.xlsx']};

    maxclust_num_radii_files = {['wt_maxClust_SizeRadiusVolConcAtoms_', temp_table{i}, '.xlsx'], ...
                   ['wt_unpack_maxClust_SizeRadiusVolConcAtoms_', temp_table{i}, '.xlsx']};
    
    
    for j = 1:2 % Loop over compact (1) and unpacked (2) datasets
   
            % Initialize storage for 40-column means across trajectories
            mean_radii_temp = zeros(40, num_of_traj);
            mean_conc_temp = zeros(40, num_of_traj);
            maxclust_num = zeros(100, num_of_traj);
            maxclust_radii = zeros(100, num_of_traj);
            for traj = 1:num_of_traj
                % Read radii and concentration data from the current sheet (trajectory)
                radii = readmatrix(radii_files{j}, 'Sheet', traj);
                radii_data = radii(501:1000,:);
                conc = readmatrix(conc_files{j}, 'Sheet', traj);
                conc_data = conc(501:1000,:);
                maxclust_num_radii = readmatrix(maxclust_num_radii_files{j}, 'Sheet', traj);
                maxclust_num_radii_data = maxclust_num_radii(501:1000,:);
                % Compute row means (100 rows x 40 columns) and then column-wise mean
                mean_radii_temp(:, traj) = mean(radii_data, 1, 'omitnan');
                mean_conc_temp(:, traj) = mean(conc_data, 1, 'omitnan');

                mean_maxclust_num_temp(:, traj) = round(mean(maxclust_num_radii_data(:,2), 2, 'omitnan'));
                mean_maxclust_radii_temp(:, traj) = mean(maxclust_num_radii_data(:,3), 2, 'omitnan');
                mean_maxclust_conc_temp(:, traj) = mean(maxclust_num_radii_data(:,6), 2, 'omitnan');

                % Store per-trajectory mean concentrations in the cell array
                mean_conc_temp_traj{i, j, traj} = mean_conc_temp(:, traj);

                maxclust_num_temp_traj{i, j, traj}  = round(mean_maxclust_num_temp(:, traj));
                maxclust_radii_temp_traj{i, j, traj}  = mean_maxclust_radii_temp(:, traj);
                maxclust_conc_temp_traj{i, j, traj}  = mean_maxclust_conc_temp(:, traj);

                % Ensure that the data fits within the 101 rows for maxclust_num and maxclust_radii
            % Adjust data size by truncating or padding as necessary
            num_frames = size(maxclust_num_radii_data, 1);
            
            maxclust_num(1:num_frames, traj) = maxclust_num_radii_data(:, 2);
            maxclust_radii(1:num_frames, traj) = maxclust_num_radii_data(:, 3);

            % Calculate the mean of maxclust_num and maxclust_radii across all trajectories
            % mean_trajmaxclust_num = round(mean(maxclust_num, 2, 'omitnan'));  % Column-wise mean for maxclust_num
            % mean_trajmaxclust_radii = round(mean(maxclust_radii, 2, 'omitnan'));  % Column-wise mean for maxclust_radii
     
             % Store the final means in the predefined cell arrays
             %mean_trajmaxclust_num_alltemp{i, j} = mean_trajmaxclust_num;
             %mean_trajmaxclust_radii_alltemp{i, j} = mean_trajmaxclust_radii;

            end

           mean_trajmaxclust_num = round(mean(maxclust_num, 2, 'omitnan'));  % Column-wise mean for maxclust_num
           std_trajmaxclust_num = std(maxclust_num, 0, 2, 'omitnan');  % Column-wise standard deviation for maxclust_num

           mean_trajmaxclust_radii = round(mean(maxclust_radii, 2, 'omitnan'));  % Column-wise mean for maxclust_radii
           std_trajmaxclust_radii = std(maxclust_radii, 0, 2, 'omitnan');  % Column-wise standard deviation for maxclust_radii
           
           mean_trajmaxclust_num_alltemp{i, j} = mean_trajmaxclust_num;
           mean_trajmaxclust_radii_alltemp{i, j} = mean_trajmaxclust_radii;

           std_trajmaxclust_num_alltemp{i, j} = std_trajmaxclust_num;
           std_trajmaxclust_radii_alltemp{i, j} = std_trajmaxclust_radii;


            % Compute mean and standard deviation across all trajectories column-wise
            mean_radii_filtered = mean(mean_radii_temp, 2, 'omitnan');
            std_radii_filtered = std(mean_radii_temp, 0, 2, 'omitnan');
            mean_conc_filtered = mean(mean_conc_temp, 2, 'omitnan');
            std_conc_filtered = std(mean_conc_temp, 0, 2, 'omitnan');
            
           mean_maxclust_num_temp_final = round(mean(mean_maxclust_num_temp, 2, 'omitnan'));
           std_maxclust_num_temp_final = std(mean_maxclust_num_temp, 0, 2, 'omitnan');
           mean_maxclust_radii_temp_final = mean(mean_maxclust_radii_temp, 2, 'omitnan');
           std_maxclust_radii_temp_final = std(mean_maxclust_radii_temp, 0, 2, 'omitnan');
           %disp(mean_maxclust_radii_temp_final);

           mean_maxclust_num_alltemp{i, j} = mean_maxclust_num_temp_final;
           std_maxclust_num_alltemp{i, j} = std_maxclust_num_temp_final;
           mean_maxclust_radii_alltemp{i, j} = mean_maxclust_radii_temp_final;
           std_maxclust_radii_alltemp{i, j} = std_maxclust_radii_temp_final;

            % Apply filtering: and keep only mean_radii ≤ 20 nm (since Å → nm)
            valid_indices = (mean_radii_filtered <= 75);
            mean_radii_final{i, j} = mean_radii_filtered(valid_indices);
            std_radii_final{i, j} = std_radii_filtered(valid_indices);
            mean_conc_final{i, j} = mean_conc_filtered(valid_indices);
            std_conc_final{i, j} = std_conc_filtered(valid_indices);
            
            % Store temperature and labels
            temp_values_expanded = [temp_values_expanded; temperature];
            packed_labels = [packed_labels; j]; % 1 for packed, 2 for unpacked
         
    end
end
%%
% ================================================
% After all your existing loops: now save mean values in a table
% ================================================

% Initialize empty arrays to store the final summary
Summary_Temperature = [];
Summary_Packed_Unpacked = [];
Summary_Mean_MaxClustNum = [];
Summary_Mean_MaxClustConc = [];

for i = 1:num_of_temp
    temperature = temp_values(i);

    for j = 1:2 % 1: packed, 2: unpacked
        % From previously computed means
        mean_maxclust_num = round(mean(maxclust_num_temp_traj{i, j, 1})); % Using traj 1, you could average across trajs if you want
        mean_maxclust_conc = mean(maxclust_conc_temp_traj{i, j, 1});

        % Append the summary
        Summary_Temperature = [Summary_Temperature; temperature];
        if j == 1
            Summary_Packed_Unpacked = [Summary_Packed_Unpacked; "Packed"];
        else
            Summary_Packed_Unpacked = [Summary_Packed_Unpacked; "Unpacked"];
        end
        Summary_Mean_MaxClustNum = [Summary_Mean_MaxClustNum; mean_maxclust_num];
        Summary_Mean_MaxClustConc = [Summary_Mean_MaxClustConc; mean_maxclust_conc*1000];
    end
end

% Create table
SummaryTable = table(Summary_Temperature, Summary_Packed_Unpacked, Summary_Mean_MaxClustNum, Summary_Mean_MaxClustConc);

% Save the table to an Excel file
writetable(SummaryTable, 'wt_MaxCluster_Mean_SizeConc_Summary_alltemp.xlsx');

disp('Summary table saved as wt_MaxCluster_Mean_SizeConc_Summary_alltemp.xlsx');

% Parameters
Name = {'wt', 'wt_unpack'};
temps = {'270K', '280K', '290K', '300K', '310K', '320K', '330K', '340K'};
temp_values = 270:10:340; % numeric temp

num_of_temp = numel(temps);
Names = numel(Name);
numTrajs = 3;
numShells = 9;

% Load Data
mergedConcs = cell(Names, num_of_temp);

for j = 1:Names
    for i = 1:num_of_temp
        concFile = strcat(Name{j}, '_concentrations_', temps{i}, '.xlsx');
        combinedConcs = [];

        for trajIdx = 1:numTrajs
            concs = xlsread(concFile, trajIdx);
            meanConcs = mean(concs(501:end, 1:numShells), 1) * 1000; % Convert to µM
            combinedConcs = [combinedConcs, meanConcs];
        end

        mergedConcs{j, i} = combinedConcs;
    end
end

%% K-means clustering at 270K reference
all_dense_conc = [];
all_dense_temp = [];
all_dilute_conc = [];
all_dilute_temp = [];

for j = 1:Names
    % 270K reference clustering
    refData = mergedConcs{j, 1};
    k = 2;
    [idx_ref, centers_ref] = kmeans(refData', k);
    [sortedCenters, refOrder] = sort(centers_ref);

    for i = 1:num_of_temp
        concData = mergedConcs{j, i};

        if isempty(concData)
            continue;
        end

        % Classification based on reference
        distances = abs(concData' - sortedCenters');
        [~, assignedClusters] = min(distances, [], 2);

        dilutePhase = concData(assignedClusters == 1);
        densePhase = concData(assignedClusters == 2);

        % Apply filters: dilutePhase between 1 and 13 µM
        %dilutePhase = dilutePhase(dilutePhase > 0.8 & dilutePhase < 14);

        % Store all points
        all_dense_conc = [all_dense_conc; densePhase(:)];
        all_dense_temp = [all_dense_temp; temp_values(i) * ones(length(densePhase), 1)];

        all_dilute_conc = [all_dilute_conc; dilutePhase(:)];
        all_dilute_temp = [all_dilute_temp; temp_values(i) * ones(length(dilutePhase), 1)];
    end
end

%% Smooth Fit
% Sort first
[sorted_dense_temp, idx_dense] = sort(all_dense_temp);
sorted_dense_conc = all_dense_conc(idx_dense);

[sorted_dilute_temp, idx_dilute] = sort(all_dilute_temp);
sorted_dilute_conc = all_dilute_conc(idx_dilute);

% Smooth data using LOESS
dense_fit = smooth(sorted_dense_temp, sorted_dense_conc, 0.3, 'loess');
dilute_fit = smooth(sorted_dilute_temp, sorted_dilute_conc, 0.3, 'loess');


% %

colors = lines(4); % just auto colors
for j = 1:Names
    figure('Position', [100, 100, 1000, 350]);

    % Tiled layout
    t = tiledlayout(2, num_of_temp, 'TileSpacing', 'compact', 'Padding', 'tight');

    % Get 270K reference clustering
    k = 2; % 2 clusters: dilute, dense
    refData = mergedConcs{j, 1}; % 270K data
    [idx_ref, centers_ref] = kmeans(refData', k);
    [sortedCenters, refOrder] = sort(centers_ref);

    % Apply reference clustering to all temperatures
    for i = 1:num_of_temp
        concData = mergedConcs{j, i};

        if isempty(concData)
            continue;
        end

        % Histogram plot
        nexttile(i);
        histogram(concData, 12, 'FaceColor', 'b', 'EdgeColor', 'k');
        xlabel('$\mathit{Conc}$ ($\mu$M)', 'FontSize', 10, 'Interpreter', 'latex');
        if i == 1
            ylabel('Frequency', 'Interpreter', 'latex', 'FontSize', 10);
        else
            yticklabels([]);
        end
        title(sprintf('%s', temps{i}), 'Interpreter', 'latex', 'FontSize', 10);
        grid off;
        xlim([0 25]);
        ylim([0 45]);

        % Classify based on distance to 270K centroids
        distances = abs(concData' - sortedCenters');
        [~, assignedClusters] = min(distances, [], 2);

        % Extract dilute and dense phases
        dilutePhasee = concData(assignedClusters == 1);
        densePhase = concData(assignedClusters == 2);

        % Apply filters to dilute phase (only 0 < conc < 15 µM)
        dilutePhase = dilutePhasee(dilutePhasee < 13 & dilutePhasee > 1);
        %densePhase = dilutePhase(dilutePhase > 10);

        % Scatter plot
        nexttile(num_of_temp + i);
        scatter(1:length(dilutePhase), dilutePhase, 20, 'r', 'filled'); hold on;
        scatter((length(dilutePhase)+1):(length(dilutePhase)+length(densePhase)), densePhase, 25, 'b', 'filled');
        hold off;
        xlabel('Shell Index', 'Interpreter', 'latex', 'FontSize', 10);
        if i == 1
            ylabel('$\mathit{Conc}$ ($\mu$M)', 'FontSize', 10, 'Interpreter', 'latex');
        else
            yticklabels([]);
        end
        if i == num_of_temp
            legend('Dilute Phase', 'Dense Phase', 'Location', 'north', 'Box', 'off', 'Interpreter', 'latex', 'FontSize', 10);
        end
        grid off;
        ylim([0 25]);

        % Calculate means and standard errors
        if ~isempty(densePhase)
            dense_means(j, i) = mean(densePhase);
            dense_errors(j, i) = std(densePhase) / sqrt(numel(densePhase));
        else
            dense_means(j, i) = NaN;
            dense_errors(j, i) = NaN;
        end

        if ~isempty(dilutePhase)
            dilute_means(j, i) = mean(dilutePhase);
            dilute_errors(j, i) = std(dilutePhase) / sqrt(numel(dilutePhase));
        else
            dilute_means(j, i) = NaN;
            dilute_errors(j, i) = NaN;
        end
    end

    % Save figure
    print(gcf, strcat(Name{j}, '_ShellConcs_histo_refKmeans2_temp.tif'), '-dtiff', '-r300');
end

% Plot data
figure('Position', [100, 100, 800, 600]); % [left, bottom, width, height]
tiledlayout(2,2, 'TileSpacing', 'tight', 'Padding', 'compact'); % Tighter layout

traj_colors = {[1, 0, 0], [0, 0, 1], [1, 0.5, 0]}; % Red, Blue, Orange

% ---- Plot: avg maxclust per traj over time, grouped by temp ----
nexttile;
hold on;
cmap = jet(num_of_temp); % Colormap for temperatures
symbols = {'o', '^'};     % 'o' for packed, '^' for unpacked

% Define actual temperature values
tempss = [270, 280, 290, 300, 310, 320, 330, 340];  % Update if needed

% Loop over packed (1) and unpacked (2) datasets
for j = 1:2  % 1 = packed, 2 = unpacked
    for i = 1:num_of_temp

        % Compute mean across time for each trajectory
        traj_means = zeros(1, num_of_traj);
        for traj = 1:num_of_traj
            time_series = mean_trajmaxclust_num_alltemp{i, j}(traj, :);
            traj_means(traj) = mean(time_series);
        end

        % Mean and SEM across trajectories
        y_mean = mean(traj_means);
        y_sem  = std(traj_means) / sqrt(num_of_traj);

        % Plot one point per temperature
        errorbar(tempss(i), y_mean, y_sem, symbols{j}, ...
            'Color', cmap(i, :), ...
            'MarkerFaceColor', cmap(i, :), ...
            'CapSize', 5, 'MarkerSize', 6, 'LineWidth', 1.2);
    end
end

xlabel('Temperature (K)', 'Interpreter', 'latex', 'FontSize', 12);

box off;

xlim([265 345]);
ylim([3, 60]);
ylabel('num of proteins in condensate', 'Interpreter', 'latex', 'FontSize', 12);
yticks_values = yticks;  % Get current y-tick values

% Convert tick values to LaTeX format
%xticklabels(arrayfun(@(x) sprintf('$\\mathrm{%g}$', x), xticks_values, 'UniformOutput', false));
yticklabels(arrayfun(@(y) sprintf('$\\mathrm{%g}$', y), yticks_values, 'UniformOutput', false));

% Set LaTeX interpreter
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

grid off;
hold off;



% ---- Plot 2: num of max clust vs max clust Radii ----
nexttile; 
hold on;
for i = 1:num_of_temp
    for j = 1:2 % 1 = packed, 2 = unpacked
        errorbar(mean_maxclust_radii_alltemp{i, j}/10, mean_maxclust_num_alltemp{i, j}, std_maxclust_num_alltemp{i, j}/5, std_maxclust_num_alltemp{i, j}/5, ...
                 std_maxclust_radii_alltemp{i, j}/90, std_maxclust_radii_alltemp{i, j}/90, symbols{j}, 'Color', cmap(i, :), 'MarkerFaceColor', cmap(i, :), 'CapSize', 2, 'MarkerSize', 3.5, 'LineWidth', 0.5);
    end
end
xlim([1.8 4.1]); 
ylim([3, 60]);
% Formatting the plot
xlabel('$\mathit{R_c} $(nm)', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('num of proteins in codensate ', 'Interpreter', 'latex', 'FontSize', 12);
xticks_values = xticks;  % Get current x-tick values
yticks_values = yticks;  % Get current y-tick values

% Convert tick values to LaTeX format
xticklabels(arrayfun(@(x) sprintf('$\\mathrm{%g}$', x), xticks_values, 'UniformOutput', false));
yticklabels(arrayfun(@(y) sprintf('$\\mathrm{%g}$', y), yticks_values, 'UniformOutput', false));

% Set LaTeX interpreter
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

hold on;

% Add dummy plots for black legend symbols
dummy_compact = plot(nan, nan, 'ko', 'LineStyle', '-', 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', 'compact to PS', 'MarkerFaceColor', 'k');
dummy_dispersed = plot(nan, nan, 'k^', 'LineStyle', '-', 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', 'dispersed to PS', 'MarkerFaceColor', 'k');

% Add legend using the dummy plots
legend([dummy_compact, dummy_dispersed], {'compact to PS', 'dispersed to PS'}, ...
       'Location', 'northwest', 'Box', 'off', 'FontSize', 10, 'Interpreter', 'latex');
% legend({'Compact to PS', 'Dispersed to PS'}, 'Location', 'northeast', 'Box','off');

grid off;
hold off;

% ---- Plot 3: Mean shell conc vs. shell radii ----
nexttile; 
hold on;
for i = 1:num_of_temp
    for j = 1:2 % 1 = packed, 2 = unpacked
        errorbar(mean_radii_final{i, j}/10, mean_conc_final{i, j}*1000, std_conc_final{i, j}*1000, std_conc_final{i, j}*1000, ...
                 std_radii_final{i, j}/20, std_radii_final{i, j}/20, symbols{j}, 'Color', cmap(i, :), 'MarkerFaceColor', cmap(i, :), 'CapSize', 2.5, 'MarkerSize', 4, 'LineWidth', 0.5);
    end
end
xlim([0.85 8]); 
 ylim([-0.5, 25]);
% Formatting the plot
xlabel('$\mathit{R_{sh}} $(nm)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('concentration ($\mu$M)', 'Interpreter', 'latex', 'FontSize', 12);
xticks_values = xticks;  % Get current x-tick values
yticks_values = yticks;  % Get current y-tick values

% Convert tick values to LaTeX format
xticklabels(arrayfun(@(x) sprintf('$\\mathrm{%g}$', x), xticks_values, 'UniformOutput', false));
yticklabels(arrayfun(@(y) sprintf('$\\mathrm{%g}$', y), yticks_values, 'UniformOutput', false));

% Set LaTeX interpreter
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

% Generate formatted tick labels in LaTeX math mode (non-italic)
tick_values = linspace(min(temp_values_expanded), max(temp_values_expanded), 8);
tick_labels = arrayfun(@(t) sprintf('$\\mathrm{%g}$', t), tick_values, 'UniformOutput', false);

% Add colorbar and adjust its position
cb = colorbar('Ticks', linspace(0, 1, 8), ... % Tick positions
              'Box', 'off', ... % Remove border around the colorbar
              'FontSize', 9);  % Adjust font size

% Generate formatted tick labels in LaTeX math mode (non-italic)
tick_values = linspace(min(temp_values_expanded), max(temp_values_expanded), 8);
tick_labels = arrayfun(@(t) sprintf('$\\mathrm{%g}$', t), tick_values, 'UniformOutput', false);

% Apply LaTeX interpreter for the colorbar tick labels
cb.TickLabels = tick_labels;
set(cb, 'TickLabelInterpreter', 'latex');

%clim([min(temp_values_expanded), max(temp_values_expanded)]);

cb.Position = cb.Position + [-0.08, 0.18, -0.01, -0.19]; % Shift colorbar left to be closer to the plot
colormap(cmap); % Apply the colormap

% Optional: Customize the colorbar label
ylabel(cb, 'temperature (K)', 'FontSize', 12, 'Interpreter', 'latex');

grid off;
hold off;

% ---- Plot 4: Concentration vs. Temperature ----
nexttile;
hold on;
cmap = parula(num_of_traj); % Assign different colors for each trajectory
symbols = {'o', '^'}; % 'o' for packed, 's' for unpacked
legend_entries = {}; % Store legend labels

% Loop over packed (1) and unpacked (2) datasets
for j = 1:2  % 1 = packed, 2 = unpacked
    for traj = 1:num_of_traj
        mean_conc_per_temp = zeros(num_of_temp, 1);
        std_conc_per_temp = zeros(num_of_temp, 1);

        for i = 1:num_of_temp
            if isempty(mean_radii_final{i, j})
                continue; % Skip if no valid data
            end

            % Extract valid concentration values where radii ≤ 160 Å
            valid_indices = (mean_radii_final{i, j} <= 30);

            % Ensure that concentration data matches valid indices
            conc_values = mean_conc_temp_traj{i, j, traj};
            if numel(conc_values) >= numel(valid_indices)
                conc_values = conc_values(valid_indices);
            end
            conc_valuess = conc_values(conc_values > 0);
            % Compute mean and standard deviation across valid data points
            mean_conc_per_temp(i) = mean(conc_valuess, 'omitnan');
            std_conc_per_temp(i) = std(conc_valuess, 0, 'omitnan');
        end

        % Plot mean concentration vs. temperature with error bars
        % Define fixed colors for the three trajectories


         errorbar(temp_values, mean_conc_per_temp * 1000, std_conc_per_temp * 500, symbols{j}, ...
         'Color', traj_colors{traj}, 'MarkerFaceColor', 'none', 'MarkerSize', 4, ...
         'LineWidth', 0.5, 'CapSize', 2.5);


        % Store legend labels
        if j == 1
            legend_entries{end+1} = sprintf('compact to PS %d C', traj);
        else
            legend_entries{end+1} = sprintf('dispersed to PS %d C', traj);
        end
    end
end
xlim([265 345]); 
ylim([5 25]);
% Formatting the plot
xlabel('temperature (K)', 'Interpreter', 'latex', 'FontSize', 12);
%ylabel('$\mathit{C}$ ($\mu$M)', 'Interpreter', 'latex', 'FontSize', 12);
xticks_values = xticks;  % Get current x-tick values
yticks_values = yticks;  % Get current y-tick values

% Convert tick values to LaTeX format
xticklabels(arrayfun(@(x) sprintf('$\\mathrm{%g}$', x), xticks_values, 'UniformOutput', false));
yticklabels(arrayfun(@(y) sprintf('$\\mathrm{%g}$', y), yticks_values, 'UniformOutput', false));

% Set LaTeX interpreter
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

%title('Concentration vs. Temperature (Filtered for Radii ≤ 160 Å)');
% legend(legend_entries, 'Location', 'southwest', 'Box','off', 'FontSize', 8, 'Interpreter', 'latex');
% Define the legend manually inside the plot
legend_handle = legend(legend_entries, 'Box', 'off', 'FontSize', 9, 'Interpreter', 'latex', 'NumColumns', 1);

% Adjust the position of the legend manually [left, bottom, width, height]
legend_handle.Position = [0.45, 0.1, 0.35, 0.17]; % Adjust these values as needed

grid on;
hold on;

% Save the figure as a high-quality TIFF image (300 DPI)
 % %print(gcf, 'wt_conc_maxclust_vs_radii_and_temp_lst5mics.tif', '-dtiff', '-r300');

xlabel('temperature (K)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('concentration ($\mu$M)', 'Interpreter', 'latex', 'FontSize', 14);
%legend({'wt Dilute', 'wt Dense', 'wt-unpack Dilute', 'wt-unpack Dense'}, ...
%    'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12);
%scatter(all_dilute_temp, all_dilute_conc, 35, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
%scatter(all_dense_temp, all_dense_conc, 35, 'r', 'filled', 'MarkerFaceAlpha', 0.5);

plot(sorted_dilute_temp, dilute_fit*2.3, 'r--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
plot(sorted_dense_temp, dense_fit, 'r--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

%xlabel('Temperature (K)', 'Interpreter', 'latex', 'FontSize', 16);
%ylabel('Concentration ($\mu$M)', 'Interpreter', 'latex', 'FontSize', 16);
%legend({'wt Dilute', 'wt Dense', 'wt-unpack Dilute', 'wt-unpack Dense'}, ...
%    'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12);
%title('Phase Separation from Concentration Profiles', 'Interpreter', 'latex', 'FontSize', 18);

%grid on;
%ylim([5 25]);
%xlim([265 345]);
%box on;

hold off;
%%
% Save figure (optional)
  print(gcf, 'wt_conc_maxclust_vs_radii_and_temp_DenseDilute_SmoothFit2.tif', '-dtiff', '-r300');
















clear;

% Define file names for the cases
Names = {'wt', 'wt_unpack'};
figname = 'wt';
% Temperatures to analyze
temperatures = 270:10:340;
num_of_temp = numel(temperatures);

% Colormap for temperatures
cmap = jet(num_of_temp);

% Initialize storage for MSD and cluster data
msd_data = struct();
cluster_data = struct();

% Process each case
for k = 1:numel(Names)
    Name = Names{k};
    msd_data.(Name).mean = [];
    msd_data.(Name).std = [];
    cluster_data.(Name) = [];

    for idx = 1:num_of_temp
        i = temperatures(idx);
        temp_str = int2str(i);

        % % --- MSD Data ---
        % Construct the MSD file name
        fileName = strcat('Aggregated_MSD_', Name, '_maxClust_Fixlag_', temp_str, 'K.csv');

        % Read the CSV file excluding the header row
        try
            data = readmatrix(fileName, 'NumHeaderLines', 1);
        catch
            warning(['MSD File not found: ', fileName]);
            continue;
        end

        % Extract MSD components
        time = data(501:end, 1);
        msd_x = data(501:end, 2);
        msd_y = data(501:end, 3);
        msd_z = data(501:end, 4);
        combined_msd = sqrt(msd_x + msd_y + msd_z);

        % Mean and STD of combined MSD
        mean_combined = mean(combined_msd, 'omitnan');
        std_combined = std(combined_msd, 'omitnan');

        % Store MSD results
        msd_data.(Name).mean = [msd_data.(Name).mean; i, mean_combined];
        msd_data.(Name).std = [msd_data.(Name).std; i, std_combined];

        % % --- MaxCluster Data ---
        cluster_filename = strcat(Name, '_maxClust_SizeRadiusVolConcAtoms_', temp_str, 'K.xlsx');
        %wt_unpack_maxClust_SizeRadiusVolConcAtoms_320K.xlsx

        cluster_means = [];
        try
            [~, sheets] = xlsfinfo(cluster_filename);
            for s = 1:length(sheets)
                sheet_data = readmatrix(cluster_filename, 'Sheet', sheets{s});
                if size(sheet_data, 2) >= 2
                    cluster_values = round(sheet_data(:, 2)); % Column 2 = MaxCluster
                    cluster_means = [cluster_means; mean(cluster_values, 'omitnan')];
                end
            end
        catch
            warning(['Cluster file not found or unreadable: ', cluster_filename]);
        end

        % Average cluster size across sheets (trajectories)
        mean_cluster = round(mean(cluster_means, 'omitnan'));
        cluster_data.(Name) = [cluster_data.(Name); i, mean_cluster];
    end
end

% %

% % Plotting MSD vs time in separate plot

% Define y-axis limits for the plot
common_ylim = [0, 7];

figure('Position', [100, 100, 400, 300]);
t = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact'); % Create single plot
% Iterate over cases
for k = 1:numel(Names)
    Name = Names{k};

    % Extract data
    case_data = msd_data.(Name);
    if isempty(case_data.mean)
        continue; % Skip if no data exists for this case
    end

    temperatures = case_data.mean(:, 1);
    means_combined = smooth(case_data.mean(:, 2));
    stds_combined = case_data.std(:, 2);

    % Plot filled region for standard deviation
    fill([temperatures; flipud(temperatures)], ...
         [means_combined - stds_combined; flipud(means_combined + stds_combined)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    % Plot smoothed mean line
    hold on;
    plot(temperatures, means_combined, '-o', 'DisplayName', Name, 'LineWidth', 1.5);

end
xlim([269 341]);
% Customize the plot
xlabel('temperature (K)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\mathrm{MSD} \ (\mathrm{nm}^2)$', 'Interpreter', 'latex', 'FontSize', 12);
ylim(common_ylim);
legend('', 'compact to PS', '', 'dispersed to PS', 'Interpreter', 'latex', 'FontSize', 12, 'box', 'off', 'location', 'best');
grid off;
box off;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);

figName1 = strcat(figname, '_maxClust_combined_MSD_vs_temp_lag.tif');
exportgraphics(gcf, figName1, 'Resolution', 300);


%% --- Plot: Mean MSD vs Mean MaxCluster with colormap for temperature ---

% Define unique marker symbols for each case
marker_styles = {'o', '^'};

% --- Plot: Mean MSD vs Mean MaxCluster with colormap for temperature ---
figure('Position', [100, 100, 350, 275]);
t = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:numel(Names)
    Name = Names{k};
    marker = marker_styles{k};

    % Extract MSD and cluster data
    msd_vals = msd_data.(Name).mean;
    msd_stds = msd_data.(Name).std(:, 2);
    cluster_vals = cluster_data.(Name);

    if isempty(msd_vals) || isempty(cluster_vals)
        continue;
    end

    % Match temperatures
    temps_common = intersect(msd_vals(:,1), cluster_vals(:,1));
    [~, idx_msd] = ismember(temps_common, msd_vals(:,1));
    [~, idx_cluster] = ismember(temps_common, cluster_vals(:,1));

    mean_msd = msd_vals(idx_msd, 2);
    std_msd = msd_stds(idx_msd);

    mean_cluster = round(cluster_vals(idx_cluster, 2));

    % Calculate STD for cluster values across sheets
    std_cluster = zeros(size(mean_cluster));
    for i = 1:length(temps_common)
        temp_str = int2str(temps_common(i));
        cluster_filename = strcat(Name, '_maxClust_SizeRadiusVolConcAtoms_', temp_str, 'K.xlsx');
        cluster_values_all = [];

        [~, sheets] = xlsfinfo(cluster_filename);
        for s = 1:length(sheets)
            sheet_data = readmatrix(cluster_filename, 'Sheet', sheets{s});
            if size(sheet_data, 2) >= 2
                cluster_values = sheet_data(:, 2); % Column 2 = MaxCluster
                cluster_values_all = [cluster_values_all, cluster_values];
            end
        end
        std_cluster(i) = std(mean(cluster_values_all, 1), 'omitnan');
    end

    % Colormap mapping
    cmap_idx = round(rescale(temps_common, 1, num_of_temp));
    cmap_idx = min(max(cmap_idx, 1), num_of_temp);

    % Plot with error bars
    for i = 1:length(mean_msd)
        errorbar(mean_cluster(i), mean_msd(i), std_msd(i)/2, std_msd(i)/2, ...
                 std_cluster(i), std_cluster(i), marker, ...
                 'Color', cmap(cmap_idx(i), :), ...
                 'MarkerFaceColor', cmap(cmap_idx(i), :), ...
                 'MarkerEdgeColor', 'none', ...
                 'CapSize', 3, 'LineWidth', 1.1, ...
                 'DisplayName', Name); 
        hold on;
    end
end
ylim([-1 8])
xlabel('num of proteins in condensate', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$\mathrm{MSD} \ (\mathrm{nm}^2)$', 'Interpreter', 'latex', 'FontSize', 13);
%title('MSD vs MaxCluster Size with STD Error Bars', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 13);
%legend('Interpreter', 'latex', 'Location', 'best');
grid off;
box off;

colormap(cmap);
% c = colorbar;
% c.Ticks = linspace(0, 1, 8);
% c.TickLabels = arrayfun(@(x) sprintf('%dK', x), temperatures, 'UniformOutput', false);
% c.Label.String = 'Temperature';
% c.Label.Interpreter = 'latex';
% %

figName2 = strcat(figname, '_MSD_vs_MaxCluster_colored_by_temp.tif');
exportgraphics(gcf, figName2, 'Resolution', 300);

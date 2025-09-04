clear;

% File names and labels
fileNames = {'wt', 'mt', 'wt_unpack', 'mt_unpack'};
systNames = {'compact to PS (WT)', 'compact to PS (MT)', 'dispersed to PS (WT)', 'dispersed to PS (MT)'};
%file_tag = {'A','B','C','D'};

% Define all temperatures
temps = 270:10:340;

% Prepare figure
numTemps = numel(temps);
numSystems = numel(fileNames);
figure('Position', [100, 100, 200*numSystems, 150*numTemps]);
tiledLayout = tiledlayout(numTemps, numSystems, 'TileSpacing', 'tight', 'Padding', 'compact');

for t = 1:numTemps
    temp = int2str(temps(t));

    for k = 1:numSystems
        % Combine data from all 3 trajectories
        combinedData = [];

        for j = 1:3
            traj = int2str(j);
            %Name = 'alpha0.925_60chain';
            fileName = strcat(fileNames{k}, '_', traj, 'C_maxclustSeleChainIDs_1k_', temp, 'K.csv');

            if ~isfile(fileName)
                disp(['File not found: ', fileName]);
                continue;
            end

            data = readmatrix(fileName, 'NumHeaderLines', 1);
            frameIndices = data(1:1000, 1);
            logicalMatrix = data(1:1000, 2:end);
            combinedData = cat(3, combinedData, logicalMatrix);
        end

        if isempty(combinedData)
            continue;
        end

        finalData = zeros(size(combinedData, 1), size(combinedData, 2));
        for i = 1:size(finalData, 2)
            for j = 1:size(finalData, 1)
                activeCount = sum(combinedData(j, i, :) == 1);
                if activeCount >= 2
                    finalData(j, i) = 1;
                end
            end
        end

        % Start plotting
        nexttile;
        imagesc(finalData');
        colormap([1 1 1; 0 0.2 1]);

        ax = gca;
        ax.TickLabelInterpreter = 'latex';   % <-- Enables LaTeX for ticks
        xticks(ax, 0:100:1000);

        % Read header
        fileID = fopen(fileName, 'r');
        headerLine = fgetl(fileID);
        fclose(fileID);
        headers = strsplit(headerLine, ',');
        chainIDs = headers(2:end);

        yticks(1:size(finalData, 2));
        customYLabels = repmat({''}, 1, length(chainIDs));
        customYLabels(1:2:end) = chainIDs(1:2:end);

        ax.XAxis.FontSize = 10;
        ax.YAxis.FontSize = 6;
        ytickangle(90);

        if k == 1
            yticklabels(ax, customYLabels);
            ylabel('peptides index', 'Interpreter', 'latex', 'FontSize', 10);

            % Label
         text(-0.18, 0.5, [temp 'K'], ...
        'FontSize', 8, 'FontWeight', 'bold', ...
        'Color', 'black', 'BackgroundColor', 'none', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'Interpreter', 'latex', 'Units', 'normalized', ...
        'Rotation', 90);
        else
            yticklabels([]);
        end

        if t == numTemps
            xlabel('time ($\mu$s)', 'Interpreter', 'latex', 'FontSize', 10);
            xticklabels({'0', '', '2', '', '4', '', '6', '', '8', '', '10'});
        else
            xticklabels([]);
        end
       if t == 1
         text(0.5, 1.02, ['\textbf{', systNames{k}, '}'], ...
        'FontSize', 10, ...
        'Color', 'black', 'Interpreter', 'latex', ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom');
       end
        

        % Right Y-axis plot
        hold on;
        totalChains = sum(finalData, 2);
        smoothTotalChains = round(smooth(totalChains, 25));
        yyaxis right;
        plot(frameIndices, smoothTotalChains, 'r', 'LineWidth', 1.0);

        % Stats from last 500 frames
        clusterCount = sum(finalData > 0, 2);
        last500_chains = round(smoothTotalChains(end-499:end));
        last500_clusters = clusterCount(end-499:end);
        avg_maxclust_last500 = round(mean(last500_chains(last500_chains > 0)));
        avg_totalclusters_last500 = round(mean(last500_clusters));

        % stat_text = {
        %     sprintf('avg max clust: %d', avg_maxclust_last500), ...
        %     sprintf('avg total clust: %d', avg_totalclusters_last500)
        % };
        % text(0.98, 0.90, stat_text, ...
        %     'FontSize', 9, 'Color', 'k', ...
        %     'BackgroundColor', 'white', 'EdgeColor', 'none', ...
        %     'Interpreter', 'latex', 'FontWeight', 'bold', ...
        %     'Units', 'normalized', 'HorizontalAlignment', 'right');

        ylim([0 60]);
        if k == numSystems
            ylabel('total peptides', 'Interpreter', 'latex', 'FontSize', 8);
        else
            ax.YAxis(2).TickLabels = [];
        end
        ax.YAxis(2).Color = 'k';
        hold off;
    end
end
% %
% Save the figure
exportgraphics(gcf, 'ChainEvolu_vs_time_Heatmaps_AllTemps.tif', 'Resolution', 300);

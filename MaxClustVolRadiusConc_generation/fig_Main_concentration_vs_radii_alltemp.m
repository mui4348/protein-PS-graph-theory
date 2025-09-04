% Temperature range
temps = 270:10:340;
nTemps = numel(temps);

% Base file name patterns
systems = {'wt', 'mt'};
file_types = {'radii', 'concentrations'};
prefixes = {'', 'unpack_'};
label_prefixes = {'compact', 'dispersed'};

% Define color maps
compact_colors = [0 0 1; 0.2 0.6 1];     % Blue tones for compact
unpack_colors = [1 0 0; 1 0.6 0.2];      % Red/Orange tones for dispersed

% Create tiled layout
figure('Position', [100, 100, 800, 700]);
t = tiledlayout(3, 3, 'TileSpacing', 'tight', 'Padding', 'tight');  % 3x3 grid

% Loop over temperatures
for tIdx = 1:nTemps
    T = temps(tIdx);
    Tstr = sprintf('%dK', T);
    nexttile;
    hold on;
    legendEntries = {};
    
    % Loop over both systems (wt, mt)
    for systemIdx = 1:2
        system_name = systems{systemIdx};
        
        % Loop over compact and unpacked
        for trajType = 1:2
            if trajType == 1
    color = compact_colors(systemIdx,:);
           else
    color = unpack_colors(systemIdx,:);
            end

            label_prefix = label_prefixes{trajType};
            file_prefix = prefixes{trajType};
            
            % Build filenames
            radii_file = sprintf('%s_%sradii_%s.xlsx', system_name, file_prefix, Tstr);
            conc_file  = sprintf('%s_%sconcentrations_%s.xlsx', system_name, file_prefix, Tstr);
            
            if ~isfile(radii_file) || ~isfile(conc_file)
                warning('Missing files for %s at %s. Skipping.', system_name, Tstr);
                continue;
            end
            
            % Read sheet names
            [~, sheets_radii] = xlsfinfo(radii_file);
            [~, sheets_conc] = xlsfinfo(conc_file);
            
            % Check sheet counts
            if length(sheets_radii) ~= length(sheets_conc)
                warning('Sheet count mismatch in %s and %s', radii_file, conc_file);
                continue;
            end
            
            % Plot all sheets
            for sheetIdx = 1:length(sheets_radii)
                radii_data = readmatrix(radii_file, 'Sheet', sheets_radii{sheetIdx});
                conc_data  = readmatrix(conc_file, 'Sheet', sheets_conc{sheetIdx});
                
                if isempty(radii_data) || isempty(conc_data)
                    continue;
                end
                
                % Compute mean and std
                mean_radii = mean(radii_data(501:end,:), 'omitnan') / 10;
                std_radii = std(radii_data(501:end,:), 0, 'omitnan') / 20;
                mean_conc = mean(conc_data(501:end,:), 'omitnan') * 1000;
                std_conc = std(conc_data(501:end,:), 0, 'omitnan') * 500;
                
                % Plot
                errorbar(mean_radii(1:10), mean_conc(1:10), ...
                    std_conc(1:10), std_conc(1:10), ...
                    std_radii(1:10), std_radii(1:10), ...
                    '-', 'Color', color, 'LineWidth', 1, 'CapSize', 3);
                
                %legendEntries{end+1} = sprintf('%s - %s (%s)', system_name, label_prefix, Tstr);
            end
        end
    end
    
    % Axes formatting
    xlim([0.85, 8]);
    ylim([-0.5, 25]);
    if tIdx == 6 ||  tIdx == 7 || tIdx == 8
    xlabel('$\mathrm{R}_{\mathrm{sh}}$ (nm)', 'Interpreter', 'latex', 'FontSize', 12);
    end
    if tIdx == 1 ||  tIdx == 4 || tIdx == 7
    ylabel('$\mathrm{concentration}$ ($\mu$M)', 'Interpreter', 'latex', 'FontSize', 12);
    end
    % title(sprintf('%d K', T), 'Interpreter', 'latex', 'FontSize', 12);
    text(7.2, 24, sprintf('%d K', T), ...
    'Units', 'data', 'FontSize', 12, 'Interpreter', 'latex', ...
    'HorizontalAlignment', 'right');

    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
    grid off;
    
    % if tIdx == 1
    %     legend(legendEntries, 'Location', 'best', 'FontSize', 7, ...
    %         'Interpreter', 'latex', 'Box', 'off');
    % end
    hold off;
end

legendEntries = {' 1C compact to PS (WT)','2C','3C','1C dispersed to PS (WT)','2C','3C','1C compact to PS (MT)','2C','3C','1C dispersed to PS (MT)','2C','3C'};
% %
% Store legend info for later
if tIdx == nTemps
    legend_handle = legend(legendEntries, ...
        'Interpreter', 'latex', 'FontSize', 10, ...
        'Box', 'off');
    
    % Manually position the legend near the lower-right of the last column
    legend_position = [0.75, 0.05, 0.1, 0.26];  % [left, bottom, width, height]
    set(legend_handle, 'Units', 'normalized', 'Position', legend_position);
end



% Save the plot as a high-resolution TIFF file
print(gcf, 'Radii_vs_Concentration_AllTemps_Subplots.tif', '-dtiff', '-r300');

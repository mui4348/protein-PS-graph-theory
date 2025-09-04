
temps = 270:10:340;
num_temps = length(temps);

% Prepare storage
all_degrees = [];
all_edgeweights = [];
all_node_strengths = [];
all_datasets = cell(num_temps, 4); % [Temp × Systems]

% Loop to collect all data first
for t = 1:num_temps
    temp = int2str(temps(t));

    % Read degree tables
    filenames_deg = {
        ['wt_degList_chain60_3k_end5mics_', temp, 'K.xlsx'];
        ['mt_degList_chain60_3k_end5mics_', temp, 'K.xlsx'];
        ['wt_unpack_degList_chain60_3k_end5mics_', temp, 'K.xlsx'];
        ['mt_unpack_degList_chain60_3k_end5mics_', temp, 'K.xlsx'];
    };

    filenames_edge = {
        ['wt_edgeList_chain60_3k_end5mics_', temp, 'K.xlsx'];
        ['mt_edgeList_chain60_3k_end5mics_', temp, 'K.xlsx'];
        ['wt_unpack_edgeList_chain60_3k_end5mics_', temp, 'K.xlsx'];
        ['mt_unpack_edgeList_chain60_3k_end5mics_', temp, 'K.xlsx'];
    };

    for i = 1:4
        deg_tbl = readtable(filenames_deg{i});
        edge_tbl = readtable(filenames_edge{i});
        deg_tbl.Properties.VariableNames = {'Nodes', 'degree'};
        all_datasets{t, i} = {deg_tbl, edge_tbl};

        all_degrees = [all_degrees; deg_tbl.degree];
        all_edgeweights = [all_edgeweights; edge_tbl.Weight];

        G_temp = graph(edge_tbl.EndNodes_1, edge_tbl.EndNodes_2);
        node_strengths_temp = zeros(numnodes(G_temp),1);
        for n = 1:numnodes(G_temp)
            connectedEdges = find(G_temp.Edges.EndNodes(:,1) == n | G_temp.Edges.EndNodes(:,2) == n);
            node_strengths_temp(n) = sum(edge_tbl.Weight(connectedEdges));
        end
        all_node_strengths = [all_node_strengths; node_strengths_temp];
    end
end
%%

% Compute global max high-weight degree
% highweight_threshold = 0.5 * max_edgeweight;  % Adjust threshold as needed
% all_highweight_degrees = [];
% 
% for t = 1:num_temps
%     for s = 1:4
%         edge_tbl = all_datasets{t, s}{2};
%         G = graph(edge_tbl.EndNodes_1, edge_tbl.EndNodes_2);
%         n = numnodes(G);
% 
%         node_highweight_degree = zeros(n,1);
%         for e = 1:numedges(G)
%             if edge_tbl.Weight(e) >= highweight_threshold
%                 n1 = G.Edges.EndNodes(e,1);
%                 n2 = G.Edges.EndNodes(e,2);
%                 node_highweight_degree(n1) = node_highweight_degree(n1) + 1;
%                 node_highweight_degree(n2) = node_highweight_degree(n2) + 1;
%             end
%         end
%         all_highweight_degrees = [all_highweight_degrees; node_highweight_degree];
%     end
% end
% 
% max_highweight_degree = max(all_highweight_degrees);
% 
% 
% % Global normalization % Find GLOBAL max values
% max_degree = max(all_degrees); % <-- for node degree
% max_edgeweight = max(all_edgeweights); % <-- for edge weights
% max_node_strength = max(all_node_strengths); % <-- for node strength

% === Global normalization (based on high-weight edges) ===

% Threshold for what counts as "high-weight" edge
highweight_threshold = 0.5 * max(all_edgeweights);  % Adjust threshold as needed

% Compute global maximum number of high-weight edges per node
max_edgeweight = max(all_edgeweights);  % Still used for edge coloring
max_node_strength = max(all_node_strengths);  % Still used for color

all_highweight_degrees = [];

for t = 1:num_temps
    for s = 1:4
        edge_tbl = all_datasets{t, s}{2};
        G = graph(edge_tbl.EndNodes_1, edge_tbl.EndNodes_2);
        n = numnodes(G);

        node_highweight_degree = zeros(n,1);
        for e = 1:height(edge_tbl)
            if edge_tbl.Weight(e) >= highweight_threshold
                n1 = edge_tbl.EndNodes_1(e);
                n2 = edge_tbl.EndNodes_2(e);
                node_highweight_degree(n1) = node_highweight_degree(n1) + 1;
                node_highweight_degree(n2) = node_highweight_degree(n2) + 1;
            end
        end

        all_highweight_degrees = [all_highweight_degrees; node_highweight_degree];
    end
end

max_highweight_degree = max(all_highweight_degrees);  % Used to normalize node size



% Setup colormaps
cmap_nodes = turbo(256); % Choose colormap, e.g., parula, jet, hot, etc.
cmap_edges = cool(256);

% Edge width limits
min_edge_width = 0.8;
max_edge_width = 10;

% Edge transparency and threshold
min_alpha = 0.1; % ↓ lower minimum transparency
max_alpha = 0.3; % ↓ maximum transparency (not fully opaque)
hide_threshold = 0.20; % optional: lower it a bit too if you want more edges visible
%
% Create figure
figure('Position', [50, 50, 500, 200*4]);  % Width for 4 cols, height for 8 rows
tl = tiledlayout(8, 4, ...
    'TileSpacing', 'none', ...
    'Padding', 'compact');

system_labels = {'compact (WT)', 'compact (MT)', 'dispersed (WT)', 'dispersed (MT)'};
temps = [270, 280, 290, 300, 310, 320, 330, 340]; % Ensure this matches your data
num_temps = numel(temps);

text_handles = gobjects(0); % For storing text label handles (if needed)

for t = 1:num_temps
    temp_val = temps(t);
    for s = 1:4
        tile_idx = (t - 1) * 4 + s;  % Row-major order (8 rows × 4 columns)
        nexttile(tile_idx);
        hold on

        data = all_datasets{t, s}; % Make sure all_datasets{temp_idx, system_idx} is valid
        degree_tbl = data{1};
        edge_tbl = data{2};
        G = graph(edge_tbl.EndNodes_1, edge_tbl.EndNodes_2);

        p = plot(G, 'Layout', 'force3', 'Iterations', 300, ...
                 'WeightEffect', 'direct', 'UseGravity', false);
        X = p.XData; Y = p.YData; Z = p.ZData; delete(p);

        % Apply spatial spread factor
        spread_factor = 1.0;
        X = spread_factor * X; Y = spread_factor * Y; Z = spread_factor * Z;

        % Compute node high-weight degrees
        node_highweight_degree = zeros(numnodes(G),1);
        for e = 1:height(edge_tbl)
            if edge_tbl.Weight(e) >= highweight_threshold
                n1 = edge_tbl.EndNodes_1(e);
                n2 = edge_tbl.EndNodes_2(e);
                node_highweight_degree(n1) = node_highweight_degree(n1) + 1;
                node_highweight_degree(n2) = node_highweight_degree(n2) + 1;
            end
        end

        gamma = 2.0;
        normalized_degrees = node_highweight_degree / max_highweight_degree;
        normalized_degrees_gamma = normalized_degrees .^ gamma;

        % Node strengths
        node_strengths = zeros(numnodes(G),1);
        for n = 1:numnodes(G)
            connectedEdges = find(G.Edges.EndNodes(:,1) == n | G.Edges.EndNodes(:,2) == n);
            node_strengths(n) = sum(edge_tbl.Weight(connectedEdges));
        end
        normalized_strengths = node_strengths / max_node_strength;
        color_idx_nodes = max(1, round(normalized_strengths * 255) + 1);
        node_colors = cmap_nodes(color_idx_nodes, :);

        % Node sizes
        pos = [X(:), Y(:), Z(:)];
        D = pdist2(pos, pos); D(D==0) = inf;
        min_dist = min(D(:));
        overlap_fraction = 0.3;
        max_allowed_radius = overlap_fraction * min_dist;
        min_node_radius = 0.5 * max_allowed_radius;
        max_node_radius = 2.5 * max_allowed_radius;
        node_radii = min_node_radius + (max_node_radius - min_node_radius) * normalized_degrees_gamma;

        % Edge rendering
        normalized_edgeweights = edge_tbl.Weight / max_edgeweight;
        color_idx_edges = max(1, round(normalized_edgeweights * 255) + 1);
        edge_colors = cmap_edges(color_idx_edges, :);
        edge_widths = min_edge_width + (max_edge_width - min_edge_width) * normalized_edgeweights;
        edge_alphas = min_alpha + (max_alpha - min_alpha) * normalized_edgeweights.^0.7;

        for e = 1:numedges(G)
            n1 = G.Edges.EndNodes(e,1);
            n2 = G.Edges.EndNodes(e,2);
            if normalized_edgeweights(e) >= hide_threshold
                patch([X(n1) X(n2)], [Y(n1) Y(n2)], [Z(n1) Z(n2)], ...
                      'k', 'EdgeColor', edge_colors(e,:), ...
                      'EdgeAlpha', edge_alphas(e), ...
                      'FaceColor', 'none', ...
                      'LineWidth', edge_widths(e));
            end
        end

        % Draw node spheres
        [sx, sy, sz] = sphere(20);
        for n = 1:numnodes(G)
            surf( ...
                node_radii(n)*sx + X(n), ...
                node_radii(n)*sy + Y(n), ...
                node_radii(n)*sz + Z(n), ...
                'FaceColor', node_colors(n,:), ...
                'EdgeColor', 'none', ...
                'FaceLighting', 'gouraud', ...
                'AmbientStrength', 0.2, ...
                'SpecularStrength', 1.0, ...
                'DiffuseStrength', 0.8, ...
                'SpecularExponent', 30);
        end

        axis equal off

        % Add column titles only at top row
        if t == 1
            title(sprintf('\\textbf{%s}', system_labels{s}), ...
                  'Interpreter', 'latex', ...
                  'FontSize', 10, ...
                  'FontWeight', 'bold');
        end

        % Add temperature labels only at first column
        if s == 1
            th = text(-5.5, 0, 0, sprintf('%dK', temp_val), ...
                'FontSize', 12, ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Rotation', 90, ...
                'Interpreter', 'latex');
            text_handles(end+1) = th;
        end
    end
end

% Billboard label updater
hax = gca;
addlistener(hax, 'CameraPosition', 'PostSet', @(src,evt) billboard_labels(text_handles, hax));
%%
% Create a dedicated axis for the colorbar on the far right
cb_ax = axes('Position', [0.92, 0.75, 0.01, 0.2]); % [left, bottom, width, height]
colormap(cb_ax, cmap_nodes);
clim(cb_ax, [0 1]);

% Create the colorbar
cb = colorbar(cb_ax, 'Location', 'eastoutside');
cb.Ticks = [0 0.25 0.5 0.75 1];
cb.TickLabels = {'0','0.25','0.5','0.75','1'};
cb.Label.String = 'node strength';
cb.Label.FontSize = 10;
cb.Label.Interpreter = 'latex';
cb.Box = 'off';
axis(cb_ax, 'off'); % Hide axis around colorbar

%%
% Optional: save the figure
 print(gcf, 'network_graph_allTemps.tif', '-dtiff', '-r300');

% Billboard label function
function billboard_labels(text_handles, hax)
    campos = hax.CameraPosition;
    camtarget = hax.CameraTarget;
    viewvec = campos - camtarget;
    az = atan2d(viewvec(2), viewvec(1));
    for k = 1:numel(text_handles)
        if isvalid(text_handles(k))
            set(text_handles(k), 'Rotation', -az);
        end
    end
end

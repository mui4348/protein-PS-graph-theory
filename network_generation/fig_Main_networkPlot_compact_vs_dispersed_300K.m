temps = 300;
num_temps = 1;

% Prepare storage
all_degrees = [];
all_edgeweights = [];
all_node_strengths = [];
all_datasets = cell(num_temps, 4); % [1 temp × 4 systems]

% Collect data for 300 K
temp = int2str(temps);

% File lists for the 4 systems at 300K
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
    all_datasets{1, i} = {deg_tbl, edge_tbl};

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

% === Global normalization ===
highweight_threshold = 0.5 * max(all_edgeweights);
max_edgeweight = max(all_edgeweights);
max_node_strength = max(all_node_strengths);

all_highweight_degrees = [];

for s = 1:4
    edge_tbl = all_datasets{1, s}{2};
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

max_highweight_degree = max(all_highweight_degrees);

% Colormaps
cmap_nodes = turbo(256);
cmap_edges = cool(256);
min_edge_width = 0.1;
max_edge_width = 1;
min_alpha = 0.1;
max_alpha = 0.3;
hide_threshold = 0.20;
%%
% === Plot ===
figure('Position', [100, 100, 600, 600]);
tl = tiledlayout(2, 2, 'TileSpacing', 'tight', 'Padding', 'compact');
%tl.TileSpacing = 'tight';  % or 'compact', 'none'

system_labels = {'compact to PS (WT)', 'compact to PS (MT)', 'dispersed to PS (WT)', 'dispersed to PS (MT)'};

for s = 1:4
    nexttile(s);
    hold on

    data = all_datasets{1, s};
    degree_tbl = data{1};
    edge_tbl = data{2};
    G = graph(edge_tbl.EndNodes_1, edge_tbl.EndNodes_2);

    p = plot(G, 'Layout', 'force3', 'Iterations', 300, 'WeightEffect', 'direct', 'UseGravity', false);
    X = p.XData; Y = p.YData; Z = p.ZData; delete(p)

    % Set title above plot
        title(system_labels{s}, ...
              'FontSize', 14, ...
              'Interpreter', 'latex');

    spread_factor = 1.0;
    X = spread_factor * X; Y = spread_factor * Y; Z = spread_factor * Z;

    % Count high-weight edges per node
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

    % Node strengths and coloring
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

    % Draw nodes
    for n = 1:numnodes(G)
        [sx, sy, sz] = sphere(20);
        surf(node_radii(n)*sx + X(n), ...
             node_radii(n)*sy + Y(n), ...
             node_radii(n)*sz + Z(n), ...
             'FaceColor', node_colors(n,:), ...
             'EdgeColor', 'none');
    end

    axis off
    %title(system_labels{s}, 'FontWeight', 'bold', 'Interpreter', 'latex');
end
%
%
cb_ax = axes('Position', [0.89 0.8 0.013 0.15]);  % Dummy axis (still needed for colormap)
colormap(cb_ax, cmap_nodes);
clim(cb_ax, [0 1]);

cb = colorbar(cb_ax, 'Location', 'eastoutside');
cb.Ticks = [0 0.25 0.5 0.75 1];
cb.TickLabels = {'0','0.25','0.5','0.75','1'};
cb.Label.String = 'node strength';
cb.Label.FontSize = 12;
cb.Label.Interpreter = 'latex';

% Turn off dummy axis
axis(cb_ax, 'off');

% Resize colorbar manually by accessing its position
cb_pos = cb.Position;
cb_pos(3) = 0.013;  % Increase width (default is ~0.01)
cb.Position = cb_pos;
cb.Box = 'off';
%%
%
% Save as high-resolution TIFF
save_filename_tif = 'network_graph_300K.tif';
print(gcf, save_filename_tif, '-dtiff', '-r300');  % 600 DPI

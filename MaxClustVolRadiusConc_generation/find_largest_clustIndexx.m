function [longest_path_arr, largest_components] = find_largest_clustIndexx(fileName, threshold_1)
    % Call ChainPairDFF function to get ChainPair and significant node count
    ChainPair = ChainPairDFF(fileName, threshold_1);
    
    % Initialize variables to store adjacency matrices
    Chain_Adjacency_matrices = {};

    % Convert ChainPair to adjacency matrices based on the threshold
    for i = 1:numel(ChainPair)
        logical_indices{i} = ChainPair{i} < threshold_1 & ChainPair{i} ~= 0;
        Chain_Adjacency_matrices{i} = double(logical_indices{i});
    end 

    % Initialize arrays to store results
    longest_path_arr = zeros(numel(Chain_Adjacency_matrices), 2); % To store network number and largest component size
    largest_components = cell(numel(Chain_Adjacency_matrices), 1); % To store node indices of largest components

    % Loop over adjacency matrices to find largest connected component
    for i = 1:numel(Chain_Adjacency_matrices)
        longest_path_arr(i, 1) = i; % Store network number
        largestComponentNodes = find_largest_component(Chain_Adjacency_matrices{i}); % Find largest component nodes
        
        % Store the size of the largest component
        longest_path_arr(i, 2) = numel(largestComponentNodes); % Size of the largest component
        
        % Store the node indices of the largest connected component
        largest_components{i} = largestComponentNodes;
    end
end

% Helper function to find the largest connected component in the adjacency matrix
function largestComponentNodes = find_largest_component(adjacency_matrix)
    % Use the same logic from previous versions to calculate the largest component
    num_nodes = size(adjacency_matrix, 1);
    Links = {};
    
    % Build the connection links for each node
    for k = 1:num_nodes
        nodes = 0;
        for i = 1:num_nodes
            if adjacency_matrix(k, i) == 1
                nodes = nodes + 1;
                Links{k}(1, nodes) = i;
            end
        end
    end
    
    % Initialize connections
    Connections = {};
    for i = 1:numel(Links)
        Connections{i}(1, 1) = i;
        [~, col] = size(Links{i});
        for j = 1:col
            Connections{i}(1, j + 1) = Links{i}(1, j);
        end
    end

    % Merge overlapping connections
    for i = 1:numel(Connections)
        for j = i+1:numel(Connections)
            combined = intersect(Connections{i}, Connections{j});
            if ~isempty(combined)
                Connections{i} = union(Connections{i}, Connections{j});
                Connections{j} = [];
            end
        end
    end

    % Find the largest connected component
    largestComponentNodes = [];
    max_size = 0;
    for i = 1:numel(Connections)
        if ~isempty(Connections{i})
            component_size = numel(Connections{i});
            if component_size > max_size
                largestComponentNodes = Connections{i};
                max_size = component_size;
            end
        end
    end
end

function network_density = calculate_network_density(A)
    % Calculate the number of nodes in the graph
    n = size(A, 1);
    
    % For undirected graph
    % Calculate the total number of possible edges
    total_possible_edges = n * (n - 1) / 2;
    
    % For directed graph
    % total_possible_edges = n * (n - 1);
    
    % Calculate the number of edges present in the graph
    num_edges = nnz(A) / 2; % For undirected graph
    
    % Calculate network density
    network_density = num_edges / total_possible_edges;
end
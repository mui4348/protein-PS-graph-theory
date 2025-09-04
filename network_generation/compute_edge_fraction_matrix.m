function exy = compute_edge_fraction_matrix(A, x, y)
    % Compute the total number of edges in the graph
    m = nnz(A) / 2; % For undirected graphs
    
    % Initialize edge fraction matrix
    exy = zeros(length(x), length(y));
    
    % Iterate over all pairs of degrees in x and y
    for i = 1:length(x)
        for j = 1:length(y)
            % Find nodes with degree x(i) and y(j)
            nodes_x = find(sum(A, 2) == x(i));
            nodes_y = find(sum(A, 2) == y(j));
            % Count the number of edges connecting nodes with degrees in x(i) and y(j)
            exy(i, j) = nnz(A(nodes_x, nodes_y) == 1);
        end
    end
    
    % Normalize the edge fraction matrix
    exy = exy / m;
end
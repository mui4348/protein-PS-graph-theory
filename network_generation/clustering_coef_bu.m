function cc = clustering_coef_bu(A)
    % Number of nodes in the graph
    n = size(A, 1);
    
    % Initialize clustering coefficient vector
    cc = zeros(n, 1);
    
    % Iterate over each node
    for i = 1:n
        % Find neighbors of the current node
        neighbors = find(A(i, :));
        k = numel(neighbors);
        
        % Check if the node has at least 2 neighbors
        if k >= 2
            % Calculate the number of edges between neighbors
            num_edges = sum(A(neighbors, neighbors), 'all') / 2;
            
            % Calculate the maximum possible number of edges between neighbors
            max_edges = k * (k - 1) / 2;
            
            % Calculate clustering coefficient for the current node
            cc(i) = num_edges / max_edges;
        end
    end
    
    % Calculate the average clustering coefficient
    cc = mean(cc(~isnan(cc)));
end
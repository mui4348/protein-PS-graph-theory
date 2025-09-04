function ge = global_efficiency(A)
    % Number of nodes in the graph
    n = size(A, 1);
    
    % Initialize total inverse shortest path length
    total_inv_shortest_path_length = 0;
    
    % Iterate over each pair of nodes
    for i = 1:n
        for j = 1:n
            if i ~= j
                % Compute shortest path length between node i and node j
                shortest_path_length = dijkstra(A, i, j);
                
                % If there is a path between nodes i and j
                if ~isinf(shortest_path_length)
                    % Add the inverse shortest path length to the total
                    total_inv_shortest_path_length = total_inv_shortest_path_length + 1 / shortest_path_length;
                end
            end
        end
    end
    
    % Calculate global efficiency
    ge = total_inv_shortest_path_length / (n*(n-1));
end
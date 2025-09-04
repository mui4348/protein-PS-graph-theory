function length = dijkstra(A, start, dest)
    % Initialize variables
    n = size(A, 1);
    visited = false(1, n);
    distance = Inf(1, n);
    distance(start) = 0;
    
    % Dijkstra's algorithm
    for i = 1:n
        % Find the vertex with the minimum distance
        [~, u] = min(distance .* ~visited);
        visited(u) = true;
        
        % Update distances to neighbors of u
        for v = 1:n
            if A(u, v) ~= 0 && ~visited(v)
                distance(v) = min(distance(v), distance(u) + A(u, v));
            end
        end
    end
    
    % Return the shortest path length from start to dest
    length = distance(dest);
end
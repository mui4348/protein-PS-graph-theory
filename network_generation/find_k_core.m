function ccc = find_k_core(A)
    % Initialize k-core with the original graph
   k_core = A;
    k = 0;
    % Compute the degree of each node
    degrees = sum(A);
    
    % Keep track of nodes with degree less than k
    low_degree_nodes = find(degrees < k);
    [r,c] = size(A);
    % Remove nodes with degree less than k and update degrees
    
    
    while size(low_degree_nodes) ~= r
        % Remove low-degree nodes
        k_core(low_degree_nodes, :) = 0;
        k_core(:, low_degree_nodes) = 0;
        
        % Recompute degrees after removal
        degrees = sum(k_core);
        
        % Update low-degree nodes
        low_degree_nodes = find(degrees < k);
k = k+1;
    end
ccc = k;
end
% will return 4 things. 
% 1: no of connections, 
% 2: list of 1, 
% 3: number of components including isolated nodes and connections, 
% 4: List of 3.
function [Links] = count_links(adjacency_matrix)
    num_nodes = size(adjacency_matrix, 1);
    node_list = 1:num_nodes;
    Links = {};
    % Compute unique paths using dynamic programming
    paths = 0;
    k = 1;
    paths_local = zeros(1,num_nodes);
    cancel = zeros(1,num_nodes);
      for k = 1:num_nodes
        nodes = 0;
        for i = 1:num_nodes
            if adjacency_matrix(k,i) == 1
                    nodes = nodes + 1;
                    Links{k}(1,nodes) = i;
            end
        end
      end
end

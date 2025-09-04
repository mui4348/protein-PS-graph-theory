function [indices_row,indices_col] = calculatePairDataidices(coords,threshold,ch1,ch2)
chain_size = 33;
% ch1 = k;
% ch2 = l;
no_of_chains = size(coords, 1)/33;
% chain_distance = zeros(no_of_chains,no_of_chains);
node_dist = zeros(chain_size,chain_size);
ch1_address = chain_size*(ch1-1);
ch2_address = chain_size*(ch2-1);
        for ch1_nod = 1:chain_size
            for ch2_nod = 1:chain_size
                ch1_node = ch1_address + ch1_nod;
                ch2_node = ch2_address + ch2_nod;
                node_dist(ch1_nod,ch2_nod) = norm(coords(ch1_node, :) - coords(ch2_node, :));
            end
        end
        
        [indices_row,indices_col] = find(node_dist <= threshold);
end



function [chain_distance,count] = calculatePairData(coords,threshold)
    % % Calculate pair distances and assign corresponding residues and indices to each pair
    % 
    % numCoords = size(coords, 1);
    % numPairs = sum(1:numCoords-3); % Calculate the total number of pairs
    % 
    % pairData = cell(numPairs, 3); % Initialize pair data cell array
    % for i = 1:numCoords
    %     for j = i:numCoords
    %         pairData{k, 1} = i;
    %         % pairData{k, 2} = residues{i};
    %         pairData{k, 2} = j;
    %         % pairData{k, 4} = residues{j};
    %         pairData{k, 3} = norm(coords(i, :) - coords(j, :));
    %         k = k + 1;
    %     end
    % end
chain_size = 33;
no_of_chains = size(coords, 1)/33;
chain_distance = zeros(no_of_chains,no_of_chains);
node_dist = zeros(chain_size,chain_size);
for ch1 = 0:no_of_chains-1
    for ch2 = ch1+1:no_of_chains-1
        for i = 1:chain_size
            for j = 1:chain_size
                ch1_node = i+ch1*chain_size;
                ch2_node = j+ch2*chain_size;
                node_dist(i,j) = norm(coords(ch1_node, :) - coords(ch2_node, :));
            end
        end
        minimum_row = min(node_dist);
        chain_distance(ch1+1,ch2+1) = min(minimum_row);
        count(ch1+1,ch2+1) = sum(node_dist(:) < threshold);
    end
end
end



function[Linkss] = chainLinkss(fileName,threshold_1)
ChainPair = ChainPairDFF(fileName,threshold_1);
Chain_Adjacency_matrices = {};

for i=1:numel(ChainPair)
    logical_indices{i} = ChainPair{i} < threshold_1 & ChainPair{i} ~= 0;
    Chain_Adjacency_matrices{i} = double(logical_indices{i});
    [Links{i}] = count_links(Chain_Adjacency_matrices{i});
    nonempty_count = 1;
for j = 1:numel(Links{i})
    if ~isempty(Links{i}{j})
        Linkss{i}{nonempty_count} = [j,Links{i}{j}];
        nonempty_count = nonempty_count+1;
    end
end
end
end
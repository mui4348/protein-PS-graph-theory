function model = extractModel(lines, startIndices, endIndices,initial,final)
    % Extract atom lines for each model
    
    model = cell(numel(startIndices), 1);
    for i = initial:final
        model{i} = lines(startIndices(i)+1:endIndices(i)-1);
    end
end
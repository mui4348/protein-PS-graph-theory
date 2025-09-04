function [pairData,Sig_node_count] = ChainPairDFF(fileName,threshold)
% fileName='MaSpI1_1C_viz_prot_nopbc_101_frames.pdb';

segment=1;
k=1;
frames = 30;
initial = frames*segment - frames + 1;
% fileName='test_mmm.pdb';
% k=1;
       % Calculate pair distances from a PDB file
    file = fopen(fileName);
    lines = textscan(file, '%s', 'Delimiter', '\n', 'CommentStyle', {'REMARK'});
    lines = lines{1}; % Extract lines from cell array
    fclose(file);

    % Remove records
    lines = removeRecords(lines, 'TER');
    lines = removeRecords(lines, 'ANISOU');
    lines = removeRecords(lines, 'HETATM'); % Remove heteroatoms

    % Find indices
    [modelStart, modelEnd] = findModelIndices(lines);
    frames=numel(modelEnd);
mod = numel(modelEnd);
    final = mod;
    % frames = min(frames,numel(modelEnd)-initial+1);
    % Extract models
    if isempty(modelStart)
        [modelStart, modelEnd] = findAtomIndices(lines);
        model = extractModel(lines, modelStart, modelEnd,initial,final);
    elseif numel(modelStart) == numel(modelEnd)
        model = cell(frames, 1);
        for i = initial:final
            model{i-initial+1} = lines(modelStart(i)+1:modelEnd(i)-1);
        end
    else
        error('Failed to determine model numbers');
    end
    if frames<=10
        cellcount=1;
    else
        cellcount = 10;
    end
pairData = cell(1, numel(model));
Sig_node_count = cell(1, numel(model));
    % Initialize pair data cell array
    % pairData = {};

    % Loop through models
    for i = 1:numel(model)
        disp(lines(modelStart(i)));
        disp(i/numel(model)*100)
        atom = cell(1, numel(model{i}));
        residues = cell(1, numel(model{i}));
        coords = zeros(numel(model{i}), 3); % Initialize coordinates matrix
        for j = 1:numel(model{i}) % Loop through lines of each model
            line = model{i}{j};
            atom{j} = upper(strtrim(line(13:15))); % Extract atom identifier, remove spaces
            if strcmp(atom{j}, 'BB')
                residues{j} = upper(strtrim(line(17:20))); % Extract residue name, remove spaces
                coords(j, :) = sscanf(line(31:54), '%f', [1, 3]); % Read coordinates
            end
        end
        residues = transpose(residues(~cellfun('isempty', residues)));

        % Remove rows with all zero coordinates
        coords = coords(any(coords, 2), :);

        % Ensure residues and coords have compatible dimensions
        if size(residues, 1) ~= size(coords, 1)
            error('Residues and coordinates dimensions are not compatible.');
        end

        % Calculate pair distances and assign corresponding residues and indices to each pair
        [pairDataModel,Sig_node_count_model] = calculatePairData(coords,threshold);
        pairData{i} = pairDataModel'+pairDataModel;
        Sig_node_count{i} = Sig_node_count_model;
        % Combine data from all models
        % pairData{jjj} = [pairData{jjj}; pairDataModel];
    end
end
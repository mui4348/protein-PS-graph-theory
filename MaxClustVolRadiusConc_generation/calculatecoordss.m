%% This function also extract chainIDs along with coords and residues for each pdb file
function [coordss, residuess, chainIDs] = calculatecoordss(fileName) 
% fileName='MaSpI1_1C_viz_prot_nopbc_101_frames.pdb';

initial = 1;
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
frames = numel(modelEnd);
mod = numel(modelEnd);
final = mod;

% Extract models
if isempty(modelStart)
    [modelStart, modelEnd] = findAtomIndices(lines);
    model = extractModel(lines, modelStart, modelEnd, initial, final);
elseif numel(modelStart) == numel(modelEnd)
    model = cell(frames, 1);
    for i = initial:final
        model{i - initial + 1} = lines(modelStart(i) + 1 : modelEnd(i) - 1);
    end
else
    error('Failed to determine model numbers');
end

pairData = cell(1, numel(model));
Sig_node_count = cell(1, numel(model));

for i = 1:numel(model)
    disp(lines(modelStart(i)));
    disp(i / numel(model) * 100)
    atom = cell(1, numel(model{i}));
    residues = cell(1, numel(model{i}));
    coords = zeros(numel(model{i}), 3); % Initialize coordinates matrix
    chainIDs = cell(1, numel(model{i})); % Initialize chain ID array
    
    for j = 1:numel(model{i}) % Loop through lines of each model
        line = model{i}{j};
        atom{j} = upper(strtrim(line(13:15))); % Extract atom identifier, remove spaces
        if strcmp(atom{j}, 'BB') % Only consider backbone atoms
            residues{j} = upper(strtrim(line(17:20))); % Extract residue name, remove spaces
            coords(j, :) = sscanf(line(31:54), '%f', [1, 3]); % Read coordinates
            chainIDs{j} = strtrim(line(22)); % Extract chain ID from column 5 (position 22)
        end
    end
    
    % Clean up empty entries
    residues = transpose(residues(~cellfun('isempty', residues)));
    chainIDs = transpose(chainIDs(~cellfun('isempty', chainIDs)));
    
    % Remove rows with all zero coordinates
    coords = coords(any(coords, 2), :);

    % Ensure residues and coords have compatible dimensions
    if size(residues, 1) ~= size(coords, 1)
        error('Residues and coordinates dimensions are not compatible.');
    end
    
    coordss{i} = coords;
    residuess{i} = residues(1:33); % Extract only the first 33 residues
    chainIDs{i} = chainIDs(1:33);  % Extract the corresponding chain IDs
end
end

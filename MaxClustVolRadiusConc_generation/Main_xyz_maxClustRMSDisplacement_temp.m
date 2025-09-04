clear;

% Parameters
time_step = 10;        % Timestep between frames in ns
num_frames = 1001;     % Maximum number of frames (corresponds to 10000 ns)
num_trajectories = 3;  % Number of trajectories
lag_time = 100;        % Fixed lag time (in ns)
num_chains = 60;       % Number of chains in each trajectory
residues_per_chain = 33; % Number of residues per chain

% Define file names for the cases
Names = { 'wt', 'wt_unpack', 'mt', 'mt_unpack' };
for i = 270:10:340 % Temperature (K)
% i = 330;  % Temperature (K)
% Loop over each file
for k = 1:numel(Names)
    Name = Names{k};

    % Initialize storage for concatenated centroids across trajectories
    concatenated_centroids = struct('x', [], 'y', [], 'z', []);
    % Preallocate centroids for all chains and trajectories
total_frames = num_trajectories * num_frames;
concatenated_centroids.x = nan(num_chains, total_frames);
concatenated_centroids.y = nan(num_chains, total_frames);
concatenated_centroids.z = nan(num_chains, total_frames);

    % Loop through the trajectories and concatenate centroids
    for j = 1:num_trajectories
    % Set the file name based on temperature and trajectory
    temp = int2str(i);
    traj = int2str(j);
    fileName = strcat(Name, '_', traj, 'C_MaxClust_', temp, 'K.pdb');

    % Extract coordinates from the PDB file
    [coords, ~, ~] = calculatecoordsss(fileName); % coords is 1001x1 cell array

    % Calculate frame offset for the current trajectory
    frame_offset = (j - 1) * num_frames;

    % Loop over frames
    for frame = 1:numel(coords)
        frame_coords = coords{frame}; % Extract frame's coordinates (Nx3 matrix)

        % Dynamically compute number of chains in current frame
        num_residues = size(frame_coords, 1);
        actual_num_chains = floor(num_residues / residues_per_chain);

        for chain_idx = 1:actual_num_chains
            start_idx = (chain_idx - 1) * residues_per_chain + 1;
            end_idx = chain_idx * residues_per_chain;
            % Only compute if the indices are within bounds
            if end_idx <= size(frame_coords, 1)
               chain_coords = frame_coords(start_idx:end_idx, :);

            % Store only if within expected max chains and frames
            %if chain_idx <= num_chains && (frame_offset + frame) <= total_frames
                concatenated_centroids.x(chain_idx, frame_offset + frame) = mean(chain_coords(:, 1));
                concatenated_centroids.y(chain_idx, frame_offset + frame) = mean(chain_coords(:, 2));
                concatenated_centroids.z(chain_idx, frame_offset + frame) = mean(chain_coords(:, 3));
            %end
            end
         end
    end
    end

    % Initialize structure to store MSD
    aggregated_msd = struct('x', [], 'y', [], 'z', []);

    % Calculate MSD for concatenated data using sliding window
   frame_step = round(lag_time / time_step);
   max_shift = size(concatenated_centroids.x, 2) - frame_step;

for start_frame = 1:max_shift
    target_frame = start_frame + frame_step;

    % if target_frame > size(concatenated_centroids.x, 2)
    %     break; % Safeguard against out-of-bounds indexing
    % end

    % Compute MSD for each chain
    for chain_idx = 1:num_chains
        x1 = concatenated_centroids.x(chain_idx, start_frame);
        x2 = concatenated_centroids.x(chain_idx, target_frame);
        y1 = concatenated_centroids.y(chain_idx, start_frame);
        y2 = concatenated_centroids.y(chain_idx, target_frame);
        z1 = concatenated_centroids.z(chain_idx, start_frame);
        z2 = concatenated_centroids.z(chain_idx, target_frame);

        % Skip if any coordinate is NaN (i.e., missing chain data)
        if any(isnan([x1, x2, y1, y2, z1, z2]))
            continue;
        end

        dx = x2 - x1;
        dy = y2 - y1;
        dz = z2 - z1;

        % Square and convert from Å² to nm²
        aggregated_msd.x(chain_idx, start_frame) = (dx^2) * 0.01;
        aggregated_msd.y(chain_idx, start_frame) = (dy^2) * 0.01;
        aggregated_msd.z(chain_idx, start_frame) = (dz^2) * 0.01;
    end
end

    % Average MSD over all chains
    avg_msd_x = mean(aggregated_msd.x, 1, 'omitnan');
    avg_msd_y = mean(aggregated_msd.y, 1, 'omitnan');
    avg_msd_z = mean(aggregated_msd.z, 1, 'omitnan');

    % Populate time values
    %time_values = (1:numel(avg_msd_x)) * time_step;
    time_values = 1:numel(avg_msd_x);

    % Save aggregated MSD data to CSV file
    aggregated_table = table(time_values', avg_msd_x', avg_msd_y', avg_msd_z', ...
        'VariableNames', {'Time_ns', 'Avg_MSD_X_nm2', 'Avg_MSD_Y_nm2', 'Avg_MSD_Z_nm2'});
     writetable(aggregated_table, strcat('Aggregated_MSD_', Name, '_maxClust_Fixlag_', temp, 'K.csv'));
end
 end

%% This function also extract chainIDs along with coords and residues for each pdb file
function [coordss, residuess, chainIDs] = calculatecoordsss(fileName) 
% fileName='mute_1C_viz_prot_nopbc_101_frames.pdb';

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
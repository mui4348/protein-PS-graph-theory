clear; 

% Define file names for the two cases
Names = { 'wt', 'wt_unpack', 'mt', 'mt_unpack' };
num_trajectories = 3;
for ii =  270:10:340
for k = 1:numel(Names)
    Name = Names{k};
% Loop through the 3 trajectories
for j = 1:num_trajectories
    % Define file parameters
    temp = int2str(ii);
    traj = int2str(j);
    threshold = 8;
    fileName = strcat(Name, '_', traj, 'C_protNoPBC_1k_', temp, 'K.pdb');
	display(['Processing file: ', fileName]);

    % Call the find_largest_cluster function to get largest connected component node indices
    [~, Nodes_with_linkss] = find_largest_clustIndexx(fileName, threshold);

    % Calculate coordinates, residues, and chain IDs from the PDB file
    [coords, residues, chainIDs] = calculatecoordss(fileName);

    % Ensure all chainIDs are character vectors
    for n = 1:numel(chainIDs)
        if iscell(chainIDs{n})
            chainIDs{n} = char(chainIDs{n});
        end
    end

    % Extract unique chain IDs with case sensitivity
    chainIDschar = char(chainIDs);
    uniqueChainIDs = unique(chainIDschar, 'stable');

    % Save the unique chain IDs in a list
    uniqueIDList = cellstr(uniqueChainIDs);

    % Create the output file name for the maximum cluster
    output = strcat(Name, '_', traj, 'C_MaxClust_', temp, 'K.pdb');
    outputFile = fopen(output, 'w');
    if outputFile == -1
        error('Could not open output file for writing.');
    end

    % Initialize cell array to store the sub-cell with the maximum number of elements
    max_subCell = cell(1, numel(Nodes_with_linkss));
   all_selectedChainIDs = cell(1, numel(Nodes_with_linkss)); % To store all chainIDs for selected clusters
   num_of_maxClust = zeros(1, numel(Nodes_with_linkss)); % Initialize this array with zeros for each frame
    NumMaxClust = zeros(1, numel(coords));

    % Loop through each frame to determine the sub-cell with the maximum number of elements
    for frameIdx = 1:numel(Nodes_with_linkss)
        currentFrame = Nodes_with_linkss{frameIdx};

        % Check if currentFrame is a cell array or not
        if iscell(currentFrame)
            max_val = -Inf; % Initialize maximum number of elements to negative infinity
            index_max = 1;   % Initialize the index of the max sub-cell

            % Loop through each sub-cell in the current frame (if currentFrame is a cell array)
            for subCellIdx = 1:numel(currentFrame)
                currentSubCell = currentFrame{subCellIdx};
                if ~isempty(currentSubCell)
                    current_max = numel(currentSubCell);
                    if current_max > max_val
                        max_val = current_max;
                        index_max = subCellIdx;
                    end
                end
            end

            % Store the sub-cell with the maximum number of elements for this frame
            max_subCell{frameIdx} = currentFrame{index_max};

        elseif isnumeric(currentFrame)  % If currentFrame is numeric or another data type (e.g., a matrix)
            % Simply find the largest row/column in terms of non-zero elements
            [max_val, index_max] = max(sum(currentFrame, 2)); % Sum over rows (for clusters)
            max_subCell{frameIdx} = currentFrame(index_max, :); % Store the largest cluster row
        end

        % Store number of elements in the largest cluster for this frame
        num_of_maxClust = numel(max_subCell{frameIdx});
        NumMaxClust(frameIdx,1) = num_of_maxClust;
		
	% Map indices to chain IDs
        selectedChainIDs = uniqueIDList(max_subCell{frameIdx});
        % Initialize a 60x1 array with NaN
        frameChainIDs = strings(60, 1);
        frameChainIDs(:) = "NaN";

        % Map selectedChainIDs to the corresponding indices in the 60x1 array
        for idx = 1:numel(uniqueIDList)
            if ismember(uniqueIDList{idx}, selectedChainIDs)
                frameChainIDs(idx) = uniqueIDList{idx};
            end
        end

        % Store the frameChainIDs in all_selectedChainIDs
        all_selectedChainIDs{frameIdx} = frameChainIDs;
		
    end
	
	% Save all_selectedChainIDs to a file
    outputCSVFile = strcat(Name, '_', traj, 'C_maxclustSeleChainIDs_1k_', temp, 'K.csv');
    % Initialize a matrix to store logical values
    logicalMatrix = zeros(numel(all_selectedChainIDs), numel(uniqueIDList));

     for frameIdx = 1:numel(all_selectedChainIDs)
    for idx = 1:numel(uniqueIDList)
        % Set 1 if the ChainID is present in the current frame, otherwise 0
        logicalMatrix(frameIdx, idx) = ismember(uniqueIDList{idx}, all_selectedChainIDs{frameIdx});
    end
    end

    % Create a cell array for the header (first row)
    csvHeader = ['Frame', uniqueIDList'];

    % Write the data to a CSV file
    csvData = [num2cell((1:numel(all_selectedChainIDs))'), num2cell(logicalMatrix)];
    writecell([csvHeader; csvData], outputCSVFile);

    % Open the original PDB file for reading
    fileID = fopen(fileName, 'r');
    if fileID == -1
        error('Could not open file: %s', fileName);
    end

    % Loop through the PDB file line by line
    while ~feof(fileID)
        line = fgetl(fileID);

        % Handle REMARK, TITLE, and CRYST1 lines (write them to the output file)
        if startsWith(line, 'REMARK') || startsWith(line, 'TITLE') || startsWith(line, 'CRYST1')
            fprintf(outputFile, '%s\n', line);
        end

        % Handle MODEL lines and extract model number
        if startsWith(line, 'MODEL')
            fprintf(outputFile, '%s\n', line);
            modelNumberStr = strtrim(line(7:end)); % Extract the part after 'MODEL'
            modelNumber = str2double(modelNumberStr); % Convert string to number

            if isnan(modelNumber)
                error('Invalid MODEL number in line: %s', line);
            end

            frameIdx = modelNumber;
            if frameIdx > numel(max_subCell)
                warning('Frame index exceeds available data: %d', frameIdx);
                continue;
            end

            % Print the max clusters for the current frame
            if ~isempty(max_subCell{frameIdx})
                fprintf('Frame %d: max clusters indices = %s, Number of max clusters = %d\n', ...
                        frameIdx, mat2str(max_subCell{frameIdx}), numel(max_subCell{frameIdx}));

                % Map max_subCell indices to chain IDs
                selectedChainIDs = uniqueIDList(max_subCell{frameIdx});
            else
                selectedChainIDs = []; % Empty in case no cluster exists
            end
        end

        % Handle ATOM lines (filter based on chain ID)
        if startsWith(line, 'ATOM') && length(line) >= 22 && contains(line, 'BB'))

            chainID = strtrim(line(22)); % Extract chain ID from 22nd column

            % Check if the chain ID is in the selectedChainIDs for the current frame
            if ismember(chainID, selectedChainIDs)
                fprintf(outputFile, '%s\n', line);
            end
        end

        % Handle TER and ENDMDL lines (write them to the output file)
        if startsWith(line, 'TER') || startsWith(line, 'ENDMDL')
            fprintf(outputFile, '%s\n', line);
        end
    end

    fclose(fileID); % Close the original PDB file
    fclose(outputFile); % Close the output PDB file

    [ClustCentroids, clustRadii] = cal_Clustcentroid_radii(output);

    % Store centroids and radii for each frame
    centroidsForTraj(1:frameIdx, :) = ClustCentroids;  % Store the centroid for this frame
    radiiForTraj(1:frameIdx) = clustRadii;            % Store the radius for this frame
	original_clust_Radius = (radiiForTraj(1:frameIdx)); % Original radius from the file
    Shell_radii = (radiiForTraj(1:frameIdx));  % Original radius for shells
    Shell_radii = Shell_radii - Shell_radii + 10;
    
    %%
    no_of_shells = 40; % Defines four concentric spherical shells around the cluster.
    box_size = 200;
    radiiArray11 = zeros(frameIdx,no_of_shells);
    max_poss_dist_frm_centroid = zeros(1,frameIdx);
    V_inner = zeros(frameIdx,1);
    V_outer_of_largest_cluster_outer = zeros(frameIdx,no_of_shells);
    for frm = 1:frameIdx
    %radiiArray11(frm,:) = originalRadius11(1,frm):25:200;% Create an array from original radius up to 200 angstrom
    max_poss_dist_frm_centroid(1,frm) = farthestCornerDistance(centroidsForTraj(frm,1),centroidsForTraj(frm,2),centroidsForTraj(frm,3),box_size);
    radiiArray11(frm,:) = linspace(Shell_radii(1,frm), max_poss_dist_frm_centroid(frm), no_of_shells);
    %V_inner = (4 / 3) * pi .* (radiiArray11.^3) / 10^24; % Volume of the sphere, convert to liters, 1 liter = 10^24 cubic nm;
    V_inner (frm,1) = sphereBoxIntersection_M_C(centroidsForTraj(frm,1),centroidsForTraj(frm,2),centroidsForTraj(frm,3),Shell_radii(1,frm),box_size)/10^24;
    volumes_in = V_inner * 1000; % convert L to ml 
    for i = 1:no_of_shells
    %V_outer_of_largest_cluster(:,i) = (4 / 3) * pi .* (radiiArray11(:,i).^3 - transpose(originalRadius11).^3) / 10^24; % Volume of Outer Spherical Shell (Difference of Two Spheres)
    V_outer_of_largest_cluster_outer(frm,i) = sphereBoxIntersection_M_C(centroidsForTraj(frm,1),centroidsForTraj(frm,2),centroidsForTraj(frm,3),radiiArray11(frm,i),box_size)/10^24;
    end
    disp(['processed Frame: ', num2str(frm)]);
    end
    V_outer11 = zeros(frameIdx,no_of_shells-1);
    for i = 1:no_of_shells-1
        ccc = no_of_shells - i + 1;
        V_outer11(:,ccc-1) = V_outer_of_largest_cluster_outer(:,ccc)-V_outer_of_largest_cluster_outer(:,ccc-1);
    end
    V_outer_region = V_outer11 * 1000;

    V_clust_inner = zeros(frameIdx, 1); % Initialize inner volume storage
    
    for frm = 1:frameIdx
        % Get the radius for this frame (single value instead of shells)
        cluster_radius = original_clust_Radius(1, frm);
    
        % Compute the inner volume, considering box boundaries
        V_clust_inner(frm, 1) = sphereBoxIntersection_M_C(...
            centroidsForTraj(frm, 1), centroidsForTraj(frm, 2), centroidsForTraj(frm, 3), ...
            cluster_radius, box_size) / 10^24; % Convert Å³ to liters
    end

    % Convert to mL
    volumes_clust_in = V_clust_inner * 1000; % Convert L to mL

    %%

    % Reopen the original PDB file for reading
    fileID = fopen(fileName, 'r');
    if fileID == -1
        error('Could not open file: %s', fileName);
    end
    distance = zeros(1,numel(chainIDs)*frameIdx);
            frameiddx = 1;    
            stmcnt=0;

            %The following while loop defines the distance array, which
            %consists of the corresponding distances of all atoms from the
            %centroid of the cluster

            while ~feof(fileID) && stmcnt < numel(chainIDs)*frameIdx
                frameiddx = floor(stmcnt/numel(chainIDs))+1;
                centroid = centroidsForTraj(frameiddx, :);
                line = fgetl(fileID);
                if startsWith(line, 'ATOM') && length(line) >= 22 && contains(line, 'BB')
                %if startsWith(line, 'ATOM') && length(line) >= 22 && (contains(line, 'BB') || contains(line, 'SC'))
                    % Get the atom's coordinates (columns 31-38 for x, 39-46 for y, 47-54 for z)
                    atomX = str2double(line(31:38));
                    atomY = str2double(line(39:46));
                    atomZ = str2double(line(47:54));

                    % Calculate distance from the centroid
                    stmcnt = stmcnt + 1;
                    dist = sqrt((atomX - centroid(1))^2 + (atomY - centroid(2))^2 + (atomZ - centroid(3))^2);
                    distance(stmcnt) = dist;
                    
                end
            end



    clust_atomCount = zeros(frameIdx, 1); % Only one count per frame
    
    for frameIdx = 1:numel(Nodes_with_linkss)
        for atmcnt = (frameIdx-1)*numel(chainIDs)+1:(frameIdx)*numel(chainIDs)
            % Check if atom is within the single radius
            if distance(1, atmcnt) <= original_clust_Radius(frameIdx)
                clust_atomCount(frameIdx) = clust_atomCount(frameIdx) + 1;
            end
        end
    end



    atomCount = zeros(frameIdx,no_of_shells);
    %The following loop calculates the number of atoms in each shell
    for frameIdx = 1:numel(Nodes_with_linkss)
            for atmcnt = (frameIdx-1)*numel(chainIDs)+1:(frameIdx)*numel(chainIDs)
                shell = 1;
                while shell <= no_of_shells
                    
                    % If the distance is within the radius, include the atom
                    if distance(1,atmcnt) <= radiiArray11(frameIdx,shell)
                        % fprintf(outputFilee, '%s\n', line);
                        atomCount(frameIdx,shell) = atomCount(frameIdx,shell) + 1;  % Increment atom count
                        break
                    else 
                        shell = shell+1;
                    end
                end
                        % Check if the atom is inside the original radius
            end
   end
 N_atoms_per_shell = atomCount *2;
 % V_shell (:,2:4) = V_outer_region;
 V_shell = zeros(size(V_outer_region,1), size(V_outer_region,2)+1); % Adjusting size
 V_shell(:, 2:end) = V_outer_region; % Assigning values
 
 V_shell (:,1) = volumes_in(:,1);
 concentration_shell = N_atoms_per_shell ./ V_shell;
 concentration_mol_perML_inner = concentration_shell / (6.022e23); % Divides by Avogadro’s number (6.022e23 molecules/mol) to convert atoms to moles.
 inner_shell_conc_Mol_perML = concentration_mol_perML_inner * 1000; % Multiplies by 1000 to get mol/mL/.
 
 N_atoms_per_clust = clust_atomCount *2;
 clust_concentration = N_atoms_per_clust ./ volumes_clust_in;
 clust_concentration_Mol_perML = clust_concentration / (6.022e23); % Divides by Avogadrols number (6.022e23 molecules/mol) to convert atoms to moles.
 clust_conc_Mol_perML = clust_concentration_Mol_perML * 1000; % Multiplies by 1000 to get mol/mL/.

 str = strcat('Trajectory  ', traj);
 writematrix(inner_shell_conc_Mol_perML, strcat( Name, '_concentrations_', temp,'K.xlsx'),'Sheet',str);
 writematrix(N_atoms_per_shell,strcat( Name, '_atomCount_', temp,'K.xlsx'),'Sheet',str);
 writematrix( V_shell,strcat( Name, '_volumes_', temp,'K.xlsx'),'Sheet',str);
 writematrix( radiiArray11,strcat( Name, '_radii_', temp,'K.xlsx'),'Sheet',str);

 fclose(fileID); % Close the original PDB file
%end

% Define the number of frames dynamically
numFrames = numel(clustRadii);

% Define the Excel file name for the current temperature
outputExcelFile = strcat(Name, '_maxClust_SizeRadiusVolConcAtoms_', temp, 'K.xlsx');

    clusterData = cell(numFrames, 6); % 1 index + 5 data columns
    
    % Define the header for the current trajectory
    header = {'Frame_Index', ...
        'Largest_Cluster_Size', 'Radius_nm', 'Volume_nm_3', ...
        'Number_of_Atoms', 'Max_Clust_Cs_mM'};

    % Fill in the data for the current trajectory
    for frameIdx = 1:numFrames
        clusterData{frameIdx, 1} = frameIdx; % Frame Index
        clusterData{frameIdx, 2} = NumMaxClust(frameIdx);
        clusterData{frameIdx, 3} = clustRadii(frameIdx);
        clusterData{frameIdx, 4} = volumes_clust_in(frameIdx);
        clusterData{frameIdx, 5} = N_atoms_per_clust(frameIdx);
        clusterData{frameIdx, 6} = clust_conc_Mol_perML(frameIdx);
    end

    % Define the sheet name for this trajectory
    sheetName = strcat('Trajectory ', num2str(j));
    
    % Write header and data to the respective sheet in the Excel file
    writecell(header, outputExcelFile, 'Sheet', sheetName, 'Range', 'A1');
    writecell(clusterData, outputExcelFile, 'Sheet', sheetName, 'Range', 'A2');

end
end
 end

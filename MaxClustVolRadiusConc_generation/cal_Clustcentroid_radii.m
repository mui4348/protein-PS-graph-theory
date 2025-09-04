function [centroids, radii] = cal_Clustcentroid_radii(pdbFile)
    % Analyze clusters in a PDB file and return the centroid and radius for the largest cluster in each frame
    %
    % Inputs:
    %   pdbFile - Path to the PDB file
    % Outputs:
    %   centroids - Centroids of the largest cluster in each frame (Nx3)
    %   radii - Radii of the largest cluster in each frame (Nx1)

    % Open the PDB file for reading
    file = fopen(pdbFile, 'r');
    if file == -1
        error('Could not open the file.');
    end
    lines = textscan(file, '%s', 'Delimiter', '\n', 'CommentStyle', {'REMARK'});
    fclose(file);
    lines = lines{1};

    % Remove unwanted records like TER, ANISOU, and HETATM
    lines = removeRecords(lines, 'TER');
    lines = removeRecords(lines, 'ANISOU');
    lines = removeRecords(lines, 'HETATM');

    % Find the indices for models in the PDB file
    [modelStart, modelEnd] = findModelIndices(lines);
    frames = numel(modelEnd);

    if isempty(modelStart)
        [modelStart, modelEnd] = findAtomIndices(lines);
        model = extractModel(lines, modelStart, modelEnd, 1, frames);
    elseif numel(modelStart) == numel(modelEnd)
        model = cell(frames, 1);
        for k = 1:frames
            model{k} = lines(modelStart(k)+1:modelEnd(k)-1);
        end
    else
        error('Failed to determine model numbers');
    end
    
    % Prepare to store the results
    centroids = NaN(frames, 3);  % Ensure centroids is a matrix (frames x 3)
    radii = NaN(frames, 1);      % Ensure radii is a column vector (frames x 1)

    % Loop through each frame to calculate centroid and radius
    for k = 1:frames
        clustcoords = []; % Coordinates for the largest cluster in the current frame
        atomCount = 0;

        % Loop over the atoms in the current frame (model{k})
        for m = 1:numel(model{k})
            line = model{k}{m};
            atom = strtrim(line(13:15));  % Check atom type (e.g., BB)

            % For now, you can add more conditions for selecting atoms, 
            % e.g., if you want to select specific residues or atoms
            if strcmpi(atom, 'BB') % Backbone atoms (or change as needed)
            %if strcmpi(atom, 'BB') || contains(atom, 'SC') 
                atomCount = atomCount + 1;
                clustcoords(atomCount, :) = sscanf(line(31:54), '%f', [1, 3]);
            end
        end
        
        % Ensure clustcoords only contains valid data
        clustcoords = clustcoords(1:atomCount, :);
        
        % Check if any atoms were found for this frame
        if atomCount > 0
            % Calculate the centroid as the mean of the coordinates of atoms in the largest cluster
            centroid = mean(clustcoords, 1);
            centroids(k, :) = centroid;  % Store the centroid for this frame

            % Calculate the distances from the centroid and find the mean radius
            distances = sqrt(sum((clustcoords - centroid).^2, 2));
            % radii(k) = (mean(distances))* 0.1; % Convert radius from Å to nm;
            radii(k) = mean(distances); % Convert radius from Å to nm;
        else
            % If no atoms were found, store NaN (or handle as needed)
            centroids(k, :) = NaN;
            radii(k) = NaN;
        end
    end
end



function [startIndices, endIndices] = findModelIndices(lines)
    % Find start and end indices of each model entry
    
    startIndices = find(strncmp(lines, 'MODEL ', 6)); % Start of each model entry
    endIndices = find(strncmp(lines, 'ENDMDL', 6)); % End of each model entry
end
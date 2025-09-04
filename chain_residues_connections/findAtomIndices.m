function [startIndices, endIndices] = findAtomIndices(lines)
    % Find start and end indices of atom lines
    
    startIndices = find(strncmp(lines(1:round(end/2)), 'ATOM  ', 6), 1); % Start of atom lines
    endIndices = find(strncmp(lines(startIndices:end), 'ATOM  ', 6), 1) - 1; % End of atom lines
end
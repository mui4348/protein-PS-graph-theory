function cleanedLines = removeRecords(lines, recordType)
    % Remove records of specified type from the lines
    
    ix = strncmp(lines, recordType, numel(recordType));
    cleanedLines = lines(~ix);
end
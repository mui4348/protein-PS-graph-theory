function [x, y] = compute_degrees(A)
    % Calculate the degree of each node
    degrees = sum(A);
    
    % Construct set x (degrees of nodes)
    x = unique(degrees);
    
    % Construct set y (same as set x)
    y = x;
end
function r = degree_assortativity_coefficient(exy, ax, by, sigma_a, sigma_b)
    % Compute the sum over xy pairs
    sum_xy = sum(sum(exy));
    
    % Compute the sum of products of xy, exy, ax, and by
    sum_product = sum(sum((exy - ax' * by) .* (1:length(ax))' .* (1:length(by))));
    
    % Compute the degree assortativity coefficient
    r = sum_product / (sigma_a * sigma_b);
end
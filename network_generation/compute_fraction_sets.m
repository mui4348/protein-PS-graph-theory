function [ax, by] = compute_fraction_sets(exy)
    % Sum along rows and columns of the edge fraction matrix
    ax = sum(exy, 2);
    by = sum(exy, 1)';
end
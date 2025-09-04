function [sigma_a, sigma_b] = compute_strand_deviations(ax, by)
    sigma_a = std(ax);
    sigma_b = std(by);
end
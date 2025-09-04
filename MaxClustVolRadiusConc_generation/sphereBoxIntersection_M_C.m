% a MATLAB function to compute the intersection volume of a sphere centered at 
% (x,y,z) with a given radius r and a cubic box of size 200×200×200. 
% Monte Carlo integration setup



function V_intersect = sphereBoxIntersection_M_C(x, y, z, r, box_length)
    % Ensure the radius does not exceed the farthest corner distance
    r_max = farthestCornerDistance(x, y, z, box_length);
    r = min(r, r_max);
    
    % Define box boundaries
    x_min = 0; x_max = box_length;
    y_min = 0; y_max = box_length;
    z_min = 0; z_max = box_length;
    
    % Monte Carlo integration setup
    num_samples = 1e6; % Number of random points
    
    % Generate random points inside the sphere
    theta = rand(num_samples, 1) * 2 * pi;
    phi = acos(2 * rand(num_samples, 1) - 1);
    radius = (rand(num_samples, 1)).^(1/3) * r;
    
    % Convert spherical coordinates to Cartesian
    x_s = x + radius .* sin(phi) .* cos(theta);
    y_s = y + radius .* sin(phi) .* sin(theta);
    z_s = z + radius .* cos(phi);
    
    % Count points inside the box
    inside_box = (x_s >= x_min & x_s <= x_max) & ...
                 (y_s >= y_min & y_s <= y_max) & ...
                 (z_s >= z_min & z_s <= z_max);
    
    % Estimate intersection volume
    V_sphere = (4/3) * pi * r^3; % Full sphere volume
    V_intersect = V_sphere * (sum(inside_box) / num_samples);
end

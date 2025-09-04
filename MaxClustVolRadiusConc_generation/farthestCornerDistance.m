% Here’s a MATLAB function to compute the farthest distance from a point 
% (x,y,z) to the corners of the box:

function d_max = farthestCornerDistance(x, y, z, box_length)
    % Define the 8 corner coordinates of the cube (size 200)
    corners = [0, 0, 0;
               0, 0, box_length;
               0, box_length, 0;
               0, box_length, box_length;
               box_length, 0, 0;
               box_length, 0, box_length;
               box_length, box_length, 0;
               box_length, box_length, box_length];

    % Compute distances to all 8 corners
    distances = sqrt((corners(:,1) - x).^2 + (corners(:,2) - y).^2 + (corners(:,3) - z).^2);

    % Find the maximum distance
    d_max = max(distances);
end

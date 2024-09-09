
% This function calculates the centroid of a grain (w/ periodic BCs)
% and shifts it to the center of the 3D volume. This currently assumes
% a boolean field where array=1 inside the grain and array=0 outside.
function [shifted_array] = shift_grain_to_center_3D(array)

    % Determine array size
    Nx = size(array, 1);
    Ny = size(array, 1);
    Nz = size(array, 1);
    
    % Determine the center of the volume
    center_x = (Nx + 1) / 2;
    center_y = (Ny + 1) / 2;
    center_z = (Nz + 1) / 2;

    % Compute the centroid of the grain (considering periodic BCs)
    [centroid_x,centroid_y,centroid_z] = calculate_centroid_periodic_3D(array,Nx,Ny,Nz);

    % Compute the shifts required
    shift_x = center_x - centroid_x;
    shift_y = center_y - centroid_y;
    shift_z = center_z - centroid_z;
    
    % Create a shifted version of the array
    shifted_array = circshift(array, [round(shift_x), round(shift_y), round(shift_z)]);

end


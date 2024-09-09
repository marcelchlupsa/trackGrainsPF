
% Calculate the centroid of a grain/island with periodic BCs.
% Note: The grain should be a boolean defined by array=1 and
% array=0 outside of the grain.
% This code follows the algorithm by Bai & Breen (2008)
% 'Calculating Center of Mass in an Unbounded 2D Environment'
function [com_x,com_y,com_z] = calculate_centroid_periodic_3D(array,Nx,Ny,Nz)

    x_i = 0;
    y_i = 0;
    z_i = 0;

    x_j = 0;
    y_j = 0;
    z_j = 0;

    x_k = 0;
    y_k = 0;
    z_k = 0;

    for i = 1:Nx
        for j = 1:Ny
            for k = 1:Nz
            
                if array(i,j,k) == 1
                
                    % Map to theta
                    theta_i = i*2*pi / Nx;
                    r_i = Nx/(2*pi);
    
                    theta_j = j*2*pi / Ny;
                    r_j = Ny/(2*pi);
    
                    theta_k = k*2*pi / Nz;
                    r_k = Nz/(2*pi);
                    
    
                    x_i = x_i + r_i * cos(theta_i);
                    y_i = y_i + j;
                    z_i = z_i + r_i * sin(theta_i);
    
                    x_j = x_j + i;
                    y_j = y_j + r_j * cos(theta_j);
                    z_j = z_j + r_j * sin(theta_j);

                    x_k = x_k + r_k * cos(theta_k);
                    y_k = y_k + r_k * sin(theta_k);
                    z_k = z_k + k;
    
                end
            end
        end
    end

    volume = sum(array,'all');

    % Calculate x_bar, y_bar, z_bar for both tubes
    x_i = x_i / volume;
    y_i = y_i / volume;
    z_i = z_i / volume;

    x_j = x_j / volume;
    y_j = y_j / volume;
    z_j = z_j / volume;

    x_k = x_k / volume;
    y_k = y_k / volume;
    z_k = z_k / volume;

    theta_i = atan2(-z_i,-x_i) + pi;
    theta_j = atan2(-z_j,-y_j) + pi;
    theta_k = atan2(-y_k,-x_k) + pi;

    com_x = Nx * theta_i / (2*pi);
    com_y = Ny * theta_j / (2*pi);
    com_z = Nz * theta_k / (2*pi);
end

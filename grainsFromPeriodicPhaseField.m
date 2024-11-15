function [gid_map,buffer] = grainsFromPeriodicPhaseField(OP)
    
    % next steps: periodicity along only selected dimensions, handling 2d
    % maps, save numElement with the order parameter
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % must handle sptensor
    if isa(OP,'sptensor')
        OP = double(full(OP)); % lost of memory here (full function)
        % OP = single(double(full(OP))); % to save memory
    end
    if ~isa(OP,'double')
        warning('expected double type input\n')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % assume a 3d grain map
    if ndims(OP)==4
        [~,OP] = max(OP,[],4);
    elseif ndims(OP)~=3
        error('OP must have either 3 or 4 dims\n')
    end
    [l,m,n] = size(OP);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform initial segmentation
    fprintf('initial segmentation\n')
    gid_map = segmentMap(OP);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for each dimension, find how many layers are needed to complete the
    % grains along the periodic boundary
    fprintf('accounting for periodic boundary conditions\n')
    
    N1 = layersToAppend(gid_map,1);
    N2 = layersToAppend(gid_map,2);
    N3 = layersToAppend(gid_map,3);
    buffer = [N1 N2 N3];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % append layers, rerun segmentation
    OP = cat(1,OP,OP(1:N1,:,:));
    OP = cat(2,OP,OP(:,1:N2,:));
    OP = cat(3,OP,OP(:,:,1:N3));
    
    fprintf('improved segmentation\n')
    gid_map = segmentMap(OP);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % filter edge grains, filter inside repeat grains
    removeid = [unique(gid_map(1,:,:)); unique(gid_map(:,1,:)); unique(gid_map(:,:,1)); ...
        unique(gid_map(end,:,:)); unique(gid_map(:,end,:)); unique(gid_map(:,:,end))];
    gid_map(ismember(gid_map,removeid)) = 0;
    [gct,grp] = groupcounts(reshape(gid_map,[numel(gid_map),1]));
    
    masked = gid_map;
    masked(1:l,1:m,1:n) = 0;
    [gct2,grp2] = groupcounts(reshape(masked,[numel(masked),1]));
    removeid = grp(ismember(cat(2,grp,gct),cat(2,grp2,gct2),'rows'));
    gid_map(ismember(gid_map,removeid)) = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function [newmap] = segmentMap(map)

        unique_vals = unique(map(map > 0));  % Unique positive order parameters
        newmap = zeros(size(map));  % Output array for unique grain IDs
        num_vals = numel(unique_vals);
        CC_results = cell(num_vals, 1);  % Store results from bwconncomp

        % Parallel processing for connected components
        parfor i = 1:num_vals
            val = unique_vals(i);
            mask = (map == val);  % Create binary mask for the current order parameter
            CC_results{i} = bwconncomp(mask);  % Compute connected components
        end

        % Compute offsets for unique grain IDs
        offsets = zeros(num_vals, 1);
        for i = 2:num_vals
            offsets(i) = offsets(i-1) + CC_results{i-1}.NumObjects;
        end

        total_indices = sum(cellfun(@(cc) sum(cellfun(@numel, cc.PixelIdxList)), CC_results));
        all_indices = zeros(total_indices, 1);
        all_values = zeros(total_indices, 1);
        % Fill preallocated arrays
        current_idx = 1;
        for i = 1:num_vals
            % fprintf('%d/%d\n',i,num_vals)
            CC = CC_results{i};
            base_id = offsets(i);
        
            for j = 1:CC.NumObjects
                num_voxels = numel(CC.PixelIdxList{j});
                all_indices(current_idx:current_idx + num_voxels - 1) = CC.PixelIdxList{j};
                all_values(current_idx:current_idx + num_voxels - 1) = base_id + j;
                current_idx = current_idx + num_voxels;
            end
        end
        
        
        % Assign grain IDs in one operation
        newmap(all_indices) = all_values;

    end

    function N = layersToAppend(map, dim)
        % Extract grain IDs from the first layer along the specified dimension
        switch dim
            case 1
                gids1 = unique(map(1, :, :)); % Unique grain IDs in the first layer along dimension 1
            case 2
                gids1 = unique(map(:, 1, :)); % Unique grain IDs in the first layer along dimension 2
            case 3
                gids1 = unique(map(:, :, 1)); % Unique grain IDs in the first layer along dimension 3
            otherwise
                error('Invalid dimension specified. Choose 1, 2, or 3.');
        end
        
        % Remove background (grain ID 0) from the list
        gids1(gids1 == 0) = [];
        
        % Initialize variables
        N = 2; % Start checking from the second layer
        condition = true; % Condition to control the loop
        
        % Loop until we find a layer with no grain IDs from the first layer
        while condition
            % Extract grain IDs from the N-th layer along the specified dimension
            switch dim
                case 1
                    gidsN = unique(map(N, :, :)); % Grain IDs in layer N along dimension 1
                case 2
                    gidsN = unique(map(:, N, :)); % Grain IDs in layer N along dimension 2
                case 3
                    gidsN = unique(map(:, :, N)); % Grain IDs in layer N along dimension 3
            end
            
            % Remove background (grain ID 0) from the list
            gidsN(gidsN == 0) = [];
            
            % Check if any grain ID from the first layer is present in the N-th layer
            condition = any(ismember(gids1, gidsN));
            
            % Increment N to check the next layer
            if condition
                N = N + 1;
                
                % Stop if N exceeds the size of the map along the specified dimension
                if N > size(map, dim)
                    error('No layer found that satisfies the condition.');
                end
            end
        end
    end

end
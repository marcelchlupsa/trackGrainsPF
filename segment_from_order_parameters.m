function [map_out, numElement, varargout] = segment_from_order_parameters(map_in, boundsmask, varargin)
% Marcel Chlupsa

% v0.1(1/19/2024) split order parameters into unique grain ids
% v0.2(2/14/2024) add a size threshold varargin
% v0.3(3/20/2024) add option for 4D input, take the order param with max
% value
% v0.4(5/02/2024) add numElement as output, with column 3 for parent OP
% v0.5(6/11/2024) if an order parameter equals one in the first corner 
% voxel, assume that order parameter is air - it doesn't contain grains.

% INPUTS
%   map_in
%   boundsmask
%   varargin(=size_thresh)

% OUTPUTS
%   map_out
%   numElement
%   varargout(=map_out_size_filt), 

%%
if ndims(map_in)==4 % grayscale map of each OP, dim-4 is OPs

    % make a new map that stores the largest op for each voxel
    [~,map_in] = max(map_in,[],4);

end

if ndims(map_in)==3 % simplified volume w OP assigned for each voxel
    
    % mandate air in the first corner
    if map_in(1)~=0
        map_in(map_in==map_in(1)) = 0;
    end

    if size(map_in,1)==size(boundsmask,1) && size(map_in,2)==size(boundsmask,2) && size(map_in,3)==size(boundsmask,3)
    
        order_params = unique(map_in);
        
        gid = 1;

        map_out = zeros(size(map_in));
        for i = 1:length(order_params)
            BW = map_in==order_params(i);
            BW(boundsmask==0) = 0;
            CC = bwconncomp(BW,26);
            for j = 1:length(CC.PixelIdxList)
                map_tmp = zeros(size(map_in));
                map_tmp(CC.PixelIdxList{j}) = gid;
                map_out = map_out + map_tmp;

                numElement(gid,1) = gid;
                numElement(gid,2) = numel(CC.PixelIdxList{j});
                numElement(gid,3) = i;

                gid = gid + 1;
            end
        end
        fprintf('result contains %d grains...\n',length(unique(map_out))-1)
        
        if ~isempty(varargin)
            % apply size filter and save as a secondary output map
            map_tmp = map_out;

            rmIDs = numElement(numElement(:,2) <= varargin{1}, 1);
            
            numElement_tmp = numElement(numElement(:,2)>varargin{1},:);

            map_tmp(ismember(map_tmp,rmIDs)) = 0;

            varargout{1} = map_tmp;
            varargout{2} = numElement_tmp;
    
            fprintf('size-filtered result contains %d grains...\n',length(unique(map_tmp))-1)
        end
    
    else
        error('dimensions of map_in and boundsmask do not match!')
    
    end

else
    warning('incorrect number of array dimensions for map_in (must be either 3D or 4D)')

end
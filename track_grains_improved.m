%% Grain tracking PF - no orientation - improved using neighborhood information - Varun

function [tracked_grains] = track_grains_improved(volume1, volume2)

% INPUTS
%   volume1 - grain ID map from older state
%   volume2 - grain ID map from new state

% OUTPUTS
%   tracked_grains = Nx2 array [older grain IDs new grain IDs]

%Get unique grain IDs
labels2 = unique(volume2(volume2 > 0));

%Create output array - second state should have fewer grains in total
tracked_grains = zeros(length(labels2),2);

%Iterate and track based on volume overlap
for i = 1:length(labels2)
    trckd_gid = mode(volume1(volume2==labels2(i)));
    tracked_grains(i,1) =  trckd_gid;
    tracked_grains(i,2) =  labels2(i);   
end

disp('Initial matching complete - Trying to improve...\n')


%Improve tracking - Check for redundancy in matches

% Find the unique values and their indices in the array
[uniqueValues, ~, indices] = unique(tracked_grains(:,1));

% Count the occurrences of each unique value
counts = accumarray(indices, 1);

% Find the values that occur more than once
repeatedIndices = counts > 1;
repeatedValues = uniqueValues(repeatedIndices);

% Check grains that have been matched to these values
toMatchgID2 = [];
for i = 1:length(repeatedValues)
    repeatedgIDs2 = tracked_grains(tracked_grains(:,1)==repeatedValues(i), 2);
    checkfrac = 0;
    bestmatch = 0;
    for j = 1:length(repeatedgIDs2)
        occupiedFrac = sum(volume1(volume2==labels2(i)) == repeatedgIDs2(j))/sum(volume1(:) == repeatedValues(i));
        if occupiedFrac>checkfrac
            checkfrac = occupiedFrac;
            bestmatch = repeatedgIDs2(i);
        end
    end
    toMatchloc = repeatedgIDs2(repeatedgIDs2~=bestmatch);
    toMatchgID2 = [toMatchgID2; toMatchloc];
    tracked_grains(tracked_grains(:,1)==repeatedValues(i), :) = [];
    tracked_grains = [tracked_grains; [repeatedValues(i) bestmatch]];
end

disp('Removed redundant matches - only best possible matches remain \n')


%Try to match these untracked grains based on vol overlap
for i = 1:length(toMatchgID2)
    possible_matches = volume1(volume2==toMatchgID2(i));
    while ~isempty(possible_matches)
        gidmatch_maybe = mode(possible_matches);
        if ismember(gidmatch_maybe, tracked_grains(:,1))
            possible_matches(possible_matches==gidmatch_maybe) = [];
        else
            tracked_grains = [tracked_grains; [gidmatch_maybe toMatchgID2(i)]];
            break;
        end
    end
end

%Print
fprintf('Tracking Complete - Fraction of grains tracked: %.2f \n', length(tracked_grains)/length(labels2));

end


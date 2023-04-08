%==========================================================================
% Get longitudinal position along the pathway, starting from the tree inlet 
function [x_loc] = x_location_GrapphTree(layer,G,segm,Length)
%
%update for Strahler ordering tree
%       Updated by Vasilina, October 30, 2020
%==========================================================================
nlayer = max(layer); % number of layers
snode = segm(:,1); %start node of the segment
enode = segm(:,2);
nsegm = size(segm,1);

x_loc = zeros(nsegm,1);

x_loc(1) = Length(1); %0?
% forward loop, layers in reverse order with generations
for j=nlayer-1:-1:1 % start from the second generation
    layernnodes = find(layer(:)==j); % end nodes at jth layer
    
    for i=1:size(layernnodes,1)
        si = layernnodes(i); % start node
        % get parent segment id
        parsegid = find(enode(:)==si); 
        
        % go over the segments
        if outdegree(G,si)>0
            % get 2 daughter-segments id
            segid = find(snode(:)==si);
            
            for ii=1:2 % daughter segments
                dsegid = segid(ii);
                x_loc(dsegid) = Length(dsegid) + x_loc(parsegid);
            end
        end
    end
end


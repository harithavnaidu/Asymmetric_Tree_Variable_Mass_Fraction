function [Graph,nsegm,snode,enode,nlayer,layer,Radius,Length] =...
         readGraphInitialData
    % Load and plot graph data
    segm = load("GraphSegmentsNodes.inp");
    snode = segm(:,1); %start node of the segment
    enode = segm(:,2); %end node of the segment
    order = segm(:,3); %order
    weights = order;

    maxord = max(order);
    minord = min(order);
    norder = maxord - minord + 1;%order(1); % first segment is of the last order
    order(:) = order(:) - minord + 1; %adjust ordering
    nsegm = size(segm,1); % number of segments

    % get directed graph
    Graph = digraph(snode,enode,weights);
    figure(1);
    plotGraph = plot(Graph,'Layout','layered');
    % node <-> layer vector (#1 node <-> 6th layer)
    layer = plotGraph.YData; % each node is assigned to the layer number
    nlayer = max(layer); % number of layers
    % get outdegree:1 for inlet; 0 for outlets; 2 for other nodes
    Gdegree = outdegree(Graph); % for node numbers

    PostprocessingGraph;
    
    % Load diameters and lenghts means vs Strahler orders (in cm)
    vesselSizeData = load("StrahlerOrderDandL.inp");
    %sort on ascending order
    vesselSizeData = flip(vesselSizeData);
    % Load diameters and lenghts vs Strahler orders
    % in cm
    Rorder = 0.5.*vesselSizeData(:,2);
    Lorder = vesselSizeData(:,3);
    % %increase the length twice
    % Lorder = Lorder*2;

    %set raidus and length arrays
    Radius = zeros(nsegm,1);
    Length = zeros(nsegm,1);

    for  i = 1:nsegm
        Radius(i) = Rorder(order(i));
        Length(i) = Lorder(order(i));
    end
end
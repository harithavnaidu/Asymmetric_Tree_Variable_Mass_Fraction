% %% ========== BUILD GRAPH TREE USING STRAHLER ORDERING DATA ===============
% - Reading the connectivity matrix and getting a tree constructed
% - Algorithm and test example from [Yang,...Marsden,Feinstein(BMM-2019)]
%   for unique reconstruction of a tree
% 
% Assumption:
% - 1 or 2 elements of highest order
% - to compute # of element per order we rounding the summ to avoid
%   noninteger values
% - C(k,k)<1 to avoid too much of self branching (also consistent with 
%   Huang data)
% - Not all elements outlets are of the same order. May need to cut or make
% an assumptions to merge lower orders to the same order?
%
% Given: 
%   -Connectivity matrix (connMatrix.inp)
%   -Order elements lengths and diameters
% 
% Output:
%   -AllSegmList.dat
%   -Tree.fig
%
% Code by Vasilina Filonova
%   - created on June 18, 2019
%   - revised on July 15, 2019
%%==========================================================================
clear

%% Read Connectivity Matrix
% Test
% connMatrix = [0, 2.4, 1;
%               0, 0,   3;
%               0, 0,   0];
% % fileID = fopen("connMatrix.dat",'w');
% % fmt = [repmat('%6.3f ', 1, size(connMatrix,2)-1),'%6.3f\n'];
% % fprintf(fileID,fmt,connMatrix.');
% % fclose(fileID);

load('connMatrix.inp')
maxOrder = size(connMatrix,1);

%% Key number sets for a tree

% Get a column for an average number of all children element of order m 
% for each element of order n
numChildPerElInOrder = zeros(maxOrder,1);
for n = 1:maxOrder
    numChildPerElInOrder(n) = sum(connMatrix(:,n),1);
end
a = numChildPerElInOrder;
numChildPerElInOrder = floor(a);
remNumChildElPerElInOrder = a - numChildPerElInOrder;

% Get number of elements of each order 
% compute backward, assume ConnMatrix(m,m)=0 always
numElPerOrder = zeros(maxOrder,1);
% numElPerOrder(maxOrder) = 1; % assumption
% for m = maxOrder-1 :-1:1% if C(m,m)==0
%      numElPerOrder(m) = floor( connMatrix(m,m+1:maxOrder)*numElPerOrder(m+1:maxOrder));
% end

% %Get number of element taking into accout the reminders
numElPerOrder(maxOrder) = floor(1/(1-connMatrix(maxOrder,maxOrder))); %my estimation
% C(m,m)<1 by definition, then extend the formula for C(k,k)>0
for m = maxOrder-1 :-1:1
    numElPerOrder(m) = connMatrix(m,m+1:maxOrder)*numElPerOrder(m+1:maxOrder);
    numElPerOrder(m) = numElPerOrder(m)/(1-connMatrix(m,m));
end
% round
numElPerOrder = floor(numElPerOrder);

% number of all elements
NumAllElem = sum(numElPerOrder,1); % number of all elements

% number of all segments
NumAllSegm = 0;
for m = maxOrder:-1:1
    if numChildPerElInOrder(m)== 0
        NumAllSegm = NumAllSegm + numElPerOrder(m);
    elseif numChildPerElInOrder(m)== 1
        NumAllSegm = NumAllSegm + 2*numElPerOrder(m);
    elseif numChildPerElInOrder(m)>1
        NumAllSegm = NumAllSegm + ...
            numElPerOrder(m)*(numChildPerElInOrder(m)-1);
    end
end

%% Find a position of the last nonzero element in a column
maxk = zeros(maxOrder,1);
for n=1:maxOrder
    a = find(connMatrix(:,n),1,'last');
    if (a) maxk(n) = a; end
end

%% ---------------------- Get parent element subgraph---------------------
% assume C(:,maxOrder) are whole numbers
% k = maxk(maxOrder);

[NewSegmList,nNodes]= GetRootSubgraph...
            (maxOrder,connMatrix,numChildPerElInOrder,numElPerOrder);

% Initialize list of all segments nodes and orders
AllSegmList = NewSegmList;
% remove rows with all zeros from AllSegmList
AllSegmList( ~any(AllSegmList,2), : ) = [];  
% remove duplicated rows if needed
AllSegmList = unique(AllSegmList,'rows');

% fileID = fopen("AllSegmList-Root.dat",'w');
% fprintf(fileID,'%5d %5d %5d\n',AllSegmList');
% fclose(fileID);

nSegm = size(AllSegmList,1);

% Construct Graph
starts = AllSegmList(:,1);  % column of start nodes
ends = AllSegmList(:,2);    % column of end nodes
weights = AllSegmList(:,3); % order labels

G = graph(starts,ends,weights);

% test: asign edge radius as a function of order
LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);

%Plot graph
fig=figure(1);
% plot(G,'Marker','none','NodeLabel',{});
p = plot(G,'EdgeLabel',G.Edges.Weight,'LineWidth',LWidths,'MarkerSize',0.01);
layout(p,'layered')
% layout(p,'force')
title('Root subgraph');
colormap(jet)
p.EdgeCData=G.Edges.Weight;
c = colorbar; c.Location='southoutside';
c.Ticks = [12 13 14 15];
% c.TickLabels = {'12','13','14','15'};
savefig(fig,'Root.fig');
% layout(p,'force','WeightEffect','direct')
% close(1);

%% -------------------- Get subgraph of all other children branches --------
truncationOrder=11
for m = maxOrder-1: -1: truncationOrder
    % largest child element order
    k = maxk(m); 
    
    % check existence of a child branch
    if k>0
        %run subroutine to udpate NewSegmList and nNodes
        [AllSegmList,NewSegmList,NewNNodes] = GetChildSubgraph...
                      (m,nNodes,AllSegmList,connMatrix,...
                      numChildPerElInOrder,numElPerOrder);  
        %update   
        nSegm = size(AllSegmList,1);
        AllSegmList(nSegm+1: nSegm+size(NewSegmList,1),:) = NewSegmList;
%         nSegm = size(AllSegmList,1);
        nNodes = NewNNodes;
    end
    
%     % remove rows with all zeros from AllSegmList
%     AllSegmList( ~any(AllSegmList,2), : ) = [];  
%     % remove duplicated rows if needed
%     AllSegmList = unique(AllSegmList,'rows');
%     
%     % plot
%     GraphTree = graph(AllSegmList(:,1),AllSegmList(:,2),AllSegmList(:,3));
%     LWidthsG = 5*GraphTree.Edges.Weight/max(GraphTree.Edges.Weight);
%     figure(m+10)
%     plot(GraphTree,'EdgeLabel',GraphTree.Edges.Weight,'LineWidth',LWidthsG);
end

% remove rows with all zeros from AllSegmList
AllSegmList( ~any(AllSegmList,2), : ) = [];  

% remove duplicated rows if needed
AllSegmList = unique(AllSegmList,'rows');
% 
% sort on order
AllSegmList = sortrows(AllSegmList,3,'descend');

% % truncate by order 9
%  C = AllSegmList(AllSegmList(:,3)>9,:);

%%
% keyboard;
% output list
fileID = fopen("AllSegmList.dat",'w');
fprintf(fileID,'%5d %5d %5d\n',AllSegmList');
fclose(fileID);

% keyboard;
% %% --------------------- Plot Graph ---------------------------------------
starts = AllSegmList(:,1);  % column of start nodes
ends = AllSegmList(:,2);    % column of end nodes
weights = AllSegmList(:,3); % order labels

GraphTree = graph(starts,ends,weights);

% keyboard;
% fig = figure(2);
% plot(GraphTree,'Marker','none','NodeLabel',{})
% savefig(fig,'Network.fig');

% test: asign edge radius as a function of order
LWidthsG = 2*GraphTree.Edges.Weight/max(GraphTree.Edges.Weight);

% plot
fig = figure(2);
% p = plot(G2,'LineWidth',LWidths2);
% p = plot(GraphTree,'EdgeLabel',GraphTree.Edges.Weight,'LineWidth',LWidthsG);
% savefig(fig,'Tree.fig');
% p = plot(GraphTree,'EdgeLabel',GraphTree.Edges.Weight,'LineWidth',LWidthsG);
p = plot(GraphTree,'LineWidth',LWidthsG,'MarkerSize',0.01);
% layout(p,'force')
layout(p,'layered')
% title('Labeled Edge Orders');
colormap(jet)
p.EdgeCData=GraphTree.Edges.Weight;
c = colorbar; c.Location='southoutside';
% c.Ticks = [9 10 11 12 13 14 15];
% c.Ticks = [6 7 8 9 10 11 12 13 14 15];
savefig(fig,'Tree.fig');
% % close(fig);
% 
% %% --------------------- Graph manipulations ------------------------------
% % length?,%maximum flow?
% 
% get adjacency matrix
% tree_adj = adjacency(GraphTree);

% 
% % add edge label = number of the edge
% figure(3)
% label = 1:1:size(AllSegmList,1);
% p = plot(GraphTree,'EdgeLabel',GraphTree.Edges.Weight,'LineWidth',LWidthsG);
% title('Labeled Edge Numbering');
% labeledge(p,starts,ends,label)
% 
% % directed graph
% figure(4)
% diGraphTree = digraph(starts,ends,weights);
% LWidthsDiG = 5*diGraphTree.Edges.Weight/max(diGraphTree.Edges.Weight);
% plot(diGraphTree,'LineWidth',LWidthsDiG);
% title('Directed Graph');
% 
% % outflow degree of nodes
% outdgr = outdegree(diGraphTree);
% 
% % Shortest path distances of all node pairs
% shdistG = distances(GraphTree);
% shdistDG = distances(diGraphTree);
% 
% %% Get steady flow





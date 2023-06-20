%get directed graph
diGraphTree = digraph(starts,ends,weights);
a=outdegree(diGraphTree)==0; % logical 1-0 for all nodes
tOut=find(outdegree(diGraphTree)==0);% get end nodes ID
% sOut=predecessors(diGraphTree,tOut);

%find maximum order of outlet nodes
% segOut = findedge(diGraphTree,sOut,tOut);
for i=1:size(tOut)
    [segOut(i),sOut(i)]=inedges(diGraphTree,tOut(i));
end

wOut=diGraphTree.Edges{findedge(diGraphTree,sOut,tOut), 'Weight'};
maxW = max(wOut); 
minW = min(wOut); 

% assume number of weigths maxW
% all weights below maxW are substituted by maxW
newWeights=weights;
newWeights(newWeights<maxW)=maxW;

newG = graph(starts,ends,newWeights);
LWidthsG = 2*newG.Edges.Weight/max(newG.Edges.Weight);

fig=figure(3);
newp = plot(newG,'LineWidth',LWidthsG);
layout(newp,'layered')
title('Graph with outlets of the same order');
colormap(jet)
newp.EdgeCData=newG.Edges.Weight;
c = colorbar; c.Location='southoutside';
c.Ticks = [11 12 13 14 15];

%export data
temp=AllSegmList;
temp(:,3) = newWeights;
fileID = fopen("AllSegmList_SmalOrderOutlets.dat",'w');
fprintf(fileID,'%5d %5d %5d\n',temp');
fclose(fileID);




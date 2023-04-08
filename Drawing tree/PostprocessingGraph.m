a=Gdegree==0; % logical 1-0 for all nodes
tOut=find(Gdegree==0);% get terminal nodes ID

%find maximum order of outlet nodes
for i=1:size(tOut)
    [segOut(i),sOut(i)]=inedges(Graph,tOut(i));
end

wOut=Graph.Edges{findedge(Graph,sOut,tOut), 'Weight'};
maxW = max(wOut); 
minW = min(wOut); 

% assume number of weigths maxW
% all weights below maxW are substituted by maxW
newWeights=weights;
newWeights(newWeights<maxW)=maxW;

newG = graph(snode,enode,newWeights);
LWidthsG = 2*newG.Edges.Weight/max(newG.Edges.Weight);

figure;
newp = plot(newG,'LineWidth',LWidthsG);
layout(newp,'layered')
%title('Graph with outlets of the same order');
colormap(jet)
newp.EdgeCData=newG.Edges.Weight;
c = colorbar; c.Location='southoutside';
c.Ticks = [11 12 13 14 15];




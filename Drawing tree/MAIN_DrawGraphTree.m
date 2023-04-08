% created by Vasilina Filonova: 01-04-2021
% draw a tree, using graphs or fractals and L,D data
% angle is prescribed based on better visualization, no physical meaning
    
clc;
clear all;
close all;
format long;
% workspace

global graphFlag fractalFlag;
graphFlag = 1; 
fractalFlag = 0;

if graphFlag == 1
    %read data for graph structure from GraphSegmentsNodes.inp
    %read R&L from StrahlerOrderDandL.inp
    %show the graph
    [Graph,nsegm,snode,enode,nlayer,layer,Radius,Length] = ...
    readGraphInitialData;
    
    %draw initial tree
%     drawTree(Graph,nsegm,snode,enode,nlayer,layer,Radius,Length);
    
    %get homeostaticly optimized tree
    % get R and L (in m) vs segment number
    vesselSizeHomeostaticData = load("HomeostaticRandL.dat");
    Radius = 100.*vesselSizeHomeostaticData(:,1); %in cm
    Length = 100.*vesselSizeHomeostaticData(:,2); %in cm
    VesselID = vesselSizeHomeostaticData(:,3);

    drawTree(Graph,nsegm,snode,enode,nlayer,layer,Radius,Length,VesselID);
elseif fractalFlag == 1
    
    %create a fractal tree here
    %to do: read data from the file
    [Ngen,Radius,Length,ID] = generateFractalTree;
    Radius = 100.*Radius; %m to cm
    Length = 100.*Length; %m to cm
    
    %draw the initial fractal tree
    drawFractal(Ngen,Radius,Length,ID)    
end




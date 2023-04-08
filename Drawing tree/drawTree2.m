function drawTree2(Graph,nsegm,snode,enode,nlayer,layer,Radius,Length)

angle = pi/3;
factor = 1/4;

% initial vessel parameters
L = Length(1); %width
R = Radius(1); %height
%starting position:
x = 0; y = 0;

%RECTANGLE('Position', [x y w h])
figure
xvertices =[x x   x+L x+L x   x];
yvertices =[y y-R y-R y+R y+R y];
plot(xvertices,yvertices,'r'); %,'LineWidth',2);
%f = fill(yvertices,-xvertices,'r');
axis equal; %set x/y aspect ratio to 1
%set(f,'EdgeColor','k');
%grid on;
% set(gca,'YDir','reverse');
set(gca,'CameraUpVector',[-1 0 0]);
xlabel('cm'); ylabel('cm');
hold on;
% axis([0 60 -40 20])

segEnd = zeros(nsegm,2);
segEnd(1,:) = [L 0];
% forward loop, layers in reverse order with generations
tetaSegm = zeros(nsegm,1); % angle the segment was rotated to be plot
layerL = zeros(nlayer,1);
layerL(nlayer) = L; % 0; %L;
for j=nlayer-1:-1:1 % start from the second generation
    layernnodes = find(layer(:)==j); % end nodes at jth layer
    
%     %get segments for these nodes
%     layersegm = zeros(size(layernnodes,1),1);
%     for k=1:size(layernnodes,1)
%         layersegm(k) = find(enode(:)==layernnodes(k));
%     end
%     
%     %add min L if current gen to the L of the previous gen
%     layerL(j) = layerL(j+1) + min(Length(layersegm));
    
    for i=1:size(layernnodes,1)
        si = layernnodes(i); % start node
        %get parent segment id
        parsegid = find(enode(:)==si); 
        
        % go over the segments
        if outdegree(Graph,si)>0
            % get 2 daughter-segments id
            segid = find(snode(:)==si);
            
            layerL(j) = layerL(j+1) + min(Length(segid));
            
            randomindex = randi(2);
            for ii=1:2 % daughter segments
                dsegid = segid(ii);
                L = Length(dsegid);
                R = Radius(dsegid);
                x = 0; y = 0; %relative starting position
                xvertices =[x x   x+L x+L x   x];
                yvertices =[y y-R y-R y+R y+R y];
                
                %get angle
                if (layerL(j+1) + L > layerL(j))
                    %use law of cosines
                    cosAngle = (layerL(j+1)^2 + L^2 - layerL(j)^2)/...
                               (2*layerL(j+1)*L);
                    delta = pi - acos(cosAngle);
                else
                    delta = angle; %0;
                end   
                % chose branch
                %random branch
                
                if (ii==1)
%                 if (ii == randomindex)
                    teta = tetaSegm(parsegid) - factor*delta; 
                else %ii==2
                    teta = tetaSegm(parsegid) + factor*delta;
                end
                tetaSegm(dsegid) = teta;
                
                % rotate and translate 
                vertices(1,:) = xvertices; 
                vertices(2,:) = yvertices;
                Rotation(1:2,1:2) = [cos(teta) -sin(teta);sin(teta) cos(teta)];
                rotatedVertices = Rotation * vertices;
                
                %get center of the parent vessel joint side
                x = segEnd(parsegid,1);
                y = segEnd(parsegid,2);
                %save center of the current vessel joint side
                segEnd(dsegid,:) = transpose([x; y] + Rotation * [L; 0]);
                
                %plot next vessel
                plot(x + rotatedVertices(1,:),...
                     y + rotatedVertices(2,:),'r');
                set(gca,'CameraUpVector',[-1 0 0]);
                hold on;                
            end
        end
    end
end
hold off

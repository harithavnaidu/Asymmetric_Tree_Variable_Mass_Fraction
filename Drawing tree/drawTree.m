function drawTree(Graph,nsegm,snode,enode,nlayer,layer,Radius,Length,VesselID)

angle = pi/4;
factor = 1/2;

% initial vessel parameters
L = Length(1); %width
R = Radius(1); %height
%starting position:
x = 0; y = 0;

%RECTANGLE('Position', [x y w h])
figure;
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

segEnd = zeros(nsegm,2);
segEnd(1,:) = [L 0];
% forward loop, layers in reverse order with generations
tetaSegm = zeros(nsegm,1); % angle the segment was rotated to be plot
for j=nlayer-1:-1:1 % start from the second generation
    layernnodes = find(layer(:)==j); % end nodes at jth layer
    
    for i=1:size(layernnodes,1)
        si = layernnodes(i); % start node
        %get parent segment id
        parsegid = find(enode(:)==si);
        
        % go over the segments
        if outdegree(Graph,si)>0
            % get 2 daughter-segments id
            segid = find(snode(:)==si);
            
            randomindex = randi(2);
            for ii=1:2 % daughter segments
                dsegid = segid(ii);
                L = Length(dsegid);
                R = Radius(dsegid);
                ID = VesselID(dsegid);
                x = 0; y = 0; %relative starting position
                xvertices =[x x   x+L x+L x   x];
                yvertices =[y y-R y-R y+R y+R y];
                if ID==1
                    clr = [62 88 166]/255;
%                     randomindex = 2;
                elseif ID==2
                    clr = [37 158 49]/255;
%                     randomindex = 1;
                else
                    clr = [238 31 35]/255;
                end
                % chose branch
                %                 if (ii==1)
                if (ii == randomindex)
                    %                     teta = angle - fraction*angle*(j+1);
                    %                     teta = angle * (1 - fraction*(j+1));
                    teta = tetaSegm(parsegid) - ...
                        angle*factor*j/nlayer; %* (1 - fraction*(j+1));
                    %                     if (teta <- pi/2)
                    %                        teta = -tetaSegm(parsegid); %- angle; %-tetaSegm(parsegid);
                    %                     end
                else %ii==2
                    teta = tetaSegm(parsegid) + ...
                        angle*factor*j/nlayer; % + fraction*angle*(j+1);
                    %                     if (teta > pi/2)
                    %                         teta = tetaSegm(parsegid);%angle; %tetaSegm(parsegid);
                    %                     end
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
                if ID==0
                    plot(x + rotatedVertices(1,:),...
                        y + rotatedVertices(2,:),'color',clr);
                    set(gca,'CameraUpVector',[-1 0 0]);
                    hold on;
                else
                    a = fill(x + rotatedVertices(1,2:5),...
                        y + rotatedVertices(2,2:5),clr);
                    set(gca,'CameraUpVector',[-1 0 0]);
                    uistack(a,'top') 
                    hold on;
                end
            end
        end
    end
end
hold off
box off
axis off
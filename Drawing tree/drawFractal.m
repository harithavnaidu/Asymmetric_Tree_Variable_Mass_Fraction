% draw a fractal tree
% by Vasilina, 08-10-2018   

function drawFractal(Ngen,Radius,Length,ID)
    angle = pi/6;
    factor = Ngen/2;

    % initial vessel parameters
    L = Length(1,1); %width
    R = 2*Radius(1,1); %height
    x = 0; y = 0; 
    
    %RECTANGLE('Position', [x y w h])
    figure;
    xvertices =[x x   x+L x+L x   x];
    yvertices =[y y-R y-R y+R y+R y];
    plot(xvertices,yvertices,'r'); %,'LineWidth',2);
%     f = fill(yv,-xv,'r');
%     set(f,'EdgeColor','r');
    set(gca,'CameraUpVector',[-1 0 0]);
    xlabel('cm'); ylabel('cm');
    axis equal;
    hold on;

    segEnd = zeros(Ngen,2^(Ngen-1),2);
    segEnd(1,1,:) = [Length(1,1) 0];
    tetaSegm = zeros(Ngen,2^(Ngen-1),1);
    for k=2:Ngen
        for s=1:2^(k-1)
            i = ID(k,s); 
            j = k + 1 - i;
            % vessel geometry (i,j)
            L = Length(i,j);
            R = Radius(i,j);
            x = 0; y = 0;
            xvertices =[x x   x+L x+L x   x];
            yvertices =[y y-R y-R y+R y+R y];

            %parent
            sp = round(s/2);
            ip = ID(k-1,sp);
            jp = k - ip;

            % chose branch
            if (mod(s,2)==1)
                teta = tetaSegm(ip,jp,1) - ... %angle;% ...
                angle ;%* factor * k/Ngen;
                %angle - fraction*angle*(k-1);
            else
                teta = tetaSegm(ip,jp,1) + ... %angle; %...
                angle ;%* factor * k/Ngen;
                %-angle + fraction*angle*(k-1);
            end
            tetaSegm(i,j,1) = teta;
            
            % rotate and translate 
            vertices(1,:) = xvertices; 
            vertices(2,:) = yvertices;
            Rotation(1:2,1:2) = [cos(teta) -sin(teta);sin(teta) cos(teta)];
            rotatedVertices = Rotation * vertices;

            %get center of the parent vessel joint side
            x = segEnd(k-1,sp,1); %? segEnd(ip,jp,1)
            y = segEnd(k-1,sp,2);
            %save center of the current vessel joint side
            segEnd(k,s,:)= transpose([x; y] + Rotation * [L; 0]);

            % plot next vessel
            plot(x + rotatedVertices(1,:),...
                 y + rotatedVertices(2,:),'r');
            set(gca,'CameraUpVector',[-1 0 0]);
            hold on; 
        end
    end
    hold off
end

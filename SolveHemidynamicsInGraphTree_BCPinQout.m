%==========================================================================
% Solve algebraic system of equations for steady flow and term. pressure
% at all branches in a tree at once
function Table = SolveHemidynamicsInGraphTree_BCPinQout...
                                (G,segm,Qinp,Pinp,Table)
% BC: root terminal pressure at inlet, uniform flow at outlet terminals
% Table(NumberOfVessels,4)=[order, flow, terminal pressure, resistance];
% -------
% Code by Vasilina, December, 2020
%==========================================================================
%%
snode = segm(:,1); %start node of the segment
enode = segm(:,2); %end node of the segment
%maxord = max(order); minord = min(order);
%norder = maxord - minord + 1;
%order(:) = order(:) - minord + 1;
nsegm = size(segm,1); % number of segments
nnode = max(enode);

% get outdegree
% 1 for inlet; 0 for outlets; 2 for other nodes
Gdegree = outdegree(G); % for node numbers

%% Create system of equations
% Algorithm of computing flow
% P(end node) - terminal pressure at a segment
% Q(end node) - flow
% One block for one bifurcation with 3 segments and 4 nodes

% A - matrix for unknown pairs [qi,pi], ith segment
A = zeros(2*nsegm);
B = zeros(2*nsegm,1);
Res = zeros(nsegm,1);

% counter for the matrix row
cr = 0;

% inlet pressure BC
cr = cr + 1;
ep = enode(1);
%Res(1) = Table(1,4);
Res(1) = Table(1,4)/10^6/133.32; %test, Pa/m^3 -> mmHg/cm^3
% Res(1)Q(1)+ PT(1) = Pinp
A(cr,2*ep-1) = Res(1);
A(cr,2*ep) = 1; % P terminal at 1st segm
%B(cr) = Pinp;
B(cr)   = Pinp/133.32; %test Pa->mmHg

% pressure-flow relations 
for i =2:nsegm % 1st segm p-BC
    cr = cr + 1;
    si = snode(i); 
    ei = enode(i);
    %Res(i) = Table(i,4);
    Res(i)=Table(i,4)/10^6/133.32; %test, Pa*s/m^3 -> mmHg*s/cm^3
    
    % get coefficients for pressure-flow equation
    A(cr,ei*2)  = 1; % Pout
    A(cr,si*2)  = -1; % Pin
    A(cr,ei*2-1)= Res(i); % Q    
end

% flow conservation at bifurcations
bid = find(Gdegree==2);
for i = 1:size(bid,1)
    [eid,nid] = outedges(G,bid(i)); % gives 2 edges and 2 out endnode numbers
    cr = cr + 1;
    A(cr,2*bid(i)-1) = 1;
    A(cr,2.*nid(:)-1)= -1; % nid is 2 nodes
end

% uniform outflow BC
oid = find(Gdegree==0);
nout = size(oid,1);
for j = 1:nout
    cr = cr + 1;
    A(cr,2*oid(j)-1) = 1;
    %B(cr) = Qinp/nout;
    B(cr) = Qinp*10^6/nout; % m^3/s ->cm^3/s=m^3/10^6/\end
end

% remove fisrt two columns as they are for 1st node
 A(:,1:2)=[];
 
 % check
 %det(A)
 if det(A)== 0.0
     disp('  Warning: matrix is singular');
     return
 end
% cond(A)
 if cond(A)>10.e+16
     disp(['  Warning: problem is ill-conditioned: ']); %, num2str(cond(A))]);
     return
 end

%% Solve the system
% Asp = sparse(A); %https://www.mathworks.com/help/matlab/math/sparse-matrix-operations.html#f6-9169

X = A\B; %nodes from 2nd

P = zeros(nnode-1,1);
Q = zeros(nnode-1,1);

% solution is given at nodes, without segment info
Q(1) = Qinp;
P(1) = Pinp;
for i=1:nnode-1
    %Q(i+1) = X(i*2-1);
    %P(i+1) = X(2*i);
    Q(i+1) = X(i*2-1)/10^6; % test back cm^3/s->m^3/s
    P(i+1) = X(2*i)*133.32; %test back mmHg -> Pa
end

% Write p,q a endnode for each segment
% update Table
for j = 1:nsegm
    Table(j,2:3) = [Q(enode(j)),P(enode(j))];
end

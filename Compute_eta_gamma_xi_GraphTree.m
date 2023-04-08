%====================== SUBPROGRAMM: GET ETA, GAMMA, XI====================
% for each vessel in a tree
% ----
% Code by Hamid, 2018
%       - update by Vasilina, October 30, 2020
%==========================================================================

%% get arrays
eta = zeros(numOfBif,1);
gamma = zeros(numOfBif,1);
xi = zeros(numOfBif,1);

parent_segm = zeros(numOfBif,1);
xi_segm = zeros(nsegm,1);

c = 1;
for i = 1:numOfBif
    parentEnode = bifNodeID(i);
    
    % get 2 edges and 2 out endnode numbers
    [daughterEdges,daughterEnodes] = outedges(Graph,parentEnode); 
        
    % get Radii
    %parent_segm(i) = enode(:) == parentEnode; 
    %Rp = Radius(parent_segm(i));
    Rp = Radius(enode(:) == parentEnode);
    Rd1 = Radius(enode(:) == daughterEnodes(1));
    Rd2 = Radius(enode(:) == daughterEnodes(2));
    
    %get eta and gamma
    eta(c) = (Rd1^2 + Rd2^2) /  Rp^2;
    gamma(c) = min(Rd1,Rd2)/max(Rd1,Rd2);
    
    % get exponent xi from the Murray's law
    % but can use Matlab inbuilt function to colve nonliner equation...
    a = Rd1/Rp;
    b = Rd2/Rp;
    err = 1000;
    tol=1e-6;
    x = 2.5;
    while err > tol
        f = 1 - a^x - b^x;
        df = -a^x*log(a) -b^x*log(b);
        sol = f/df;
        err = abs(sol);
        x = x - sol;
    end
    xi(c) = x;
    xi_segm(enode(:) == parentEnode) = x; 
    c = c + 1;
       
end

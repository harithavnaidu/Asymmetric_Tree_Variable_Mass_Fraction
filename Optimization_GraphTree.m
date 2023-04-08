%==========================================================================
% This function solves the optimization problem of the arterial tree
function [Mt,Me,SMCtoCOL,R,L,TableOrder,p_mid,p_inp] = ...
          Optimization_GraphTree(q_parent,p_BC,segm,G,...
          Rorder,Lorder,RsdOrder,LsdOrder,RorderLow,RorderUp) 
% Calls functions:
%   - SolveHemidynamicsInGraphTree_BCPinQout
%   - mass_optimiz_tree_2
% PTflag=0: p_BC = pInp_parent at inlet, q at outlets
% PTflag=1: p_BC = p_terminal at outlets, q_parent at inlet
%
% WARNING: 
%   - assume that length is updated to keep the same L/R ratio! 
%       L(j) = R(j)*Lorder(order(j))/Rorder(order(j));
% NEED: 
%   -to check that initial L values are reasonable 
%   (possible do iterations for tree initializations to identify suitable L)
%   -to draw the resulted tree!
% -----------------
% Code by Hamid, 2018 
% Modified by Vasilina, December 29, 2020 
%==========================================================================

global mu PTflag

%% convert [m,s,kg] & Pa into [cm,s,gr] & mmHg
% q_parent = q_parent*10^6; % m3/s into cm3/s
% p_BC = p_BC/133.32; % Pa into mmHg
% Rorder = Rorder.*100; %cm
% Lorder = Lorder.*100; %cm

%% Initialize tree from Strahler ordering

% Strahler ordering data
order = segm(:,3); %order
maxord = max(order); minord = min(order);
norder = maxord - minord + 1;
order(:) = order(:) - minord + 1;
%norder = order(1); % first segment is of the last order
NumberOfVessels = size(segm,1); % number of segments

% ResOrder = zeros(norder,1);
% for i=1:norder
%     ResOrder(i) = 8*mu*Lorder(i)/(pi*Rorder(i)^4);
% end
%get ratio of Length to Radius

% % convert from Pa *s/cm^3 in mmHg s/cm3
% ResOrder = ResOrder./133.22;
%---------------
%% Initialization

% preallocate arrays
Mt = zeros(NumberOfVessels,1);
Me = zeros(NumberOfVessels,1);
Mtotal = zeros(NumberOfVessels,1);
SMCtoCOL = zeros(NumberOfVessels,1); 
phi_e = zeros(NumberOfVessels,1); 

R = zeros(NumberOfVessels,1);
L = zeros(NumberOfVessels,1);

Res = zeros(NumberOfVessels,1);   
TableOrder = zeros(NumberOfVessels,4);

% scale pressure and flow
% pT_parent = pInp_parent; % for initialization
p_mid = p_BC*ones(NumberOfVessels,1); 
p_term = p_BC*ones(NumberOfVessels,1); 
p_inp = p_BC*ones(NumberOfVessels,1); 
q = q_parent*ones(NumberOfVessels,1); 

% data = load("StrahlerGraphTreePinQout.mat"); %for initialization
% R = data.Radius;
% L = data.Length;

% intitalize R,L,Res and Table
for j=1:NumberOfVessels
      %version 1:
    L(j) = Lorder(order(j));
    R(j) = Rorder(order(j));
%     %Res(j) = ResOrder(order(j));
        
      %version 2:
    %use normal distribution with mean and standard deviations
    %rng('default') % randon number generator, for reproducibility
%     L(j) = normrnd(Lorder(order(j)),LsdOrder(order(j)));
%     R(j) = normrnd(Rorder(order(j)),RsdOrder(order(j)));
%     %take positive values only
%     L(j) = sign(L(j))*L(j);
%     R(j) = sign(R(j))*R(j);
%     if (R(j)== 0 || L(j) == 0)
%         L(j) = normrnd(Lorder(order(j)),LsdOrder(order(j)));
%         R(j) = normrnd(Rorder(order(j)),RsdOrder(order(j)));
%         L(j) = sign(L(j))*L(j);
%         R(j) = sign(R(j))*R(j);
%     end

    %verion 3:
%     %use normal distribution with mean and standard deviations
%     L(j) = Lorder(order(j));
%     %rng('default') % randon number generator, for reproducibility
%     R(j) = normrnd(Rorder(order(j)),RsdOrder(order(j)));
%     %take positive values only
%     R(j) = sign(R(j))*R(j);
%     if (R(j)== 0)
%         R(j) = normrnd(Rorder(order(j)),RsdOrder(order(j)));
%         R(j) = sign(R(j))*R(j);
%     end
    
    Res(j) = 8*mu*L(j)/(pi*R(j)^4);
    
    %[order, flow, terminal pressure, resistance];
    TableOrder(j,:) = [order(j),q(j),p_term(j),Res(j)];
end

% LtoR = L./R;

%remember initial values
Rinitial = R; Linitial = LengthSegmentk(R(:)); 

% Initialization of hemodynamics, update Table(:,2:3)
if PTflag==0 %for p_BC = pInp_parent
    TableOrder = SolveHemidynamicsInGraphTree_BCPinQout...
                 (G,segm,q_parent,p_BC,TableOrder);
elseif PTflag==1 %for p_BC = p_terminal
    TableOrder = SolveHemidynamicsInGraphTree_BCQinPout...
                 (G,segm,q_parent,p_BC,TableOrder);
%     disp("stop here, a function is under construction");
end

%update initalization
for j=1:NumberOfVessels
    q(j) = TableOrder(j,2);
    p_term(j) = TableOrder(j,3);
    p_inp(j) = p_term(j)+Res(j)*q(j);
    p_mid(j) = (p_term(j)+p_inp(j))/2 ;
end

Y0 = p_term;

%% Solve for terminal pressure to fit mass and geometry relations
% for lhat,R,L,Table

err = 1000;
c = 0; 
tol= 1e-6;

while err > tol %test
        
    for j=1:NumberOfVessels

        %update
        q(j) = TableOrder(j,2);
        p_term(j) = TableOrder(j,3);
        p_inp(j) = p_term(j)+Res(j)*q(j);
        p_mid(j) = (p_term(j)+p_inp(j))/2 ; 
        
        % optimization of the mass at each segment
        [Mtotal(j), R(j), phi_e(j), SMCtoCOL(j),cNR,tolNR] = ...
                          mass_optimiz_tree_2(p_mid(j),q(j),j,R(j));   
        Mt(j) = (1 - phi_e(j))*Mtotal(j);
        Me(j) = phi_e(j)*Mtotal(j);
        
        %update Strahler ordering based on numbers
        

        % keep L/R the same 
        % Warning: need to check if L/R is reasonable!
%         L(j) = R(j)*Lorder(order(j))/Rorder(order(j));
%         L(j) = LtoR(j) * R(j);

        %test fixed L
%         L(j) = Lorder(order(j));
%         L(j) = Linitial(j);
        
%         %test2: update order based on R and then update L based on new
%         %order
%         orderJ = order(j);
%         for k = norder:-1:1 
%             if (R(j) < RorderUp(k) && R(j) > RorderLow(k))
%                 %if (orderJ ~= k) %may give infinite loop of swithching
%                 %order back and forth
%                 if (orderJ > k)
%                     disp(['      update order from ', num2str(order(j)),...
%                          'to ',num2str(k)]);
%                     order(j) = k;
%                     break;                   
%                 end
%             end
%         end
%         
%         L(j) = Lorder(order(j));  %fixed L, updated order;
         %%too stip steps, too large Qs-Qd
         
%         %test3: keep L/R ratio fixed for new order
%         factor = 2.5; %to make the slope closer to LengthSegment
%         L(j) = factor*R(j)*Lorder(order(j))/Rorder(order(j)); 
        %good near liner L(R)
                
%         L(j) = normrnd(Lorder(order(j)),LsdOrder(order(j))); % diverges!!
        
        %test4: 
        L(j) = LengthSegmentk(R(j));
        
        % update resistance
        Res(j) = 8*mu*(L(j))/(pi* R(j)^4); 
        TableOrder(j,4) = Res(j);

    end

    % update steady hemodynamics for entire tree
    % with new resistance
    if PTflag==0    %for p_BC = pInp_parent
        TableOrder = SolveHemidynamicsInGraphTree_BCPinQout...
                       (G,segm,q_parent,p_BC,TableOrder);
    elseif PTflag==1 % for symm p_BC = p_terminal
%         disp("stop here, a function under construction");
        TableOrder = SolveHemidynamicsInGraphTree_BCQinPout...
                       (G,segm,q_parent,p_BC,TableOrder);
    end

    % update terminal pressure
    Y = TableOrder(:,3);
    
    % compute residual error
    err = norm(abs(Y-Y0));

    % update reference solution
    Y0 = Y;
    c = c+1;
    
    disp(['  end of iteration # ',num2str(c),', error ',num2str(err),]);

end %test

% need to draw the resulted tree!

disp('===End Optimization===');
% disp([' Radius Optimization at last vessel: Newton-Raphson iterations # ',...
%     num2str(cNR),', tollerance ',num2str(tolNR),]);
disp([' Total iterations # ', num2str(c),', tollerance ',num2str(tol)]);

end


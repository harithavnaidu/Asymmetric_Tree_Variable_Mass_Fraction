%% ======================================================================== 
% Function to build subgraph for orders next to a root element
% Algorithm follows [Yang,.. Marsden-2019] and [Spilker, 2006]
% Assumption:
%   - C(n,n)<1
%   - If C(n,n) accumulates by reminders one n-order branch, then parent 
%   element doesn't have two children in the end but instead one n-order 
%   element at midpoint
%   - update reminder remEl afeter exciding 1
% InOuts:
%   - AllSegmList - inout (update accumulated list of all segments)
%   - NewSegmList - out (Segments for this subgraph)
%   - cn - out (new total number of nodes)
% Issues:
%   - round-off error in 'floor' in double precision is used, changed loccaly
%   for single
% Last revision:
%   - by Vasilina July 26, 2019
%
%% ========================================================================
function [AllSegmList,NewSegmList,cn] = GetChildSubgraph...
                  (nOrder,nNodes,AllSegmList,connMatrix,...
                   numChildPerElInOrder,numElPerOrder)
%% ========================================================================
    %% Initialization
    
    % get local number of segments (ns) and vertices (nv)
    ns = numChildPerElInOrder(nOrder);
    nv = numChildPerElInOrder(nOrder);
    if numChildPerElInOrder(nOrder)== 0
        ns = ns + 1;
        nv = nv + 1;
    elseif numChildPerElInOrder(nOrder)== 1
        ns = ns + 2;
        nv = nv + 1;
    elseif numChildPerElInOrder(nOrder)>1
        ns = ns + numChildPerElInOrder(nOrder)-1;
    end
    % add a buffer for non-segmented parent element 
    % which later will be zero-out
    ns = ns + 1;
    
    % add a buffer for extra child element from accumalated children
    nv = nv + 1;
    
    % vertices of root elements of nOrder
    StartRoot = zeros(numElPerOrder(nOrder),1);
    EndRoot = zeros(numElPerOrder(nOrder),1);
    
    % reorder

    % find existent parent segments of nOrder
    % they will be first (np) element roots for nOrder
    % order is important: 
    ind = AllSegmList(:,3)== nOrder;
    np = sum(ind,1); %number of parent elements; =<numElPerOrder(nOrder)
    StartRoot(1:np) = AllSegmList(ind,1);
    EndRoot(1:np) = AllSegmList(ind,2);   
     
    % remove duplicates from original matrix
     AllSegmList(ind,:) = [];
      
    % initialize a global list of segments for a subgraph
    % SegmList[segm,1:3,elem]= [StartIndex EndIndex Order]
    SegmList = zeros(ns,3,numElPerOrder(nOrder)); 
    NewSegmList = zeros(ns*numElPerOrder(nOrder),3);
    %asign 1st segments to 1:np first nOrder elements
    %SegmList(1:np,1,:) = [StartRoot(1:np) EndRoot(1:np) nOrder];
    
%     SegmList(1,1,1:np) = StartRoot(1:np);
%     SegmList(1,2,1:np)= EndRoot(1:np);
%     SegmList(1,3,1:np)= nOrder;
              
    cn = nNodes; % nodes Global counter
    
    ChildPerElem = zeros(nOrder);
    % reminder for all orders of children segments
%     remEl(1:nOrder)= connMatrix(1:nOrder,nOrder)- ...
%                                        floor(connMatrix(1:nOrder,nOrder));
%     tempRem = remEl; % local reminder in each order;   
    tempRem(1:nOrder) = 0.0;      
    ElemStart = StartRoot(1);
    
    cp = 0; % count for existent parent elements
    
    %% ------------------  Main Loop -------------------------------------  
    for p = 1:numElPerOrder(nOrder)
        % recover old start nodes
        oldElemStart(p) = ElemStart;
        
        %---------- root element------------
        % local initialization of start-end vertices for element 
        StartIndex = zeros(nv,nOrder);
        EndIndex   = zeros(nv,nOrder);
        IntmIndex  = zeros(nv,1);
        
        %---------additional children--------
        % initialize a number of children
%         ChildPerElem(1:nOrder,p) = floor(connMatrix(1:nOrder,nOrder));  
        for i=1:nOrder
              if connMatrix(i,nOrder) > 0.0
                  a = connMatrix(i,nOrder) +  tempRem(i);
                  % to avoid round off error for floor use single precision
                  ChildPerElem(i) = floor(single(a));
                  tempRem(i) = a - ChildPerElem(i);
              end
        end
        
        % update a total number of children
        totalChildPerElem = sum(ChildPerElem(:),1);
        
%         if totalChildPerElem > 0
        
            % update k accoring to ChildPerElem
            k = find(ChildPerElem(:),1,'last');
            
            if isempty(k)
                k=0;
            end

            %------get start and end point of the element              
            ce = 0;  % edges local counter
            cm = 0;  % counter for midpoints

%             if k<nOrder
             if k<nOrder
                cp = cp + 1;
                if cp > np
                    disp('Warning:number of connected elements of this order is less than prescribed');
                    % possibly due to ommiting degrees higher than 4 (no trifurcations etc.)
                    break
                end
                ce = ce + 1; %existent vessel
                SegmList(1,1,p) = StartRoot(cp);
                SegmList(1,2,p) = EndRoot(cp);
                SegmList(1,3,p) = nOrder;
                % copy nodes
                StartIndex(1,nOrder) = SegmList(1,1,p);
                EndIndex(1,nOrder) =  SegmList(1,2,p);

            else %k==nOrder % add new children segment and its nodes as a
%             root; later attach start node to the parent element midpoint
%             elseif totalChildPerElem>0
                cn = cn + 1;
                StartIndex(1,nOrder) = cn;
                cn = cn + 1;
                EndIndex(1,nOrder)   = cn;
                ce = ce + 1;
                SegmList(ce,:,p)     = ...
                    [StartIndex(1,nOrder),EndIndex(1,nOrder),nOrder];
            end

            %update element start and end nodes 
            ElemStart = StartIndex(1,nOrder);
            ElemEnd = EndIndex(1,nOrder);

            % floating start of parent segm. 
            FloatMidStart = ElemStart; 

            %-------- adding segments -----------
            switch totalChildPerElem
                case 0 % for end segments of order 1
                    % only root segment exists which was added before
                    % so do nothing
                    disp('zero children');

                case 1 % add one midpoint of order k
                    
                    if k<nOrder
                        cm = cm + 1;
                        cn = cn + 1;
                        IntmIndex(cm) = cn;
                        cn = cn + 1;
                        EndIndex(1,k)   = cn;
                        ce = ce + 1;
                        SegmList(ce,:,p)  = [IntmIndex(cm),EndIndex(1,k),k];

                        % divide a parent element
                        ce = ce + 1;
                        SegmList(ce,:,p)  = [ElemStart,IntmIndex(cm),nOrder];
                        ce = ce + 1;
                        SegmList(ce,:,p)  = [IntmIndex(cm),ElemEnd,nOrder];
                        % zerout parent element index
                        SegmList(1,:,p) = 0; 
                        
                    else %k==nOrder %should define k==1?
%                     elseif k==1 % k==nOrder; attach 1st order children to this midpoint
                        % -- add a midpoint
                        % find previous element with different start node
                        pp = p; 
                        while ElemStart==oldElemStart(pp)
                            pp=pp-1;
                        end
                        oldp = pp-1;
                        if oldp==0 %should be for p>1
                            disp('Error: zero index');
                            return
                        end
                        ElemStart2 = ElemStart; % current start
                        %find old start node in previous element segment
                        st1 = SegmList(:,1,oldp)==oldElemStart(pp);  
                        %find  ElemStart2 in current segment
                        st2 = SegmList(:,1,p)==ElemStart2;
                        %swap start nodes in current and previous element
                        SegmList(st1,1,oldp)=ElemStart2;
                        SegmList(st2,1,p)=oldElemStart(pp);

                        %divide current element accordingly
                        cm = cm + 1;
                        IntmIndex(cm) = ElemStart2;
                        ce = ce + 1;
                        SegmList(ce,:,p) = [oldElemStart(pp),IntmIndex(cm),nOrder];
                        %zerout initial parent segment
                        SegmList(1,:,p) = 0;

                        % update element start and floating
                        ElemStart = oldElemStart(pp); % asign old
                        FloatMidStart = ElemStart2; 
                        
                        % connect last midpoint to the parent end
                        ce = ce + 1;
                        SegmList(ce,:,p)  = [FloatMidStart,ElemEnd,nOrder];

                    end

                otherwise 
                    % ----- add two children elements to the end
                    if (k < nOrder)
                        % at least ChildPerElem(k)>1
                        % 1st child:
                        cn = cn + 1;
                        StartIndex(1,k) = ElemEnd; 
                        EndIndex(1,k)   = cn;
                        ce = ce + 1;
                        SegmList(ce,:,p)  = [StartIndex(1,k),EndIndex(1,k),k];
                        % 2nd child:
                        cn = cn + 1;
                        ce = ce + 1;
                        if ChildPerElem(k)>=2
                            StartIndex(2,k) = ElemEnd; 
                            EndIndex(2,k)   = cn;
                            SegmList(ce,:,p)  = [StartIndex(2,k),EndIndex(2,k),k];
                        else % take next order (k-1), k>1
                            StartIndex(2,k-1) = ElemEnd; 
                            EndIndex(2,k-1)   = cn;
                            SegmList(ce,:,p)  = [StartIndex(2,k-1),EndIndex(2,k-1),k-1];
                            
                            %-- for next midpoints update k as this child was
                            % considered just above
                            k = k-1;
                        end
                    else %k==nOrder 
    %                 elseif p>1 %k==nOrder 

                        % -- add a midpoint
                        % find previous element with different start node
                        pp = p; 
                        while ElemStart==oldElemStart(pp)
                            pp=pp-1;
                        end
                        oldp = pp-1;
                        if oldp==0 %should be for p>1
                            disp('Error: zero index');
                            return
                        end
                        ElemStart2 = ElemStart; % current start
                        %find old start node in previous element segment
                        st1 = SegmList(:,1,oldp)==oldElemStart(pp);  
                        %find  ElemStart2 in current segment
                        st2 = SegmList(:,1,p)==ElemStart2;
                        %swap start nodes in current and previous element
                        SegmList(st1,1,oldp)=ElemStart2;
                        SegmList(st2,1,p)=oldElemStart(pp);

                        %divide current element accordingly
                        cm = cm + 1;
                        IntmIndex(cm) = ElemStart2;
                        ce = ce + 1;
                        SegmList(ce,:,p) = [oldElemStart(pp),IntmIndex(cm),nOrder];
                        %zerout initial parent segment
                        SegmList(1,:,p) = 0;

                        % update element start and floating
                        ElemStart = oldElemStart(pp); % asign old
                        FloatMidStart = ElemStart2; 

                        % -- add bifurcation of k-1 and lower orders
                        if k>1
                            % 1st child
                            cn = cn + 1; 
                            StartIndex(1,k-1) = ElemEnd;
                            EndIndex(1,k-1)   = cn;
                            ce = ce + 1;
                            SegmList(ce,:,p)  = [StartIndex(1,k-1),EndIndex(1,k-1),k-1];
                            % 2nd child
                            cn = cn + 1; 
                            ce = ce + 1;
                            if ChildPerElem(k-1)>=2 % includes k==2 as C(1,2)>=2
                                StartIndex(2,k-1) = ElemEnd;
                                EndIndex(2,k-1)   = cn;
                                SegmList(ce,:,p)  = [StartIndex(2,k-1),EndIndex(2,k-1),k-1];
                            elseif k>2 % take next order (k-2), for k>2
                                StartIndex(2,k-2) = ElemEnd; 
                                EndIndex(2,k-2)   = cn;
                                SegmList(ce,:,p)  = [StartIndex(2,k-2),EndIndex(2,k-2),k-2];
                            else
                                disp('warning: out of bounds');
                                break
                            end
                        end

                        %-- for next midpoints update k as this child was
                        % considered just above
                        k = k-1;
                    end            

                    %------- add other midpoints     
                    if ChildPerElem(k) > 2 
                        for j = 3:ChildPerElem(k)
                            % midpoint at parent element
                            cn = cn + 1;                
                            cm = cm + 1;
                            IntmIndex(cm)  = cn;
                            StartIndex(j,k) = cn;
                            cn = cn + 1;
                            EndIndex(j,k)   = cn;

                            % divide a parent element 
                            ce = ce + 1;
                            SegmList(ce,:,p)  = [FloatMidStart,IntmIndex(cm),nOrder];
                            % update floating start
                            FloatMidStart = IntmIndex(cm);

                            % add a child element
                            ce = ce + 1;
                            SegmList(ce,:,p)  = ...
                                   [StartIndex(j,k),EndIndex(j,k),k];
                        end
                    end

                    % add (k-1) order midpoints and children
                    if k>1 
                        for i = k-1:-1:1 
                            if (ChildPerElem(i) > 0) 
                                for j=1:ChildPerElem(i)
                                    % add midpoint
                                    cn = cn + 1;
                                    cm = cm + 1;
                                    IntmIndex(cm) = cn;
                                    % add child
                                    StartIndex(j,i) = cn;
                                    cn = cn + 1;
                                    EndIndex(j,i)   = cn;
                                    % add child segment
                                    ce = ce + 1;
                                    SegmList(ce,:,p)  = [StartIndex(j,i),EndIndex(j,i),i];

                                    % add segment to parent element
                                    ce = ce + 1;
                                    SegmList(ce,:,p)  = [FloatMidStart,IntmIndex(cm),nOrder];
                                    % update next midpoint start
                                    FloatMidStart = IntmIndex(cm);                  
                                end
                            end
                        end
                    end

                    % connect last midpoint to the parent end
                    ce = ce + 1;
                    SegmList(ce,:,p)  = [FloatMidStart,ElemEnd,nOrder];

                    % zero-out original parent element
                    SegmList(1,:,p)  = 0;

            end
            % to avoid duplicates remove this segment from global list
            for i=1:p
                ind1 = SegmList(:,1,i)==ElemStart;
                ind2 = SegmList(:,2,i)==ElemEnd;
                ind3 = SegmList(:,3,i)==nOrder;
                ind = ind1 & ind2 & ind3;
                if ind
                    SegmList(ind,:,i)=[];
                end
            end
%         end  
    end
    
    %save
    nn=0;
    for pp=1:numElPerOrder(nOrder)
        cep = size(SegmList(:,:,pp),1);
        
        % shuffle position of segments in this element 
        % to introduce ramdomness of children branches positions
        temp=randperm(size(SegmList,1));
        SegmList(:,:,pp) = SegmList(temp,:,pp);
        
        NewSegmList(nn+1:nn+cep,:)= SegmList(:,:,pp);
        nn = nn+cep;
    end
    
end

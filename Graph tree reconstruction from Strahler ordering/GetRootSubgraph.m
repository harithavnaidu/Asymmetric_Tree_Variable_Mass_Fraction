%% ======================================================================== 
% Create nodes and edges for Root subgraph of the leadgin order
% Assumptions (based on convenience and Huang's data)
%   - Maximum two elements of the leading order
%   - C(1,2)>=2 as in Huang's data
%   - numbver of segments per element is >=2
%   - C(n,n)<1
%   - C(m,n)>0  for firstNonzeroOrder < m < lastNonzeroOrder
%
%InOuts:
%   -NewSegmList - out (Segments for this subgraph)
%   -cn - out (new total number of nodes)
%
%Created:
%   - by Vasilina on June 2019
%Lastg edit 
%   - by Vasilina on July 26, 2019
%% ========================================================================
function [NewSegmList,cn] = GetRootSubgraph...
                (maxOrder,connMatrix,numChildPerElInOrder,numElPerOrder)  
%% ========================================================================
            %% Initialization
    % (maybe not accurate)get local number of segments (ns) and vertices (nv)
    % in the parent subgraph
    ns = numChildPerElInOrder(maxOrder); 
    nv = numChildPerElInOrder(maxOrder);
    if numChildPerElInOrder(maxOrder)== 0
        ns = ns + 1; 
        nv = nv + 1;
    elseif numChildPerElInOrder(maxOrder)== 1
        ns = ns + 2; 
        nv = nv + 1;
    elseif numChildPerElInOrder(maxOrder)>1
        ns = ns + numChildPerElInOrder(maxOrder)-1; 
    end
    % add a buffer for non-segmeneted parent element 
    % which later will be zero-out
    ns = ns + 1;
    % add a buffer for extra child element from accumalated children
    nv = nv + 1;
    
    % initialize a global list of segments for a subgraph
    % SegmList[a,1:3,b]= [StartIndex EndIndex Order]
    SegmList   = zeros(ns,3,numElPerOrder(maxOrder)); 
    NewSegmList = zeros(ns*numElPerOrder(maxOrder),3);
      
    % reminder for all orders of children segments
    remEl(1:maxOrder) = connMatrix(1:maxOrder,maxOrder)- ...
                                   floor(connMatrix(1:maxOrder,maxOrder));
    tempRem = remEl; % local reminder in each order;   
    
    nSegm = 0;
    
    %set 1st root element start and end; will be updated 
    ElemStart1 = 1;
    ElemEnd1 = 2;  %1st element end

    cn = 2; % nodes Global counter  
    
    %% Main Loop
    for p = 1:numElPerOrder(maxOrder)
    
        % Local initialization of start-end vertices for an element 
        StartIndex = zeros(nv,maxOrder); 
        EndIndex   = zeros(nv,maxOrder);
        IntmIndex  = zeros(nv,maxOrder); 
              
        if p==1
            ElemStart = ElemStart1;
            ElemEnd = ElemEnd1;
            StartIndex(1,maxOrder)=ElemStart;
            EndIndex(1,maxOrder)=ElemEnd;
            ce = 1;
            SegmList(1,:,p) = [ElemStart,ElemEnd,maxOrder];
 
        elseif p==2 % no more?
            cn = cn + 1;
            ElemStart = cn;
            cn = cn + 1;
            ElemEnd = cn;
            ce = 1;   
            SegmList(1,:,p) = [ElemStart,ElemEnd,maxOrder];
        end
                 
        FloatMidStart = ElemStart; %floating start of parent segm.  

        cm = 0; % local counter for midpoints 
        %---------additional children--------
        % initialize a number of children
        ChildPerElem(1:maxOrder) = floor(connMatrix(1:maxOrder,maxOrder)); 
        % check if need to add one more child segment to the current elem.       
        % update tempRem for current element
               
        for i=1:maxOrder
              if tempRem(i)>0
                 if tempRem(i)>= 1
                     ChildPerElem(i) = ChildPerElem(i) + 1;
                     tempRem(i) = tempRem(i)-1;
                 end
                 tempRem(i) = remEl(i) + tempRem(i);     
              end
        end
        
       % update a total number of children
       totalChildPerElem = sum(ChildPerElem); 
        
        % update k accoring to ChildPerElem
        k = find(ChildPerElem,1,'last');
             
        % get nodes and segments for the parent element
        switch totalChildPerElem
                            
            case 1 % one mid-point
                cm = cm + 1;
                cn = cn + 1;
                IntmIndex(cm,maxOrder) = cn;
                StartIndex(1,k) = cn;
                cn = cn + 1;
                EndIndex(1,k)   = cn;

                % add child edge
                ce = ce + 1;
                SegmList(ce,:,p)  = [StartIndex(1,k),EndIndex(1,k),k];

                % add segments to parent edges
                ce = ce + 1;
                SegmList(ce,:,p)  = [ElemStart,IntmIndex(1,maxOrder),maxOrder];
                ce = ce + 1;
                SegmList(ce,:,p)  = [IntmIndex(1,maxOrder),ElemEnd,maxOrder];
                
                % zerout initial parent segment
                SegmList(1,:,p) = 0;

            otherwise  
                if k<maxOrder % add two children elements and midpoints
                    % 1st child
                    cn = cn + 1; 
                    StartIndex(1,k) = ElemEnd;
                    EndIndex(1,k)   = cn;
                    ce = ce + 1;
                    SegmList(ce,:,p)  = [StartIndex(1,k),EndIndex(1,k),k];
                    % 2nd child
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
                    end
                else %k=maxOrder, for p==2
                     
%                     % maximum one child then add kth order element at the
%                     % midpoint of this element
%                     % also add two daughters of the lower order to the end
%                     
                    % -- add a midpoint
                    
                    % for convenience: exchange starting nodes
                    % of this elementwith the previous
                    ElemStart2=ElemStart;
                    st1 = NewSegmList(:,1)==ElemStart1;
                    st2 = SegmList(:,1,p)==ElemStart2;
                    %swap
                    NewSegmList(st1,1)=ElemStart2;
                    SegmList(st2,1,p)=ElemStart1;
                    ElemStart = ElemStart1;
                    ElemStart1 = ElemStart2;
                                        
                    % add this parent as a midpoint to previous
                    % parent element
                    cm = cm + 1;
                    IntmIndex(cm,maxOrder) = ElemStart1;
                    %divide a segment
                    ce = ce + 1;
                    SegmList(ce,:,p)  = [ElemStart,IntmIndex(cm,maxOrder),maxOrder];
                    ce = ce + 1;
                    SegmList(ce,:,p)  = [IntmIndex(cm,maxOrder),ElemEnd,maxOrder];
                    % zerout initial parent segment
                    SegmList(1,:,p) = 0;
                                                            
                    %update element start and floating point                
                    ElemStart = IntmIndex(cm,maxOrder);
                    FloatMidStart = ElemStart;
                                                          
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
                            return
                        end
                    end
                    
                    %-- for next midpoints update k as this child was 
                    % considered just above
                    k = k-1;                  
                end   
                
                % --- add other midpoints 
                if ChildPerElem(k) > 2
                    for j = 3:ChildPerElem(k)
                        % midpoint at parent element
                        cn = cn + 1;                
                        cm = cm + 1;
                        IntmIndex(cm,maxOrder)  = cn;
                        StartIndex(j,k) = cn;
                        cn = cn + 1;
                        EndIndex(j,k)   = cn;

                        % divide a parent element 
                        ce = ce + 1;
                        SegmList(ce,:,p)  = [FloatMidStart,IntmIndex(cm,maxOrder),maxOrder];
                        FloatMidStart = IntmIndex(cm,maxOrder);

                        % child element
                        ce = ce + 1;
                        SegmList(ce,:,p)  = [StartIndex(j,k),EndIndex(j,k),k];
                    end
                end

                % add (k-1) order midpoints and childs
                if k>1
                    for i = k-1:-1:1
                        if (ChildPerElem(i) > 0) 
                            for j=1:ChildPerElem(i)
                                % add midpoint
                                cn = cn + 1;
                                cm = cm + 1;
                                IntmIndex(cm,maxOrder) = cn;
                                % add child
                                StartIndex(j,i) = cn;
                                cn = cn + 1;
                                EndIndex(j,i)   = cn;
                                % add child segment
                                ce = ce + 1;
                                SegmList(ce,:,p)  = [StartIndex(j,i),EndIndex(j,i),i];

                                % divide a parent element
                                ce = ce + 1;
                                SegmList(ce,:,p)  = [FloatMidStart,IntmIndex(cm,maxOrder),maxOrder];
                                FloatMidStart = IntmIndex(cm,maxOrder);                  
                            end
                        end
                    end
                end
             

                % connect last midpoint to the parent end
                ce = ce + 1;
                SegmList(ce,:,p)  = [FloatMidStart,ElemEnd,maxOrder];

                % zero-out original parent element
                SegmList(1,:,p)  = 0;

        end
        
        % shuffle position of segments in this element 
        % to introduce ramdomness of children branches positions
        temp=randperm(ce-1);
        SegmList(2:ce,:,p) = SegmList(temp+1,:,p);
        
        % segments of the subgraph NewSegmList(numElPerOrder(maxOrder)*nsl,3)
        % ce>=2
        NewSegmList(nSegm+(ce-1)*(p-1)+1:nSegm+(ce-1)*p,:)=SegmList(2:ce,:,p);
        
        % to avoid duplicates:
        % remove this segment from global list
        ind1 = NewSegmList(:,1)==ElemStart;
        ind2 = NewSegmList(:,2)==ElemEnd;
        ind3 = NewSegmList(:,3)==maxOrder;
        ind = ind1 & ind2 & ind3;
        NewSegmList(ind,:)=[];
        nSegm = size(NewSegmList,1);
        
    end
   
end
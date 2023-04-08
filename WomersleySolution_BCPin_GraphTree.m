%==========================================================================
% Get Deformable wall Womersley solution for longitudinally constrained 
% vesses in each branch of asymmetric tree
function [InpPressureTime,TermPressureTime,qInpTime,qTermTime,...
          WaveSpeed,Alpha,zinp,zterm,zchar,InpImpedance]...
          = WomersleySolution_BCPin_GraphTree...
            (T,PnInp,pInpSteady,Table,R,L,h,Ctt,...
             Nt,NumModes,GammaVal,segm,G)
% Advanced for stiffness component Ctttt instead of isotopic (E,sigma)
% where C = [Ctttt,Cttzz;Cttzz,Czzzz]./((Mt+Me)/(0.3*rho_w))
%
% Table(j,:)=[order, flow, terminal pressure, resistance];
%
% Inlet BC: flow waveform (PnInp)
%
% Calls functions:
%   - WomersleySolutionCoeffTethered
%   - CharacteristicImpedanceOneVessel
%   - InputImpedanceOneVessel
%   - TerminalImpedanceVessel
%   - FilterAndIFFT
%   - TerminalImpedanceJoint
%   - ReflectionCoeff
%   - HforwardCoeff
%   - Qrecovery
%   - Precovery
%-------------------
% Code by Vasilina, April 16, 2020
%       - updated on October 30, 2020
%==========================================================================
global mu rho
%% Get graph layers, descending
snode = segm(:,1); %start node of the segment
enode = segm(:,2); %end node of the segment
% order = segm(:,3); %order
% maxord = max(order); minord = min(order);
% norder = maxord - minord + 1;
% order(:) = order(:) - minord + 1;

Gdegree = outdegree(G);

plotGraph = plot(G,'Layout','layered');
% node <-> layer vector (#1 node <-> 6th layer)
layer = plotGraph.YData; % each node is assigned to the layer number
nlayer = max(layer); % number of layers
% layer numberings is bottom-up

%% --------- Preallocation and initialization-------------
NumberOfVessels = size(Table,1);
qSteady = Table(:,2); 
pSteadyTerm = Table(:,3);

%Impedance is zero at first mode
GammaList = zeros(NumModes,NumberOfVessels);
CharactImpedance= zeros(NumModes,NumberOfVessels);
InpImpedance = zeros(NumModes,NumberOfVessels);
TermImpedance = zeros(NumModes,NumberOfVessels);
WaveSpeed = zeros(NumModes,NumberOfVessels);
zinp = zeros(Nt,NumberOfVessels);
zterm = zeros(Nt,NumberOfVessels);
zchar = zeros(Nt,NumberOfVessels);
Alpha = zeros(NumModes,NumberOfVessels);
omegan = zeros(NumModes,1);

%% ---------------- Backward Computation of Impedance ---------------------
%Set Frequency Modes
for m=1:NumModes
	omegan(m) = 2.0*pi*(m-1)/T;
end
 
% --- Get Zinp at terminal vessels-----------------------------------------
for j = 1:nlayer %from last gen to first
    layernnodes = find(layer(:)==j); % end nodes at jth layer
    
    for i = 1:size(layernnodes,1)
        ei = layernnodes(i); %end node
        segid = find(enode(:)==ei);% segment id
   
        % go over the segments
        if indegree(G,ei)>0 %exclude 1st starting node

            %outlets
            if (Gdegree(ei)==0) % outlets
                GammaList(:,segid)=GammaVal;

                %get Womersley coefficients
                [cn,Mn,gn,c_Rn,alphan] = WomersleySolutionCoeffTethered(...
                     NumModes,omegan,mu,rho,R(segid),Ctt(segid),h(segid));
                WaveSpeed(:,segid) = cn;  % c_Rn is real part of PWV
                Alpha(:,segid)=alphan;

                %get characteristic impedance
                [Z0] = CharacteristicImpedanceOneVessel...
                       (NumModes,R(segid),rho,cn,Mn,gn); 
                CharactImpedance(2:NumModes,segid) = Z0(2:NumModes); %oscil.

                %get input impedance
                [Zinp] = InputImpedanceOneVessel...
                         (NumModes,omegan,L(segid),Z0,cn,GammaList(:,segid)); 
                InpImpedance(2:NumModes,segid) = Zinp(2:NumModes); % oscil.         

                %get terminal impedance, for figures
                [ZT] = TerminalImpedanceVessel(NumModes,Z0,GammaList(:,segid));
                TermImpedance(1,segid) = 0.0;      
                TermImpedance(2:NumModes,segid) = ZT(2:NumModes); 

                %get time domain functions (steady part is zero)
                [zchar(:,segid)] = FilterAndIFFT...
                                  (NumModes,Nt,CharactImpedance(:,segid));
                [zinp(:,segid)] = FilterAndIFFT...
                                  (NumModes,Nt,InpImpedance(:,segid));                                   
                [zterm(:,segid)] = FilterAndIFFT...
                                  (NumModes,Nt,TermImpedance(:,segid)); 

                % add steady part to input and term impedance 
                zinp(:,segid)= zinp(:,segid) + pInpSteady(segid)/qSteady(segid);
                zterm(:,segid)= zterm(:,segid) + pSteadyTerm(segid)/qSteady(segid);

            else % other nodes
                % find duaghters segments of segid
                dsegid = find(snode(:)==ei); % two daughters

                % as Zinp at j1 and j2 are computed at preveious generation
                % compute terminal impedance 
                [ZT] = TerminalImpedanceJoint(NumModes,...
                       InpImpedance(:,dsegid(1)),InpImpedance(:,dsegid(2)));
                TermImpedance(2:NumModes,segid) = ZT(2:NumModes);   

                %get Womersley coefficients
                [cn,Mn,gn,c_Rn,alphan] = WomersleySolutionCoeffTethered(...
                       NumModes,omegan,mu,rho,R(segid),Ctt(segid),h(segid));
                WaveSpeed(:,segid) = cn; 
                Alpha(:,segid)=alphan;

                %get Z0
                [Z0] = CharacteristicImpedanceOneVessel...
                       (NumModes,R(segid),rho,cn,Mn,gn);
                CharactImpedance(2:NumModes,segid) = Z0(2:NumModes);

                %get reflection coefficient
                [Gamma] = ReflectionCoeff(NumModes,ZT,Z0);
                GammaList(:,segid) = Gamma(:);

                %get input impedance
                [Zinp] = InputImpedanceOneVessel...
                         (NumModes,omegan,L(segid),Z0,cn,Gamma); 
                InpImpedance(2:NumModes,segid) = Zinp(2:NumModes);

                %get time domain functions (steady part is zero)
                [zchar(:,segid)] = FilterAndIFFT...
                                   (NumModes,Nt,CharactImpedance(:,segid));
                [zinp(:,segid)] = FilterAndIFFT(NumModes,Nt,InpImpedance(:,segid));                                 
                [zterm(:,segid)] = FilterAndIFFT(NumModes,Nt,TermImpedance(:,segid)); 

                % add steady part to input and term impedance 
                zinp(:,segid)= zinp(:,segid) + pInpSteady(segid)/qSteady(segid);
                zterm(:,segid) = zterm(:,segid) + pSteadyTerm(segid)/qSteady(segid);            
            end   
        end
    end
end

%% --------------- Forward Computation of Total Pressure and Flow -----------
Pinp = zeros(1,NumModes);
TermPressure = zeros(NumModes,NumberOfVessels);
TermPressureTime = zeros(Nt,NumberOfVessels);
InpPressureTime = zeros(Nt,NumberOfVessels);
qInpTime = zeros(Nt,NumberOfVessels);
qTermTime = zeros(Nt,NumberOfVessels);  

InpPressureTimeConv = zeros(Nt,NumberOfVessels);
TermPressureTimeConv =  zeros(Nt,NumberOfVessels);

%----------------Get solution at the root vessel---------------------------
% Input pressure in frequency domain is given PnInp

% test: comment out Qinp, use Pinp instead:
% % Get root input flow in frequency domain
% Qinp = zeros(NumModes,1);
% % steady
% Qinp(1) = qSteady(1); % q_parent
% % oscillatory
% for m=2:NumModes
%     if (InpImpedance(m,1)>1.e-4)
%         Qinp(m) = PnInp(m)/InpImpedance(m,1);
%     else
%         disp(['Warning: Devision by small number: InpImpedance(m,1)']);
% %         ,...
% %              InpImpedance(m,1)]);
%     end
% end
% %get root total input flow in time domain
% [qInpTime(:,1)] = FilterAndIFFT(NumModes,Nt,Qinp);
 
%get Hf coefficient, for oscillations only
[Hfn] = HforwardCoeff(NumModes,omegan,L(1),...
                      PnInp,WaveSpeed(:,1),GammaList(:,1));
                  
%get root inflow in time domain, x=0, segment id is 1
[Qn]=Qrecovery(NumModes,omegan,L(1),0,Hfn,WaveSpeed(:,1),...
                CharactImpedance(:,1),GammaList(:,1));
Qn(1) = qSteady(1); 
[qInpTime(:,1)] = FilterAndIFFT(NumModes,Nt,Qn);

%get root outflow in time domain, x=L, segment id is 1
[Qn]=Qrecovery(NumModes,omegan,L(1),L(1),Hfn,WaveSpeed(:,1),...
                CharactImpedance(:,1),GammaList(:,1));
Qn(1) = qSteady(1); 
[qTermTime(:,1)] = FilterAndIFFT(NumModes,Nt,Qn);
                  
%get root terminal P in frequency domain
[PnT]=Precovery(NumModes,omegan,L(1),L(1),Hfn,WaveSpeed(:,1),GammaList(:,1));
TermPressure(1,1) = pSteadyTerm(1);   %     %steady
TermPressure(2:NumModes,1)=PnT(2:NumModes);   %oscillatory
          
%get root terminal pressure in time domain
[TermPressureTime(:,1)] = FilterAndIFFT(NumModes,Nt,TermPressure(:,1)); 

% %get input pressure,x=0 (for test)
[InpPressureTime(:,1)] = FilterAndIFFT(NumModes,Nt,PnInp); 

%--------------------- Get total p and q for all branches------------------

% forward loop, layers in reverse order with generations
for j=nlayer-1:-1:1 % start from the second generation
    layernnodes = find(layer(:)==j); % end nodes at jth layer

    for i=1:size(layernnodes,1)
        si = layernnodes(i); % start node
        % get parent segment id
        parsegid = find(enode(:)==si); 
        
        % go over the segments
        if outdegree(G,si)>0
            % get 2 daughter-segments id
            segid = find(snode(:)==si);
            
            for ii=1:2 % daughter segments
                dsegid = segid(ii);
                       
                % get input P and Q, x=0, at dsegid
                Pinp(:) = TermPressure(:,parsegid);   
                % pressure conservation
                InpPressureTime(:,dsegid) = TermPressureTime(:,parsegid);

                % get Hf coefficient, for oscillations only
                [Hfn] = HforwardCoeff(NumModes,omegan,L(dsegid),...
                        Pinp,WaveSpeed(:,dsegid),GammaList(:,dsegid));            
                % get input flow
                [Qn] = Qrecovery(NumModes,omegan,L(dsegid),0.0,Hfn,...
                    WaveSpeed(:,dsegid),CharactImpedance(:,dsegid),...
                    GammaList(:,dsegid)); 
                Qn(1) = qSteady(dsegid); 
                [qInpTime(:,dsegid)] = FilterAndIFFT(NumModes,Nt,Qn);

                %get terminal P and Q, x=L(j)
                [PnT]=Precovery(NumModes,omegan,L(dsegid),L(dsegid),Hfn,...
                               WaveSpeed(:,dsegid),GammaList(:,dsegid));
                TermPressure(1,dsegid) = pSteadyTerm(dsegid);  
                TermPressure(2:NumModes,dsegid) = PnT(2:NumModes);
                [TermPressureTime(:,dsegid)]=FilterAndIFFT(NumModes,Nt,...
                                         TermPressure(:,dsegid)); 

                [Qn]=Qrecovery(NumModes,omegan,L(dsegid),L(dsegid),Hfn,...
                    WaveSpeed(:,dsegid),CharactImpedance(:,dsegid),...
                    GammaList(:,dsegid)); 
                Qn(1) = qSteady(dsegid); 
                [qTermTime(:,dsegid)] = FilterAndIFFT(NumModes,Nt,Qn);

                % alternative p: convolution
                TermPressureTimeConv(:,dsegid)=cconv(qTermTime(:,dsegid),zterm(:,dsegid),Nt)./Nt;
                InpPressureTimeConv(:,dsegid)= cconv(qInpTime(:,dsegid),zinp(:,dsegid),Nt)./Nt; 
            end            
        end
    end
end

% test conservation of flow at bifurcation, Okay!
% % small tree: parent id=1, daughters id=2,6
% figure(101) %test
% plot(qTermTime(:,1),'b'); hold on
% plot(qInpTime(:,2),'m.-'); hold on
% plot(qInpTime(:,6),'m.'); hold on
% plot(qInpTime(:,2)+qInpTime(:,6), 'r.'); hold off
% % plot(2*qInpTime(:,2), 'r.'); hold off

end %function
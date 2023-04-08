%%=========================================================================
%-------------- ASYMMETRIC TREE, VARIABLE MASS FRACTION -------------------
% Homeostatic Optimization and Hemodynamics in arterial tree
%   - Extension of Murray's Law for homeostasis (steady flow)
%   - Bifurcating tree defined by graph from Strahler ordering system
%   - Axisymmetric straight cylindrical vessels
%   - Pulsitile hemodynamics (postprocessing) based on longitudinally
%     constrained deformable wall Womersley's theory
%
% Assumption: vessel length is updated to keep the same L/R ratio! 
%
% Boundary conditions: (PinQout)
%   - Steady state optimization problem in a tree:
%       > Inlet of the root vessel: mean pressure
%       > Terminal vessels outlets: uniform mean flow
%   - Pulsatile hemodynamics in a tree:
%       > Inlet of the root vessel: pressure waveform
%       > Terminal vessels outlets: reflection coefficient
%
% Boundary conditions: (QinPout) diverges, may need ot have extra
% constraints and another implementation of the optimization
%   - Steady state optimization problem in a tree:
%       > Inlet of the root vessel: mean flow
%       > Terminal vessels outlets: uniform mean pressure
%   - Pulsatile hemodynamics in a tree:
%       > Inlet of the root vessel: flow waveform
%       > Terminal vessels outlets: reflection coefficient
%
% Flags to switch btw symm/asymm: pruning (prunflag) and mean pressure at
% terminal outlets (PTflag)
%   - PinQout: PTflag = 0;
%              Provide mean pInp_parent; q_term (or q_parent)
%   - QinPout: PTflag = 1;
%              Provide mean q_parent; p_terminal
%
% Note: 
%   - System of Units (m,s,kg) Pa
%   - mmHg = 1333.2 dyn/cm2 = 133.32 Pa
%
% Input file: FlowBZ.dat, Pinlet.dat
%             GraphSegmentsNodes.inp
%             StrahlerOrderDandL.inp, StrahlerOrderSDforDandL.inp
% Output file: Strahler_GraphTree.mat
%             HomeostaticRandL.dat
%
% Calls functions:
%   - Optimization_GraphTree
%   - ZeroP
%   - YoungMod_2
%   - FilterAndFFT
%   - WomersleySolution_BCPin_GraphTree
% Calls subprogramms:
%   - Compute_eta_gamma_xi_GraphTree
%   - PlottingResults_GraphTree
%
%-------------------------------------------------------------------------
% Code by Hamid Gharahi, 2018 (optimization)
% and Vasilina Filonova, 2018 (tree, hemodynamics)
%   - Last updates on January, 2021
%
%% ========================================================================
clear; close all; clc;

%% =========================== PARAMETERS =================================
%------------------------ Wall Tissue Parameters --------------------------
% Assumptions: 
%   1) we have 4 collagen fiber families in 0, 90, 45, and -45 degrees. 
%   2) SMCs that are circumferentially oriented. 
%   3) isotropic elastin matrix.
%   4) The mass fraction in generations is known.
%   5) Ratio of collagen to SMC is constant throughout the arterial tree. 
%   6) The deposition stretch of SMC and Collagen fibers are constant.
%   8) Stiffness components are computed using Small on Large theory (SoL).
%   9) The metabolic cost of maintenance of collagen and SMCs is the same 
%   and constant. 
%   7) Blood viscosity is constant. 

% Set some global parameters:
global H0 rho_w rho beta gamma alpha_t mu lmax lmin S ...
       Ghe1 Ghe2 Ghc Ghm fiber_ang P_h PTflag R0 
       %Rmin pruneflag eta xi

% Metabolic costs
alpha_t = 2*746;    % W/m^3, Taber 1998 (Biophys J.), porcine carotid artery
beta = 160;         % W/m^3, Lindstrom et al (2014)
gamma = 0.00891267681;    % J*s/m^3, Lindstrom et al (2014)

% Active wall
S =  2.0e+004;
lmax = 1.2;     %  nondim., maximumu tension in active tone
lmin = 0.7;     %  nondim., minimum tension in active tone

% Wall density
rho_w = 1060;   % kg/m3 

% Mechanical_properties (homest. prestretches, angle, pressure) 
% Ghe1, Ghe2, Ghc are updated in 'mechanical_properties_PA_4' for 
% each segment
Ghe1 = 1.27; 
Ghe2 = 1.27;
Ghc = 1.154; %adjust for h/R fit??
Ghm = 1.21;
fiber_ang = 45*pi/180; % angle of helical collagen fibers
% homeostatic reference pressure needed for calibration
% of prestretches in 'mechanical_properties_PA_4'
P_h = 12*133.32;        % 12 mmHg, [ref.] 

%---------------------- Tree Structure Parameters -------------------------
% Get Strahler ordering info for a vascular tree
% create directed graph out from it

% Load and plot graph data
segm = load("GraphSegmentsNodes.inp");
snode = segm(:,1); %start node of the segment
enode = segm(:,2); %end node of the segment
order = segm(:,3); %order

maxord = max(order);
minord = min(order);
norder = maxord - minord + 1;%order(1); % first segment is of the last order
order(:) = order(:) - minord + 1; %adjust ordering
nsegm = size(segm,1); % number of segments

% get directed graph
Graph = digraph(snode,enode);
% check
%plot(Graph);

plotGraph = plot(Graph,'Layout','layered');
% node <-> layer vector (#1 node <-> 6th layer)
layer = plotGraph.YData; % each node is assigned to the layer number
nlayer = max(layer); % number of layers

% get outdegree
% 1 for inlet; 0 for outlets; 2 for other nodes
Gdegree = outdegree(Graph); % for node numbers

% Load diameters and lenghts means vs Strahler orders (in cm)
vesselSizeData = load("StrahlerOrderDandL.inp");
%sort on ascending order
vesselSizeData = flip(vesselSizeData);
% Load diameters and lenghts vs Strahler orders
% in cm
R = 0.5.*vesselSizeData(:,2);
L = vesselSizeData(:,3);
%from cm to m
L = L./100; R = R./100;

%Load dimater and length starndard deviation (in cm)
vesselSizeSD = load("StrahlerOrderSDforDandL.inp");
vesselSizeSD = flip(vesselSizeSD);
Rsd = 0.5.*vesselSizeSD(:,2);
Lsd = vesselSizeSD(:,3);
%from cm to m
Lsd = Lsd./100; Rsd = Rsd./100;
%let's scale Lsd an order down otherwise it is too large
% Lsd = Lsd./10;
% Lsd = 0; Rsd = 0;

%get R borders for order definition ( [Huang 1996] eqs.(2),(3))
Rlow = zeros(norder,1);
Rup = zeros(norder,1);
for i=1:norder
    if (norder == 1)
        disp("Error: just one order, need more");
        exit;
    end
    
    if i == 1
        Rlow(i) = (R(i) - Rsd(i))/2;
        Rup(i) = R(i) + Rsd(i) + (R(i+1) - Rsd(i+1))/2;
    elseif i == norder
        Rlow(i) = R(i-1) + Rsd(i-1) + (R(i) - Rsd(i))/2;
        Rup(i) = R(i) + Rsd(i);
    else
        Rlow(i) = R(i-1) + Rsd(i-1) + (R(i) - Rsd(i))/2;
        Rup(i) = R(i) + Rsd(i) + (R(i+1) - Rsd(i+1))/2;
    end
end

% To initialize a fractal tree [Olufsen, 2000]
% eta = 1.2;  % area ratio >1
% xi = 2.72;  % Murray's law exponent

% initial inner radius of root vessel; it affects the prunning
% R0 = 0.0036; %m % as from symmetric tree results
R0 = max(R);
% assumed initial thickness of root vessel (as disscussed with Baek)
H0 = 0.07*R0;

% Minimum vesel radius to define terminal vessel (or number of generations)
% We take ~8th generation of PAs (small arteries)
% Rmin = 0.00018; %m 

% Flag to assign terminal outlets BC: 1 - pressure, 0 - flow
PTflag = 0;

disp(['number of generations:', num2str(nlayer)]);
disp(['initial root radius:', num2str(R0*1000),'mm']);
disp(['terminal pressure BC flag:', num2str(PTflag)]);

%----------------- Reference Values for Results Comparisons ---------------
% --- Stiffness metrics reference values (Ehr)
% Stiffness parameter E11*H*Rp0 = 3/4/Lambda for isotropic incompressible
%   - Lambda - distensibility parameter
%   - Rp0 - radius at zero pressure
% (Krenz 2003) Lambda=0.02/mmHg
% (Yen 1990) according to (Qureshi 2014) Lambda=0.012 /mmHg
%   - HG: average lambda from Yen et al (1990) 0.0146 /mmHg
EhrY = 3/(4*0.012)*133.32; % Pa, Yen

% questionable assumptions:
EhrK = 37.5*133.32; % Pa, (Olufsen 2012) Eh/r0=3/4/Lambda, 
EhrQ = 195*133.32;  % Pa, (Qureshi 2014) 

% Ehr from porcine experiments data (hard-codded, Hamid), mmHg
EhrExperiment = [10530.92136,9496.887786,9311.097003,8679.807346,8432.782151]; 
EhrEx_mean = mean(EhrExperiment);
EhrEx_err = EhrExperiment - EhrEx_mean;

%----------------------- Hemodynamics Parameters --------------------------
% Blood properties
mu = 0.0035;    % Pa*s, dynamic viscosity
rho = 1060;     % kg/m^3, blood density
nu = mu/rho;    % m^2/s

% Terminal reflection coefficient
%   0   - matching impedance (Z0=ZT=Zinp), 
%   1   - closed end (ZT>>Z0, e.g. Z0=0, Zinp=0)
%   -1  - open end (Z0>>ZT, e.g. ZT=0)
% [Hollander-2001]: pulmonary circulation is of open-end type reflector
GammaVal = -1;

% BC type
if PTflag==0 % && pruneflag =1 asymm BC
    % Pressure waveform from SymmCase2 results for root vessel
    load Pinlet.dat; 
    %Pinlet = load("Pressure.inp"); %test
    time = Pinlet(:,1);
    press = Pinlet(:,2);
    %press = Pinlet(:,2)*133.32; %test
    pInp_parent =  mean(press);
%     pInp_parent=12*133.32; %Pa
    p_BC = pInp_parent;
    
    % Terminal mean flow is q_parent/number of outlets
    % it will be determined later
    q_parent = 1.1646e-05; %m^3/s 
    
    disp(['mean input flow:',num2str(q_parent*10^6),'cm^3/s']);
    disp(['mean input pressure:',num2str(pInp_parent/133.32),'mmHg']);

else %PTflag==1
    % Flow waveform data from [Zambrano 2018]
    load FlowBZ.dat;
    time = FlowBZ(:,1);
    flow = FlowBZ(:,2); 
%     Start from 4th generation (after RIA in Olufsen 2012) of symm. tree
    flow = flow./2^3; %4th gen from MPA
    
%     data = load("StrahlerGraphTreePinQout.mat"); %from other BC simulations
%     time = data.time;
%     flow = data.qInpTime(:,1); %fisrt segment
    
    q_parent =  mean(flow);
    p_terminal = 8*133.32; % use 8 as in PinQout values at terminal vessels , 10*133.32;     % Pa
    p_BC = p_terminal;
    
    disp(['mean input flow:',num2str(q_parent*10^6),'cm^3/s']);
    disp(['mean terminal pressure:',num2str(p_terminal/133.32),'mmHg']);
    
    
%%PTflag==1, (&& pruneflag =0) symm BC
%     % redefines some parameterst to be as in SymmCase2 (for verification)
%     N_gen = 19;
%     R0 = 0.0055; %m
%     Rmin = 0.00018; %m
%     
%     % Flow waveform data from [Zambrano 2018]
%     load FlowBZ.dat;
%     time = FlowBZ(:,1);
%     flow = FlowBZ(:,2); 
%     % Start from 4th generation (after RIA in Olufsen 2012) of symm. tree
%     flow = flow./2^3; %4th gen from MPA
%     q_parent =  mean(flow);
% 
%     p_terminal = 10*133.32;     % Pa
%     p_BC = p_terminal;
%     
%     disp(['mean input flow:',num2str(q_parent*10^6),'cm^3/s']);
%     disp(['mean terminal pressure:',num2str(p_terminal/133.32),'mmHg']);
end
    
% time steps
Nf = size(time,1);
T = time(Nf) - time(1);
Nt = 125;
% frequency modes, filter higher modes
NumModes = 10; 

%% =================== RUN HOMEOSTATIC OPTIMIZATION =======================
% - Optimization of Mass, Geometry and Material parameters
% - Steady solution for the asymmetric tree
%   - Steady flow and mid-pressure for Mass Optimization
%   - R>Rmin, terminal flow(Newgen) is given
%   - Mt is the mass of Collagen and SMC, Me is the mass of elastin pInpSteady
%   - SMCtoCOL is ratio of content of smc to collagen
%   - Mtotal = Mc + Ms + Me; Mtotal = Mt*(5/4+SMCtoCOL)/(SMCtoCOL+1) + Me
%-------Step 0: Optimization for initial number of generations-------------
disp('=== Start Homeostatic Optimization===');

%update segm,p,q
[Mt,Me,SMCtoCOL,Radius,Length,Table,p_mid,p_inp]=...
   Optimization_GraphTree(q_parent,p_BC,segm,Graph,R,L,Rsd,Lsd,Rlow,Rup); 

disp('End Homeostatic Optimization');

%Table(j,:)=[order, flow, terminal pressure, resistance];
    
% ----- Set arrays for postprocessing---------
% preallocate material and geometry parameters
NumberOfVessels = nsegm;
Thickness = zeros(1,NumberOfVessels);
YoungMod_tt = zeros(1,NumberOfVessels);
YoungMod_zz = zeros(1,NumberOfVessels);
nu_tz = zeros(1,NumberOfVessels);
nu_zt = zeros(1,NumberOfVessels);
StiffMatrix = zeros(2,2,NumberOfVessels);
ksi = zeros(1,NumberOfVessels-1);

% preallocate homeostatic values
ratio = zeros(1,NumberOfVessels);
sigma_h = zeros(1,NumberOfVessels);
shear = zeros(1,NumberOfVessels);
Pmid = zeros(1,NumberOfVessels);

% preallocate steady state solution
qSteady = zeros(1,NumberOfVessels);
pInpSteady = zeros(1,NumberOfVessels);
pTermSteady =zeros(1,NumberOfVessels);
HydRes = zeros(1,NumberOfVessels);
Ehr = zeros(1,NumberOfVessels);

% % assign given inlet mean-flow and outlet terminal pressure
% qSteady(1) = q_parent; 
% pTermSteady(NumberOfVessels) = Table(1,3); 

RzeroP = zeros(1,NumberOfVessels);

for k=1:NumberOfVessels
    
    %------------- Find unstressed radius ---------------------------------
    RzeroP(k)=ZeroP(Me(k),Mt(k),Radius(k),p_mid(k),0);

    %-------------- get parameters and steady solution -------------
    % pInpSteady(k) = pTermSteady(k)+ HydRes(k)*qSteady(k);
    qSteady(k) = Table(k,2);
    pTermSteady(k) = Table(k,3);
    HydRes(k) = Table(k,4);
    pInpSteady(k) = p_inp(k);

    % geometry and material parameters
    % use mid-pressure
    [YoungMod_tt(k),YoungMod_zz(k),nu_tz(k),nu_zt(k),StiffMatrix(:,:,k)] = ...
           YoungMod_2(Me(k),Mt(k),1,1,p_mid(k),SMCtoCOL(k));
       
    Thickness(k) =(1/(SMCtoCOL(k)+1)*(1.25)*Mt(k) + ...
         SMCtoCOL(k)/(SMCtoCOL(k)+1)*Mt(k) + Me(k))/(0.3*rho_w);

    % for output
    Pmid(k) = p_mid(k);
    ratio(k) = Thickness(k)/(2*Radius(k));
    sigma_h(k) = Pmid(k)/(2*ratio(k));          % homeostatic stress
    shear(k) = 4*mu*qSteady(k)/(pi*Radius(k)^3);% homeostatic shear stress
    Ehr(k) = YoungMod_tt(k)*Thickness(k)/RzeroP(k);
end  

%% Get eta, gamma, xi

%get connectivity at segments at end node:
% either 0 (outlet) or 2(bifurcation)
connectOfSegm  = zeros(nsegm,1);
outNodeID = find(Gdegree==0);
bifNodeID = find(Gdegree==2); 

numOfOut = size(outNodeID,1);
numOfBif = size(bifNodeID,1);

for i=1:numOfOut
    connectOfSegm(enode(:) == outNodeID(i)) = 0;
end

for i=1:numOfBif
    connectOfSegm(enode(:) == bifNodeID(i)) = 2;
end

% get distributions of eta, gamma and xi accross the tree
%eta_compute;
Compute_eta_gamma_xi_GraphTree

%% ======================= RUN PULSATILE FLOW =============================
% Deformable wall Womersley solution (Womersley 1957)
 
disp('Start Pulsatile Solutions');

% for longitudinally constrained vessel (tethered)
Ctt(1:NumberOfVessels)=StiffMatrix(1,1,:);

if PTflag==0 %asymm BC
    [PnInp] = FilterAndFFT(press,NumModes,Nf-1);
    
    [InpPressureTime,TermPressureTime,qInpTime,qTermTime,...
    cn,Alpha,zinp,zterm,zchar,InpImpedance]...
    = WomersleySolution_BCPin_GraphTree...
      (T,PnInp,pInpSteady,Table,Radius,Length,Thickness,Ctt,...
       Nt,NumModes,GammaVal,segm,Graph);
  
else % PTflag==1
%     disp('Stop here: pulsatile solution at this BC is under construction');
    [QnInp] = FilterAndFFT(flow,NumModes,Nf-1);
    
    [InpPressureTime,TermPressureTime,qInpTime,qTermTime,...
     cn,Alpha,zinp,zterm,zchar,InpImpedance]...
     = WomersleySolution_BCQin_GraphTree...
       (T,QnInp,pInpSteady,Table,Radius,Length,Thickness,Ctt,...
       Nt,NumModes,GammaVal,segm,Graph);
end

% %check input root pressure, almost Okay!
% figure(1)
% plot(InpPressureTime(:,1),'r'); hold on;
% plot(press,'b'); hold off;

% real wave speed for 1st frequency (n=2)
c_R2 = zeros(1,NumberOfVessels); c_R10 = zeros(1,NumberOfVessels);
for k=1:NumberOfVessels
    c_R2(k) = 1.0/real(cn(2,k)^(-1));
    c_R10(k) = 1.0/real(cn(10,k)^(-1));
end

% Moens-Korteweg pulse wave velocity
c0_MK = zeros(1,NumberOfVessels); 
for k=1:NumberOfVessels
    c0_MK(k) = sqrt(Thickness(k)*YoungMod_tt(k)/(2*rho_w*Radius(k)));
end

% Estimate delta parameter 
delta1 = zeros(1,NumberOfVessels);delta9 = zeros(1,NumberOfVessels); 
deltaMK = zeros(1,NumberOfVessels);
for k=1:NumberOfVessels
    maxqpulse = max(qInpTime(:,k)) - qSteady(k);
    delta1(k) = maxqpulse/(pi*Radius(k)^2*c_R2(k));
    delta9(k) = maxqpulse/(pi*Radius(k)^2*c_R10(k));
    deltaMK(k) = maxqpulse/(pi*Radius(k)^2*c0_MK(k));
end

% Warning: the more assymetric tree the more complient shortest path and
% delta(short) exceeds 1 which means that pulsatile flow there cannot be
% computed using Womerslety (linear) theory (need to use nonlinear 1D th)
[a,i]=max(delta1);
disp(['NOTE: maximum delta(',num2str(i),')=',num2str(a)]);
if (max(delta1)>1)
    disp(['WARNING:delta(',num2str(i),')=',num2str(a),...
        '>1 - violates Womersley assumptions!'])
end

disp('End Pulsatile Solutions');

%% ========================== FIGURES =====================================
% Output for the presentatoins and paper

PlottingResults_GraphTree;

% pause
% close all

%% ====================== SAVE RESULTS ====================================
% if PTflag==0 %asymm BC
%     save('Strahler_GraphTree'); %.mat workspace
% else % PTflag==1 
%     %
% end                
% 
fileID = fopen('HomeostaticRandL.dat','w');
for i=1:nsegm
    if sum(ismember(longpathsegm,i))
        fprintf(fileID,'%12.8f %12.8f %d\n',Radius(i),Length(i),1);
    elseif sum(ismember(shortpathsegm,i))
        fprintf(fileID,'%12.8f %12.8f %d\n',Radius(i),Length(i),2);
    else
        fprintf(fileID,'%12.8f %12.8f %d\n',Radius(i),Length(i),0);
    end
end
fclose(fileID);
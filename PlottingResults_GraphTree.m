%% ======================= SUBPROGRAMM: GET PLOTS =========================
% Plotting some results for Strahler ordering graph
% Call:
%    - x_location_GrapphTree
% -----
% By Hamid and Vasilina, 2018
%   - Updated on October 26 2020
%==========================================================================
%%
%global PTflag  pruneflag % for the testing symm case


%% get location along the branches (pathway)
[x_loc] = x_location_GrapphTree(layer,Graph,segm,Length);

%% ------------------- graph and path----------
% Get longest path from given graph segments
% get last layer nlayer
longpathnodes = zeros(nlayer,2);
longpathsegm = zeros(nlayer,1);
lnode = find(layer(:)==1);
longpathnodes(nlayer,2)=lnode(randi(size(lnode,1),1)); %end node
longpathsegm(nlayer) = find(enode(:)==longpathnodes(nlayer,2));%inedges(G,longpathnodes(nlayer,2));
longpathnodes(nlayer,1)=snode(longpathsegm(nlayer));
for i = nlayer-1:-1:2 %reverse order
    longpathnodes(i,2)=longpathnodes(i+1,1);
    longpathsegm(i) = find(enode(:)==longpathnodes(i,2));%inedges(G,longpathnodes(i,2));
    longpathnodes(i,1)=snode(longpathsegm(i));
end
longpathnodes(1,:)=[];
longpathsegm(1)=[];

% find one of the short paths
% find outlet node
for i=1:nlayer %reverse order
    layernnodes = find(layer(:)==i);
    for j=1:size(layernnodes,1)
        if Gdegree(layernnodes(j))==0
            shrnode = layernnodes(j);
            sgen = nlayer - i +1; %i;
            break
        end
    end
end
shortpathnodes = zeros(sgen,2);
shortpathsegm = zeros(sgen,1);
shortpathnodes(sgen,2)=shrnode; %end node
shortpathsegm(sgen) = find(enode(:)==shrnode); %inedges(G,shrnode);
shortpathnodes(sgen,1)=snode(shortpathsegm(sgen));
if sgen>2
    for i = sgen-1:-1:2
        shortpathnodes(i,2)=shortpathnodes(i+1,1);
        shortpathsegm(i) = find(enode(:)==shortpathnodes(i,2)); %inedges(G,shortpathnodes(i,2));
        shortpathnodes(i,1)=snode(shortpathsegm(i));
    end
end
shortpathnodes(1,:)=[];
shortpathsegm(1)=[];

h6 = figure;
p = plot(Graph,'EdgeColor','k','MarkerSize',0.01);
layout(p,'layered');
highlight(p,longpathnodes(:,1),longpathnodes(:,2),'EdgeColor','r','LineWidth',2);
highlight(p,shortpathnodes(:,1),shortpathnodes(:,2),'EdgeColor','b','LineWidth',2);

%% Omptimization results

%Terminal Pressure along the short and long paths
% h7 = figure; nid(1)=1;
% %include input pressure
% Plong(1) = pInpSteady(1); Pshort(1) = pInpSteady(1);
% Plong(2:size(longpathsegm,1)+1) = pTermSteady(longpathsegm);
% Pshort(2:size(shortpathsegm,1)+1) = pTermSteady(shortpathsegm);
% plot(0:1:size(longpathsegm,1),Plong./133.32,'ro--','LineWidth',2); hold on
% plot(0:1:size(shortpathsegm,1),Pshort./133.32,'bs--','LineWidth',2); 
% xlabel('Generation','FontSize',14);
% ylabel('Steady Terminal Pressure (mmHg)','FontSize',14);
% legend({'Longest path','Shortest path'},'FontSize',14);
% set(gca,'xtick',0:1:size(longpathsegm,1));
% box on; hold off;

% %Steady Flow along the short and long paths
% h8 = figure; nid(1)=1;
% plot(1:size(longpathsegm,1),qSteady(longpathsegm).*10^6,'ro--','LineWidth',2); hold on
% plot(1:size(shortpathsegm,1),qSteady(shortpathsegm).*10^6,'bs--','LineWidth',2);
% xlabel('Generation','FontSize',14);
% ylabel('Steady Inflow (ml/s)','FontSize',14);
% legend({'Longest path','Shortest path'},'FontSize',14);
% set(gca,'xtick',1:size(longpathsegm,1));
% box on; hold off;

%Terminal Pressure vs Distance
h9 = figure; hold on;
scatter(x_loc*100,Table(:,3)/133.32,'.','MarkerEdgeColor',[238 31 35]/255);
xlabel({'X (cm)'},'FontSize',14);
ylabel({'Steady Terminal Pressure (mmHg)'},'FontSize',14);
plot(x_loc(longpathsegm)*100,Table(longpathsegm,3)./133.32,'o--','LineWidth',2,'Color',[62 88 166]/255);
plot(x_loc(shortpathsegm)*100,Table(shortpathsegm,3)./133.32,'s--','LineWidth',2,'Color',[37 158 49]/255);
box on;
set(gca,'FontSize',14)
print(h9, '-dmeta', ['PressX','.emf']);
saveas(h9, ['PressX','.png']);

%log(Flow) vs Distance
h11 = figure; hold on;
scatter(x_loc*100,log10(Table(:,2)*10^6),'.','MarkerEdgeColor',[238 31 35]/255);
xlabel({'X (cm)'},'FontSize',14);
ylabel({' Log of steady flow (cm^3/s)'},'FontSize',14);
plot(x_loc(longpathsegm)*100,log10(Table(longpathsegm,2).*10^6),'o--','LineWidth',2,'Color',[62 88 166]/255);
plot(x_loc(shortpathsegm)*100,log10(Table(shortpathsegm,2).*10^6),'s--','LineWidth',2,'Color',[37 158 49]/255);
box on;
set(gca,'FontSize',14)
print(h11, '-dmeta', ['FlowX','.emf']);
saveas(h11, ['FlowX','.png']);

%WSS vs Distance
h12 = figure; hold on;
scatter(x_loc*100,shear,'.','MarkerEdgeColor',[238 31 35]/255);
xlabel({'X (cm)'},'FontSize',14);
ylabel({'\tau (Pa)'},'FontSize',14);
plot(x_loc(longpathsegm)*100,shear(longpathsegm),'o--','LineWidth',2,'Color',[62 88 166]/255);
plot(x_loc(shortpathsegm)*100,shear(shortpathsegm),'s--','LineWidth',2,'Color',[37 158 49]/255);
box on;
set(gca,'FontSize',14)
print(h12, '-dmeta', ['WSSX','.emf']);
saveas(h12, ['WSSX','.png']);

%Eh/r0 vs Distance
h13 = figure; hold on;
scatter(x_loc*100,Ehr/1000,'.','MarkerEdgeColor',[238 31 35]/255);
xlabel({'X (cm)'},'FontSize',14); 
plot(x_loc(longpathsegm)*100,Ehr(longpathsegm)/1000,'o--','LineWidth',2,'Color',[62 88 166]/255);
plot(x_loc(shortpathsegm)*100,Ehr(shortpathsegm)/1000,'s--','LineWidth',2,'Color',[37 158 49]/255);
XLOC = min(x_loc):0.01:max(x_loc);
EhrKall(1:length(XLOC))=EhrK; EhrYall(1:length(XLOC))=EhrY;
EhrQall(1:length(XLOC))=EhrQ;
plot(XLOC*100, EhrKall./1000,'--k','LineWidth',2);
plot(XLOC*100, EhrYall./1000,'-.k','LineWidth',2);
legend({'Vessels','Shortest path','Longest path',...
 'Eh/R0 [Krenz-2003]','Eh/R0 [Yen-1990]'},'Location','SouthEast','FontSize',14);
xlabel('X (cm)','FontSize',14);ylabel('E_{\theta\theta}h/R0 (kPa)','FontSize',14);
ylim([0 12]);
box on; hold off;
set(gca,'FontSize',14)
print(h13, '-dmeta', ['EHRX','.emf']);
saveas(h13, ['EHRX','.png']);

%xi along short and long path, bifurcation at the end of segment
h14 = figure; hold on;
set(gca,'Units','normalized','FontUnits','points','FontWeight','normal',...
    'FontSize',12,'FontName','Calibri');
plot(1:size(longpathsegm,1),xi_segm(longpathsegm),'ro','LineWidth',2);
plot(1:size(shortpathsegm,1),xi_segm(shortpathsegm),'bs','LineWidth',2);
legend('Longest path','Shortest path','location','SouthEast');
xlabel('Generation','FontSize',14);
ylabel('xi','FontSize',14);
ylim([2 3.5]);
box on; hold off;
set(gca,'FontSize',14)
print(h14, '-dmeta', ['Xi','.emf']);
saveas(h14, ['Xi','.png']);

%Youngs modulus vs Distance
% h15 = figure; hold on;
% scatter(x_loc*100,YoungMod_tt./1000,'.');
% %xlabel({'2R (cm)'},'FontSize',14); 
% plot(x_loc(longpathsegm).*100,YoungMod_tt(longpathsegm)./1000,'ro--','LineWidth',2);
% plot(x_loc(shortpathsegm).*100,YoungMod_tt(shortpathsegm)./1000,'bs--','LineWidth',2);
% legend({'Optimization Results','Shortest path','Longest path'},...
%        'Location','NorthEast','FontSize',14);
% xlabel('X (cm)','FontSize',14);ylabel('E_{\theta\theta}(kPa)','FontSize',14);
% box on; hold off;

%D for long and short paths
% h16 = figure; hold on;
% scatter(x_loc*100,Radius.*200,'.');
% plot(x_loc(longpathsegm,1).*100,Radius(longpathsegm).*200,'ro--','LineWidth',2);
% plot(x_loc(shortpathsegm,1).*100,Radius(shortpathsegm).*200,'bs--','LineWidth',2);
% legend({'Optimization Results','Shortest path','Longest path'},...
%        'Location','NorthEast','FontSize',14);
% xlabel('X (cm)','FontSize',14);ylabel('D(cm)','FontSize',14);
% box on; hold off;

%L for long and short paths
% h17 = figure; hold on;
% scatter(x_loc*100,Length.*100,'.');
% plot(x_loc(longpathsegm,1).*100,Length(longpathsegm).*100,'ro--','LineWidth',2);
% plot(x_loc(shortpathsegm,1).*100,Length(shortpathsegm).*100,'bs--','LineWidth',2);
% legend({'Optimization Results','Shortest path','Longest path'},...
%        'Location','NorthEast','FontSize',14);
% xlabel('X (cm)','FontSize',14);ylabel('L(cm)','FontSize',14);
% box on; hold off;

% L(R)
% h18 = figure; hold on;
% scatter(Radius*100,Length.*100,'.');
% plot(Radius(longpathsegm,1).*100,Length(longpathsegm).*100,'ro--','LineWidth',2);
% plot(Radius(shortpathsegm,1).*100,Length(shortpathsegm).*100,'bs--','LineWidth',2);
% legend({'Optimization Results','Shortest path','Longest path'},...
%        'Location','SouthEast','FontSize',14);
% xlabel('R(cm)','FontSize',14);ylabel('L(cm)','FontSize',14);
% box on; hold off;

% %h/D for long and short paths
% h19 = figure; hold on;
% scatter(2.*Radius*.100,transpose(Thickness)./(2.*Radius),'.');
% plot(2.*Radius(longpathsegm).*100,ratio(longpathsegm),'ro--','LineWidth',2);
% plot(2.*Radius(shortpathsegm).*100,ratio(shortpathsegm),'bs--','LineWidth',2);
% legend({'Optimization Results','Shortest path','Longest path'},...
%        'Location','NorthEast','FontSize',14);
% xlabel('D (cm)','FontSize',14);ylabel('h/D','FontSize',14);
% set(gca, 'XDir','reverse');
% box on; hold off;

%% Pulsatile Hemodynamics

% %Delta
% q1 = figure; set(gca,'xscale'); 
% plot(1:size(longpathsegm,1),delta1(longpathsegm),'ro--','LineWidth',2);
% hold on;
% plot(1:size(shortpathsegm,1),delta1(shortpathsegm),'bs--','LineWidth',2);
% legend({'Longest path','Shortest path'},'FontSize',14);%'Location','NorthEast','FontSize',14);
% xlabel('Generation','FontSize',14);
% ylabel('{\delta} = maxFlowVel/WaveSpeed','FontSize',14);
% set(gca,'xtick',1:1:size(longpathsegm,1));
% box on; hold off; 

%PWV
q10=figure; set(gca,'xscale'); hold on;
% Banks (1978): Moens-Korteweg PVW for 20-30yo humans c0=2.24m/s
% Milnor (1969): c0=1.68m/s;
PwvB = 2.24; PwvM = 1.68;
plot(x_loc(longpathsegm)*100, c_R2(longpathsegm),'ro--','LineWidth',2,'Color',[62 88 166]/255);
plot(x_loc(shortpathsegm)*100, c_R2(shortpathsegm),'s--','LineWidth',2,'Color',[37 158 49]/255);
plot(x_loc(longpathsegm)*100, c0_MK(longpathsegm),'.-k','LineWidth',2);
plot(x_loc(longpathsegm(1))*100 ,PwvB,'dk','LineWidth',2);
plot(x_loc(longpathsegm(1))*100 ,PwvM,'^k','LineWidth',2);
legend({'PWV, long path','PWV, short path',...
        'c^{MK} with E_{\theta\theta}','[Bank-1978]',...
        '[Milnor-1969]'},'Location','NorthEast','FontSize',10);
xlabel('X (cm)','FontSize',14);
ylabel('Pulse wave velocity (m/s)','FontSize',14);
axis([0 size(longpathsegm,1) 0 3.5]);
box on;
set(gca,'FontSize',14)
print(q10, '-dmeta', ['PWVX','.emf']);
saveas(q10, ['PWVX','.png']);

% Pulsatile P and Q along the paths
% time points
t = linspace(0,T-T/Nt,Nt);

%long path
q2 = figure; col=jet(nlayer); hold on;
for k=1:nlayer-1
    plot(t',qInpTime(:,longpathsegm(k)).*10^6,'color',col(k,:),'LineWidth',1);
    if k==nlayer-1
         plot(t',qTermTime(:,longpathsegm(k)).*10^6,'color',col(nlayer,:),'LineWidth',1);
    end
    hold on;
end
% title({'Input Flow along the long path'},'FontSize',14);
xlabel('Time (s)','FontSize',14); ylabel('Flow (cm^3/s)','FontSize',14); 
box on; hold off;
set(gca,'FontSize',14)
print(q2, '-dmeta', ['FlowPulsLong','.emf']);
saveas(q2, ['FlowPulsLong','.png']);

q3 = figure; col=jet(nlayer); hold on;
for k=1:nlayer-1
    plot(t',InpPressureTime(:,longpathsegm(k))./133.32,'color',col(k,:),'LineWidth',1);
    if k==nlayer-1
         plot(t',TermPressureTime(:,longpathsegm(k))./133.32,'color',col(nlayer,:),'LineWidth',1);
    end
    hold on
end
hold off;
% title({'Terminal Pressure along the long path'},'FontSize',14);
xlabel('Time (s)','FontSize',14); ylabel('Pressure (mmHg)','FontSize',14); 
box on;
set(gca,'FontSize',14)
print(q3, '-dmeta', ['PressPulsLong','.emf']);
saveas(q3, ['PressPulsLong','.png']);

%short path
q4 = figure; col=jet(sgen); hold on;
for k=1:sgen-1
    plot(t',qInpTime(:,shortpathsegm(k)).*10^6,'color',col(k,:),'LineWidth',1);
    if k==nlayer-1
         plot(t',qTermTime(:,shortpathsegm(k)).*10^6,'color',col(sgen,:),'LineWidth',1);
    end
    hold on;
end
hold off;
% title({'Input Flow along the short path'},'FontSize',14);
xlabel('Time (s)','FontSize',14); ylabel('Flow (cm^3/s)','FontSize',14); 
box on;
set(gca,'FontSize',14)
print(q4, '-dmeta', ['FlowPulsShort','.emf']);
saveas(q4, ['FlowPulsShort','.png']);

q5 = figure; col=jet(sgen); hold on;
for k=1:sgen-1
    plot(t',InpPressureTime(:,shortpathsegm(k))./133.32,'color',col(k,:),'LineWidth',1);
    if k==nlayer-1
         plot(t',TermPressureTime(:,shortpathsegm(k))./133.32,'color',col(sgen,:),'LineWidth',1);
    end
    hold on;
end
hold off;
% title({'Term Pressure along the short path'},'FontSize',14);
xlabel('Time (s)','FontSize',14); ylabel('Pressure (mmHg)','FontSize',14); 
box on;
set(gca,'FontSize',14)
print(q5, '-dmeta', ['PressPulsShort','.emf']);
saveas(q5, ['PressPulsShort','.png']);

%% -----------
% %% Hemodynamics
% 
% %PWV
% q1=figure; set(gca,'xscale'); hold on;
% % Banks (1978): Moens-Korteweg PVW for 20-30yo humans c0=2.24m/s
% % Milnor (1969): c0=1.68m/s;
% PwvB = 2.24; PwvM = 1.68;
% plot(Table(y,1), c_R2(y),'o-','LineWidth',2); 
% plot(Table(x,1), c_R2(x),'^b--','LineWidth',2); 
% plot(Table(y,1), c0_MK(y),'.-k','LineWidth',2);
% plot(1 ,PwvB,'*','LineWidth',2);
% plot(1 ,PwvM,'*','LineWidth',2);
% legend({'PWV(1st mode), longest path','PWV(1st mode), shortest path',...
%         'c^{MK} with E_{\theta\theta}','[Bank-1978]',...
%         '[Milnor-1969]'},'Location','NorthEast','FontSize',14);
% xlabel('gen.','FontSize',14);
% ylabel('Pulse wave velocity (m/s)','FontSize',14);
% axis([0 max(Table(y,1)) 0 3]);
% box on;
% 
% %Delta
% q2=figure; set(gca,'xscale'); hold on;
% plot(Table(y,1),delta1(y),'o-','LineWidth',2);
% plot(Table(x,1),delta1(x),'^b--','LineWidth',2);
% legend({'{\delta}(1st mode), longest path',...
%       '{\delta}(1st mode), shortest path'},...
%      'Location','NorthEast','FontSize',14);
% xlabel('gen.','FontSize',14);
% ylabel('{\delta}=max(v^f_{z})/c','FontSize',14);
% axis([0 max(Table(y,1)) 0 1]);
% box on;
% 
% % time points
% t = linspace(0,T-T/Nt,Nt);
% 
% %Flow and Pressure along the shortest path
% q3 = figure; 
% col=jet(size(x,2)); 
% 
% subplot(2,1,1);
% hold on;
% for k=1:size(x,2)
%     if k==1 
%         plot(t,qInpTime(:,x(k))*10^6,'color',col(k,:),'LineWidth',2);
%     else
%         plot(t,qInpTime(:,x(k))*10^6,'color',col(k,:));
%     end
%     hold on
% end
% title(['Input Flow along shortest path'],'FontSize',14);
% xlabel('Time (s)','FontSize',14); ylabel('Flow (cm^3/s)','FontSize',14); 
% set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
% axis([0 T -10 55]);
% 
% subplot(2,1,2);
% hold on;
% for k=1:size(x,2)
%     if k==1
%         plot(t,TermPressureTime(:,x(k))/133.32,'color',col(k,:),'LineWidth',2);    
%     else
%         plot(t,TermPressureTime(:,x(k))/133.32,'color',col(k,:));
%     end
%     hold on
% end
% hold off
% title(['Terminal Pressure along shortest path'],'FontSize',14);
% xlabel('Time (s)','FontSize',14);  ylabel('Presure (mmHg)','FontSize',14);   
% set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
% axis([0 T 5 25]);
% box on;
% 
% %Flow and Pressure along the longest path
% q4 = figure; 
% col=jet(size(y,1)); 
% 
% subplot(2,1,1);
% hold on;
% for k=1:size(y,1)
%     if k==1 
%         plot(t,qInpTime(:,y(k))*10^6,'color',col(k,:),'LineWidth',2);
%     else
%         plot(t,qInpTime(:,y(k))*10^6,'color',col(k,:));
%     end
%     hold on
% end
% title('Input Flow along longest path','FontSize',14);
% xlabel('Time (s)','FontSize',14); ylabel('Flow (cm^3/s)','FontSize',14); 
% set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
% axis([0 T -10 55]);
% 
% subplot(2,1,2);
% hold on;
% for k=1:size(y,1)
%     if k==1
%         plot(t,TermPressureTime(:,y(k))/133.32,'color',col(k,:),'LineWidth',2);    
%     else
%         plot(t,TermPressureTime(:,y(k))/133.32,'color',col(k,:));
%     end
%     hold on
% end
% hold off
% title('Terminal Pressure along longest path','FontSize',14);
% xlabel('Time (s)','FontSize',14);  ylabel('Presure (mmHg)','FontSize',14);   
% set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
% axis([0 T 5 25]);
% box on;
% 
% %Root pressure
% q5=figure;hold on
% % plot(t,press/133.32,'-','LineWidth',2); % test
% plot(t,InpPressureTime(:,1)/133.32,'-','LineWidth',2);
% plot(t,TermPressureTime(:,1)/133.32,'-','LineWidth',2);
% % legend({'Input p as BC','Input p obtained at 1st gen.','Terminal p at 1st gen.'},'FontSize',14);
% legend({'Input p at 1st gen.','Terminal p at 1st gen.'},'FontSize',14);
% xlabel('Time (s)','FontSize',14);  ylabel('Presure (mmHg)','FontSize',14);   
% set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
% axis([0 T 5 25]);
% box on;
% 
% %Root flow
% q6=figure;hold on
% plot(t,qInpTime(:,1)*10^6,'-','LineWidth',2);
% plot(t,qTermTime(:,1)*10^6,'-','LineWidth',2);
% legend({'Input q at 1st gen.','Terminal q at 1st gen.'},'FontSize',14);
% xlabel('Time (s)','FontSize',14);  ylabel('Flow (cm^3/s)','FontSize',14);   
% set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
% axis([0 T -10 55]);
% box on;
% 
% %Impedance in time domain
% q7=figure;hold on
% plot(t,zinp(:,1)./(133.32*10^6),'-','LineWidth',2);  
% plot(t,zterm(:,1)./(133.32*10^6),'-.','LineWidth',1);  
% plot(t,zchar(:,1)./(133.32*10^6),'--','LineWidth',1);  
% set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
% legend({'Z_{inp}','Z_{term}','Z_c-(pulsatile)'},'FontSize',14);
% xlabel('Time (s)','FontSize',14); ylabel('Impedance(mmHg*s/cm^3)','FontSize',14);
% axis([0 T -5 10]);
% box on;
% 
% %Impedance in frequency domain
% q8=figure;
% kk = 0:1:NumModes-1;
% 
% subplot(1,2,2);hold on;
% plot(kk,angle(InpImpedance(kk+1,1)),'o','LineWidth',2);  
% xlabel('Frequency modes','FontSize',14); ylabel('Phase Angle','FontSize',14);
% 
% subplot(1,2,1); hold on;
% plot(kk,abs(InpImpedance(kk+1,1))./(133.32*10^6),'o','LineWidth',2);
% set(gca,'xtick'); xlabel('Frequency modes','FontSize',14);
% ylabel('Modulus (mmHg*s/cm^3)','FontSize',14);
% legend('Asymm.');
% box on;
% 
% %% Pulmonary Vascular Resistance (Wood Units = mmHG Litre/min)
% % longest path, y
% pl = pInpSteady(1,y); ql = qSteady(1,y);
% PVR(:) = (pl(:)-pl(end))./ql(:);
% PVR = PVR./(133.32*10^6*60); % to mmHg min/L
% % shortest path, x
% ps = pInpSteady(1,x); qs = qSteady(1,x);
% PVRs(:) = (ps(:)-ps(end))./qs(:);
% PVRs = PVRs./(133.32*10^6*60); % to mmHg min/L
% 
% q9=figure; hold on
% plot(PVR,'ro'); plot(PVRs,'bs');
% legend('longest path','shortest path','Location','northwest');
% xlabel('gen.'); ylabel('PVR (mmHg min/cm^3)');
% box on;
% hold off;
% 
% %% Verify Asymm (symm case) results with Var
% 
% if (PTflag==1 && pruneflag==0)
%     Var=load('Var.mat');
% 
%     %Root pressure
%     q9=figure;hold on
%     plot(t,InpPressureTime(:,1)/133.32,'-','LineWidth',2);
%     plot(t,TermPressureTime(:,1)/133.32,'-','LineWidth',2);
%     plot(t,Var.pInpTime(:,1)/133.32,'o','LineWidth',2);
%     plot(t,Var.pTermTime(:,1)/133.32,'o','LineWidth',2);
%     legend({'Input Asymm.p at 1st gen.','Terminal Asymm.p at 1st gen.',...
%         'Input Var.p at 1st gen.','Terminal Var.p at 1st gen.'},'FontSize',14);
%     xlabel('Time (s)','FontSize',14);  ylabel('Presure (mmHg)','FontSize',14);   
%     set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
%     % axis([0 T 5 25]);
%     box on;
% 
%     %Root flow
%     q10=figure;hold on
%     plot(t,qInpTime(:,1)*10^6,'-','LineWidth',2);
%     plot(t,qTermTime(:,1)*10^6,'-','LineWidth',2);
%     plot(t,Var.qInpTime(:,1)*10^6,'o','LineWidth',2);
%     plot(t,Var.qTermTime(:,1)*10^6,'o','LineWidth',2);
%     legend({'Input Asymm.q at 1st gen.','Terminal Asymm.q at 1st gen.',...
%         'Input Var.q at 1st gen.','Terminal Var.q at 1st gen.'},'FontSize',14);
%     xlabel('Time (s)','FontSize',14);  ylabel('Flow (cm^3/s)','FontSize',14);   
%     set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
%     % axis([0 T -10 55]);
%     box on;
% 
%     %Impedance in time domain
%     q11=figure;hold on
%     plot(t,zinp(:,1)./(133.32*10^6),'-','LineWidth',2);  
%     plot(t,zterm(:,1)./(133.32*10^6),'-.','LineWidth',1);  
%     plot(t,zchar(:,1)./(133.32*10^6),'--','LineWidth',1);   %Ok
%     plot(t,Var.zinp_total(:,1)./(133.32*10^6),'o','LineWidth',2);  
%     plot(t,Var.zterm_total(:,1)./(133.32*10^6),'*','LineWidth',1);  
%     plot(t,Var.zchar(:,1)./(133.32*10^6),'^','LineWidth',1);  %Ok
%     set(gca,'xtick',[0 T/5 2*T/5 3*T/5 4*T/5 T]);
%     legend({'Asymm.Z_{inp}','Asymm.Z_{term}','Asymm.Z_c-(pulsatile)',...
%             'Var.Z_{inp}','Var.Z_{term}','Var.Z_c-(pulsatile)'},'FontSize',14);
%     xlabel('Time (s)','FontSize',14); ylabel('Impedance(mmHg*s/cm^3)','FontSize',14);
%     axis([0 T -5 10]);
%     box on;
% 
%     %Root steady pressure
%     disp({'root steady inp p:','asymm:',pInpSteady(1,1)/133.32,'symm.var',Var.pInpSteady(1,1)/133.32})
%     disp({'root steady term p:','asymm',pTermSteady(1,1)/133.32,'symm.var',Var.pTermSteady(1,1)/133.32})
% end
movefile('*.emf','.\Figs\');
movefile('*.png','.\Figs\');
function output = main_SW(theta,y0,tspan,ver,mutant)
% close all
%1/24/2021
% clear all
global T_e1
T_e1=35;
if isempty(theta)
load('T_6.mat');%load parameters
end
if isempty(y0)
%load y0
 load('y0_74.mat');%SW IC - first cell cycle
end

global p;
parameters(1,theta,ver,mutant);%all parameters



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Integration parameters
t0  =  0;       % Start time
tf  = tspan;%120;%150;     % End time  %120min Z-ring closed

options  =  odeset('Events',@podJ_event,'RelTol',1e-4,'AbsTol',1e-6);

tout = t0;
y0 = y0.';
yout = y0;

teout  =  [];
yeout  =  [];
ieout  =  [];

while t0<tf
%  [t,y,te,ye,ie] = ode15s(@ODE_CtrA,[t0 tf],y0,options);
[t,y,te,ye,ie] = ode15s(@ODE_DivK,[t0 tf],y0,options);
 nt = length(t);
%     
 tout = [tout;t(2:nt)]; 
 yout = [yout;y(2:nt,:)]; 
 teout  =  [teout;te];  
 yeout  =  [yeout;ye]; ieout  =  [ieout;ie];
    y0  =  y(nt,:);
     if isscalar(ie)  ==  0
        ie  =  0;
     end
    if ie  ==  1%DNA replication initiates
        T_e1 = te;
%         disp(sprintf('T_e1= %8.5f',T_ei))
    elseif ie == 2%Fork passes ctrA
y0(70)=1;%SctrA changed from 0 to 1
    elseif ie==3
y0(71)=1;%SpleC
    elseif ie==4
        y0(72)=1;%SperP
    elseif ie ==5
        y0(69)=1;%SpodJ
         elseif ie==6
       y0(69:72)=0;
    end
    t0 = t(nt);
    if t0 >= tf
        break;
    end
end


%%
    
 %% construct grid for plotting
L=length(yout(:,1));
M(:,4)=yout(:,74);
M(:,5)=yout(:,73)+yout(:,74);
M(:,3)=zeros(1,L);
M(:,1)=-yout(:,73)-yout(:,74);
M(:,2)=-yout(:,74);
    

PodJ(:,1:4) = yout(:,1:4);%PodJm
PodJL(:,1:4)=yout(:,5:8);%PodJL polymer
PodJS(:,1:4)=yout(:,9:12);

% SpmX(:,1:4)=yout(:,13:16)+yout(:,17:20);%SpmX total
SpmXm(:,1:4)=yout(:,13:16);
SpmXp(:,1:4)=yout(:,17:20);

PopZm(:,1:4)=yout(:,21:24);
PopZp(:,1:4)=yout(:,25:28);

CtrA(:,1:4)=yout(:,29:32);
CtrAP(:,1:4)=yout(:,33:36);

PleCf(:,1:4)=yout(:,37:40);
PleCb(:,1:4)=yout(:,41:44);

DivJf(:,1:4)=yout(:,45:48);
DivJb(:,1:4)=yout(:,49:52);

DivK(:,1:4)=yout(:,53:56);
DivKPT(:,1:4)=yout(:,57:60)+yout(:,61:64);

PerP(:,1:4)=yout(:,65:68);
%%%%%


PodJ = fliplr(PodJ);
PodJL=fliplr(PodJL);
PodJS=fliplr(PodJS);
SpmXm = fliplr(SpmXm);
SpmXp = fliplr(SpmXp);
PopZm = fliplr(PopZm);
PopZp = fliplr(PopZp);
CtrA = fliplr(CtrA);
CtrAP = fliplr(CtrAP);
PleCf = fliplr(PleCf);
PleCb = fliplr(PleCb);
DivJf = fliplr(DivJf);
DivJb = fliplr(DivJb);
DivK = fliplr(DivK);
DivKPT = fliplr(DivKPT);
PerP=fliplr(PerP);



PodJ = PodJ.';
M = M.';
PodJL = PodJL.';
PodJS = PodJS.';
PerP=PerP.';
SpmXm=SpmXm.';
SpmXp=SpmXp.';
PopZm=PopZm.';
PopZp=PopZp.';
CtrA=CtrA.';
CtrAP=CtrAP.';
PleCf=PleCf.';
PleCb=PleCb.';
DivJf=DivJf.';
DivJb=DivJb.';
DivK=DivK.';
DivKPT=DivKPT.';
PerP=PerP.';

output.PodJ = PodJ;
output.PodJL = PodJL;
output.PodJS = PodJS;
output.SpmXm = SpmXm;
output.SpmXp = SpmXp;
output.PopZm = PopZm;
output.PopZp = PopZp;
output.CtrA= CtrA;
output.CtrAP= CtrAP;
output.PleCf = PleCf;
output.PleCb= PleCb;
output.DivJf = DivJf;
output.DivJb= DivJb;
output.DivK= DivK;
output.DivKPT= DivKPT;
output.PerP = PerP;
output.time = tout;
output.grid = M;
 output.yout = yout;



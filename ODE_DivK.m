%1/21/2021
function dydt = ODE_DivK(t,y)

global p;
   
%% Equations

dydt=zeros(74,1);
m=0.01;%0.15;
mplec=0.15;%0.15;
%% PodJL_m
m_podj=((1-m)*y(69)+m);
% Bin 1
dydt(1) = 0*p.syn_podJ*m_podj+0*p.syn_podJ2*p.Ji_PodJCtrA^p.n_PodJCtrA/(p.Ji_PodJCtrA^p.n_PodJCtrA+y(1+32)^p.n_PodJCtrA)- (p.deg_podJ1+p.deg_podJm*y(1+64))*y(1)... % synthesis & degradation*PerP (65-68); y(61)=SpodJ 
    - p.dnv_podJ*y(1)...% denovo polymerization; SpmX (13-16&17-20)
    - p.aut1_podJ/(1+p.alpha_PodJSpmX*(y(1+12)+y(1+16)))*y(1+4)^p.podj*y(1)...% autocatalytic polymerization at poles
+ p.depol_podJ*y(1+4)...%
    + 4*p.D_podJm*(y(2)-y(1))/((y(73)+y(74))^2)-p.mu*y(1); %
% Bin 2 
dydt(2) = p.syn_podJ*m_podj+p.syn_podJ2*p.Ji_PodJCtrA^p.n_PodJCtrA/(p.Ji_PodJCtrA^p.n_PodJCtrA+y(2+32)^p.n_PodJCtrA)- (p.deg_podJ1+p.deg_podJm*y(2+64))*y(2)... % synthesis & degradation*PerP (65-68); y(61)=SpodJ 
    - p.dnv_podJ*y(2)...% denovo polymerization; SpmX (13-16&17-20)
    - p.aut1_podJ1/(1+p.alpha_PodJSpmX*(y(2+12)+y(2+16)))*y(2+4)^p.podj*y(2)...% autocatalytic polymerization at centers
+ p.depol_podJ*y(2+4)...%
    +4*p.D_podJm*(y(1)-y(2))/((y(73)+y(74))^2)+p.D_podJm*(y(3)-y(2))/(y(74)^2)- p.mu*y(2);
% Bin 3
dydt(3) = p.syn_podJ*m_podj+p.syn_podJ2*p.Ji_PodJCtrA^p.n_PodJCtrA/(p.Ji_PodJCtrA^p.n_PodJCtrA+y(3+32)^p.n_PodJCtrA)- (p.deg_podJ1+p.deg_podJm*y(3+64))*y(3)... % synthesis & degradation*PerP (65-68); y(61)=SpodJ 
    - p.dnv_podJ*y(3)...% denovo polymerization; SpmX (13-16&17-20)
    - p.aut1_podJ1/(1+p.alpha_PodJSpmX*(y(3+12)+y(3+16)))*y(3+4)^p.podj*y(3)...% autocatalytic polymerization at poles
+ p.depol_podJ*y(3+4)...%
        +4*p.D_podJm*(y(4)-y(3))/((y(73)+y(74))^2)+p.D_podJm*(y(2)-y(3))/(y(74)^2)- p.mu*y(3);
% Bin 4
dydt(4) = 0*p.syn_podJ*m_podj+0*p.syn_podJ2*p.Ji_PodJCtrA^p.n_PodJCtrA/(p.Ji_PodJCtrA^p.n_PodJCtrA+y(4+32)^p.n_PodJCtrA)- (p.deg_podJ1+p.deg_podJm*y(4+64))*y(4)... % synthesis & degradation*PerP (65-68); y(61)=SpodJ 
    - p.dnv_podJ*y(4)...% denovo polymerization; SpmX (13-16&17-20)
    - p.aut1_podJ/(1+p.alpha_PodJSpmX*(y(4+12)+y(4+16)))*y(4+4)^p.podj*y(4)...% autocatalytic polymerization at poles
+ p.depol_podJ*y(4+4)...%
    + 4*p.D_podJm*(y(3)-y(4))/((y(73)+y(74))^2)- p.mu*y(4); 

%% PodJL_p
% Bin 5
dydt(5) =-(p.deg_podJ1+p.deg_podJp*y(5+60))*y(5)...%PerP:65-68
    + p.dnv_podJ*y(5-4)...%SpmX13-16&17-20
    + p.aut1_podJ/(1+p.alpha_PodJSpmX*(y(5+8)+y(5+12)))*y(5)^p.podj*y(5-4)...% autocatalytic polymerization at poles
    - p.depol_podJ*y(5)...%
    + 4*p.D_podJL*(y(6)-y(5))/((y(73)+y(74))^2)- p.mu*y(5); % polymer diffusion % dilution
% Bin 6
dydt(6) =-(p.deg_podJ1+p.deg_podJp*y(6+60))*y(6)...%PerP:65-68
    + p.dnv_podJ*y(6-4)...%SpmX13-16&17-20
    + p.aut1_podJ1/(1+p.alpha_PodJSpmX*(y(6+8)+y(6+12)))*y(6)^p.podj*y(6-4)...% autocatalytic polymerization at centers
    - p.depol_podJ*y(6)...%
    + 4*p.D_podJL*(y(5)-y(6))/((y(73)+y(74))^2)+p.D_podJL*(y(7)-y(6))/(y(74)^2) - p.mu*y(6); 
%Bin 7
dydt(7) =-(p.deg_podJ1+p.deg_podJp*y(7+60))*y(7)...%PerP:65-68
    + p.dnv_podJ*y(7-4)...%SpmX13-16&17-20
    + p.aut1_podJ1/(1+p.alpha_PodJSpmX*(y(7+8)+y(7+12)))*y(7)^p.podj*y(7-4)...% autocatalytic polymerization at centers
    - p.depol_podJ*y(7)...%
        + 4*p.D_podJL*(y(8)-y(7))/((y(73)+y(74))^2)+p.D_podJL*(y(6)-y(7))/(y(74)^2) - p.mu*y(7); 
% Bin 8
dydt(8) =-(p.deg_podJ1+p.deg_podJp*y(8+60))*y(8)...%PerP:65-68
    + p.dnv_podJ*y(8-4)...%SpmX13-16&17-20
    + p.aut1_podJ/(1+p.alpha_PodJSpmX*(y(8+8)+y(8+12)))*y(8)^p.podj*y(8-4)...% autocatalytic polymerization at poles
    - p.depol_podJ*y(8)...%
    + 4*p.D_podJL*(y(7)-y(8))/((y(73)+y(74))^2)-p.mu*y(8); % polymer diffusion % dilution
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PodJS
% no diffussion assumed
for i=9:12
dydt(i)=(p.deg_podJ1+p.deg_podJm*y(i+56))*(y(i-8)+y(i-4))...%
    -(p.mu+p.deg_s)*y(i);
end

%% SpmXm
%Bin 13
q=1;
dydt(13) = 0*p.syn_spmx*y(13+20)^q/(y(13+20)^q+p.Ja_SpmXCtrA^q)-(p.deg_spmx+p.mu)*y(13)...%CtrAP 33-36
    -p.dnv_spmx*y(13)+p.depol_spmx*y(13+4)...
    -p.aut_spmx*(1+p.alpha_SpmXPopZ*(y(13+8)+y(13+12)))*y(13)*y(13+4)^2 ...%PopZp 25-28
+4*p.D_spmx*(y(14)-y(13))/((y(73)+y(74))^2);
%Bin 14
dydt(14) = p.syn_spmx*y(14+20)^q/(y(14+20)^q+p.Ja_SpmXCtrA^q)-(p.deg_spmx+p.mu)*y(14)...%CtrAP 33-36 
    -p.dnv_spmx*y(14)+p.depol_spmx*y(14+4)...
  -p.aut_spmx*(1+p.alpha_SpmXPopZ*(y(14+8)+y(14+12)))*y(14)*y(14+4)^2 ...%PopZp 25-28
    +4*p.D_spmx*(y(13)-y(14))/((y(73)+y(74))^2)+p.D_spmx*(y(15)-y(14))/(y(74)^2);
%Bin 15
dydt(15) = p.syn_spmx*y(15+20)^q/(y(15+20)^q+p.Ja_SpmXCtrA^q)-(p.deg_spmx+p.mu)*y(15)...%CtrAP 33-36 
    -p.dnv_spmx*y(15)+p.depol_spmx*y(15+4)...
  -p.aut_spmx*(1+p.alpha_SpmXPopZ*(y(15+8)+y(15+12)))*y(15)*y(15+4)^2 ...%PopZp 25-28
    +4*p.D_spmx*(y(16)-y(15))/((y(73)+y(74))^2)+p.D_spmx*(y(14)-y(15))/(y(74)^2);
%Bin 16
dydt(16) = 0*p.syn_spmx*y(16+20)^q/(y(16+20)^q+p.Ja_SpmXCtrA^q)-(p.deg_spmx+p.mu)*y(16)...%CtrAP 33-36 
    -p.dnv_spmx*y(16)+p.depol_spmx*y(16+4)...
  -p.aut_spmx*(1+p.alpha_SpmXPopZ*(y(16+8)+y(16+12)))*y(16)*y(16+4)^2 ...%PopZp 25-28
    +4*p.D_spmx*(y(15)-y(16))/((y(73)+y(74))^2);

%% SpmXp
%Bin 17
dydt(17)=-(p.deg_spmx+p.mu)*y(17)...
    +p.dnv_spmx*y(17-4)-p.depol_spmx*y(17)...
      +p.aut_spmx*(1+p.alpha_SpmXPopZ*(y(17+4)+y(17+8)))*y(17-4)*y(17)^2 ...%PopZp 25-28
     +4*p.D_spmxp*(y(18)-y(17))/((y(73)+y(74))^2);
 %Bin 18
dydt(18)=-(p.deg_spmx+p.mu)*y(18)...
    +p.dnv_spmx*y(18-4)-p.depol_spmx*y(18)...
      +p.aut_spmx*(1+p.alpha_SpmXPopZ*(y(18+4)+y(18+8)))*y(18-4)*y(18)^2 ...%PopZp 25-28
      +4*p.D_spmxp*(y(17)-y(18))/((y(73)+y(74))^2)+p.D_spmxp*(y(19)-y(18))/(y(74)^2);
%Bin 19
dydt(19)=-(p.deg_spmx+p.mu)*y(19)...
    +p.dnv_spmx*y(19-4)-p.depol_spmx*y(19)...
      +p.aut_spmx*(1+p.alpha_SpmXPopZ*(y(19+4)+y(19+8)))*y(19-4)*y(19)^2 ...%PopZp 25-28
        +4*p.D_spmxp*(y(20)-y(19))/((y(73)+y(74))^2)+p.D_spmxp*(y(18)-y(19))/(y(74)^2);
%Bin 20
dydt(20)=-(p.deg_spmx+p.mu)*y(20)...
    +p.dnv_spmx*y(20-4)-p.depol_spmx*y(20)...
      +p.aut_spmx*(1+p.alpha_SpmXPopZ*(y(20+4)+y(20+8)))*y(20-4)*y(20)^2 ...%PopZp 25-28
      +4*p.D_spmxp*(y(19)-y(20))/((y(73)+y(74))^2);

%% PopZm
%Bin 21
dydt(21) =0*p.syn_popz  - (p.mu+p.deg_popzm)*y(21)...
    -p.dnv_popz*y(21)-p.aut_popz*(1+p.alpha_PopZPodJ*(y(21-20)+y(21-16)))*y(21+4)^2*y(21)...%PodJp 5-8
    +p.depol_popz*y(21+4)...
    +4*p.D_popzm*(y(22)-y(21))/((y(73)+y(74))^2);
%Bin 22
dydt(22) =p.syn_popz  - (p.mu+p.deg_popzm)*y(22)...
    -p.dnv_popz*y(22)-p.aut_popz1*(1+p.alpha_PopZPodJ*(y(22-20)+y(22-16)))*y(22+4)^2*y(22)...%PodJp 5-8
    +p.depol_popz*y(22+4)...
+4*p.D_popzm*(y(21)-y(22))/((y(73)+y(74))^2)+p.D_popzm*(y(23)-y(22))/(y(74)^2);
%Bin 23
dydt(23) = p.syn_popz -( p.mu+p.deg_popzm)*y(23)...
  -p.dnv_popz*y(23)-p.aut_popz1*(1+p.alpha_PopZPodJ*(y(23-20)+y(23-16)))*y(23+4)^2*y(23)...%PodJp 5-8
    +p.depol_popz*y(23+4)...
+4*p.D_popzm*(y(24)-y(23))/((y(73)+y(74))^2)+p.D_popzm*(y(22)-y(23))/(y(74)^2);
% Bin 24
dydt(24) = 0*p.syn_popz -( p.mu+p.deg_popzm)*y(24)...
  -p.dnv_popz*y(24)-p.aut_popz*(1+p.alpha_PopZPodJ*(y(24-20)+y(24-16)))*y(24+4)^2*y(24)...%PodJp 5-8
    +p.depol_popz*y(24+4)...
+4*p.D_popzm*(y(23)-y(24))/((y(73)+y(74))^2);

%% PopZp
%Bin 25
dydt(25) = -(p.mu+p.deg_popzp)*y(25)...
    +p.dnv_popz*y(25-4)+p.aut_popz*(1+p.alpha_PopZPodJ*(y(25-24)+y(25-20)))*y(25)^2*y(25-4)...%PodJp 5-8
    -p.depol_popz*y(25)...
 +4*p.D_popzp*(y(26)-y(25))/((y(73)+y(74))^2);
%Bin 26
dydt(26) = -(p.mu+p.deg_popzp)*y(26)...
    +p.dnv_popz*y(26-4)+p.aut_popz1*(1+p.alpha_PopZPodJ*(y(26-24)+y(26-20)))*y(26)^2*y(26-4)...%PodJp 5-8
    -p.depol_popz*y(26)...
   +4*p.D_popzp*(y(25)-y(26))/((y(73)+y(74))^2)+p.D_popzp*(y(27)-y(26))/(y(74)^2);
%Bin 27
dydt(27) = -(p.mu+p.deg_popzp)*y(27)...
    +p.dnv_popz*y(27-4)+p.aut_popz1*(1+p.alpha_PopZPodJ*(y(27-24)+y(27-20)))*y(27)^2*y(27-4)...%PodJp 5-8
    -p.depol_popz*y(27)...
  +4*p.D_popzp*(y(28)-y(27))/((y(73)+y(74))^2)+p.D_popzp*(y(26)-y(27))/(y(74)^2);
%Bin 28
dydt(28) = -(p.mu+p.deg_popzp)*y(28)...
    +p.dnv_popz*y(28-4)+p.aut_popz*(1+p.alpha_PopZPodJ*(y(28-24)+y(28-20)))*y(28)^2*y(28-4)...%PodJp 5-8
    -p.depol_popz*y(28)...
   +4*p.D_popzp*(y(27)-y(28))/((y(73)+y(74))^2);

%% CtrA
m_ctrA = (1-m)*y(70)+m; %S_ctrA=y(70)
%Bin 29
dydt(29)=0*p.syn_ctrA1*(1-y(29+4)/(y(29+4)+y(29)+p.Ji_CtrACtrA))*m_ctrA...
    +0*p.syn_ctrA2*y(29+4)/(y(29+4)+y(29)+p.Ja_CtrACtrA) ...
    -(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(29+28)+y(29+32))^2/(p.Jd_CtrA^2+(y(29+28)+y(29+32))^2))*y(29) ...%DivKPT (57~60)+(61~64)
    -p.phoCtrA*y(29)+p.dephoCtrA*(y(29+28)+y(29+32))*y(29+4) ...
+4*p.D_CtrA*(y(30)-y(29))/((y(73)+y(74))^2);
%Bin 30
dydt(30)=p.syn_ctrA1*(1-y(30+4)/(y(30+4)+y(30)+p.Ji_CtrACtrA))*m_ctrA...
    +p.syn_ctrA2*y(30+4)/(y(30+4)+y(30)+p.Ja_CtrACtrA) ...
    -(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(30+28)+y(30+32))^2/(p.Jd_CtrA^2+(y(30+28)+y(30+32))^2))*y(30) ...%DivKP 57-60
    -p.phoCtrA*y(30)+p.dephoCtrA*(y(30+28)+y(30+32))*y(30+4) ...
+4*p.D_CtrA*(y(29)-y(30))/((y(73)+y(74))^2)+p.D_CtrA*(y(31)-y(30))/(y(74)^2);
%Bin 31
dydt(31)=p.syn_ctrA1*(1-y(31+4)/(y(31+4)+y(31)+p.Ji_CtrACtrA))*m_ctrA...
    +p.syn_ctrA2*y(31+4)/(y(31+4)+y(31)+p.Ja_CtrACtrA) ...
    -(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(31+28)+y(31+32))^2/(p.Jd_CtrA^2+(y(31+28)+y(31+32))^2))*y(31) ...%DivKP 57-60
    -p.phoCtrA*y(31)+p.dephoCtrA*(y(31+28)+y(31+32))*y(31+4) ...
+4*p.D_CtrA*(y(32)-y(31))/((y(73)+y(74))^2)+p.D_CtrA*(y(30)-y(31))/(y(74)^2);
%Bin 32
dydt(32)=0*p.syn_ctrA1*(1-y(32+4)/(y(32+4)+y(32)+p.Ji_CtrACtrA))*m_ctrA...
    +0*p.syn_ctrA2*y(32+4)/(y(32+4)+y(32)+p.Ja_CtrACtrA) ...
    -(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(32+28)+y(32+32))^2/(p.Jd_CtrA^2+(y(32+28)+y(32+32))^2))*y(32) ...%DivKP 57-60
    -p.phoCtrA*y(32)+p.dephoCtrA*(y(32+28)+y(32+32))*y(32+4) ...
+4*p.D_CtrA*(y(31)-y(32))/((y(73)+y(74))^2);

%% CtrAP
%Bin 33
 dydt(33)=-(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(33+24)+y(33+28))^2/(p.Jd_CtrA^2+(y(33+24)+y(33+28))^2))*y(33) ...%DivKP 57-60
    +p.phoCtrA*y(33-4)-p.dephoCtrA*(y(33+24)+y(33+28))*y(33) ...
    + 4*p.D_CtrAP*(y(34)-y(33))/((y(73)+y(74))^2);
%Bin 34
 dydt(34)=-(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(34+24)+y(34+28))^2/(p.Jd_CtrA^2+(y(34+24)+y(34+28))^2))*y(34) ...%DivKP 57-60
    +p.phoCtrA*y(34-4)-p.dephoCtrA*(y(34+24)+y(34+28))*y(34) ...
    + 4*p.D_CtrAP*(y(33)-y(34))/((y(73)+y(74))^2)+p.D_CtrAP*(y(35)-y(34))/(y(74)^2);
%Bin 35
 dydt(35)=-(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(35+24)+y(35+28))^2/(p.Jd_CtrA^2+(y(35+24)+y(35+28))^2))*y(35) ...%DivKP 57-60
    +p.phoCtrA*y(35-4)-p.dephoCtrA*(y(35+24)+y(35+28))*y(35) ...
    + 4*p.D_CtrAP*(y(36)-y(35))/((y(73)+y(74))^2)+p.D_CtrAP*(y(34)-y(35))/(y(74)^2);
%Bin 36
 dydt(36)=-(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(36+24)+y(36+28))^2/(p.Jd_CtrA^2+(y(36+24)+y(36+28))^2))*y(36) ...%DivKP 57-60
    +p.phoCtrA*y(36-4)-p.dephoCtrA*(y(36+24)+y(36+28))*y(36) ...
    + 4*p.D_CtrAP*(y(35)-y(36))/((y(73)+y(74))^2);
%% PleCf
m_pleC = (1-mplec)*y(71)+mplec;%S_pleC=y(71)
%Bin 37
dydt(37)=0*p.syn_pleC*m_pleC -(p.mu+p.deg_pleC)*y(37) ...
    -p.fb_PleC*(1+p.alpha_PleCPodJ*(y(37-36)+y(37-32)))*y(37)+p.bf_PleC*y(37+4) ...%PodJp 5-8
    + 4*p.D_PleC*(y(38)-y(37))/((y(73)+y(74))^2);
%Bin 38
dydt(38)=p.syn_pleC*m_pleC -(p.mu+p.deg_pleC)*y(38) ...
    -p.fb_PleC*(1+p.alpha_PleCPodJ*(y(38-36)+y(38-32)))*y(38)+p.bf_PleC*y(38+4) ...%PodJp 5-8
    + 4*p.D_PleC*(y(37)-y(38))/((y(73)+y(74))^2)+p.D_PleC*(y(39)-y(38))/(y(74)^2);
%Bin 39
dydt(39)=p.syn_pleC*m_pleC -(p.mu+p.deg_pleC)*y(39) ...
    -p.fb_PleC*(1+p.alpha_PleCPodJ*(y(39-36)+y(39-32)))*y(39)+p.bf_PleC*y(39+4) ...%PodJp 5-8
    + 4*p.D_PleC*(y(40)-y(39))/((y(73)+y(74))^2)+p.D_PleC*(y(38)-y(39))/(y(74)^2);
%Bin 40
dydt(40)=0*p.syn_pleC*m_pleC -(p.mu+p.deg_pleC)*y(40) ...
    -p.fb_PleC*(1+p.alpha_PleCPodJ*(y(40-36)+y(40-32)))*y(40)+p.bf_PleC*y(40+4) ...%PodJp 5-8
    + 4*p.D_PleC*(y(39)-y(40))/((y(73)+y(74))^2);

%% PleCb
for i=41:44
    dydt(i)=-(p.mu+p.deg_pleC)*y(i) ...
        +p.fb_PleC*(1+p.alpha_PleCPodJ*(y(i-40)+y(i-36)))*y(i-4)-p.bf_PleC*y(i) ;%PodJp 5-8
end

%% DivJf
%Bin 45
dydt(45)=0*p.syn_divJ -(p.mu+p.deg_divJ)*y(45) ...
    -p.fb_DivJ*(1+p.alpha_DivJSpmX*(y(45-32)+y(45-28)))*y(45)+p.bf_DivJ*y(45+4) ...%SpmXp 17-20
    + 4*p.D_DivJ*(y(46)-y(45))/((y(73)+y(74))^2);
%Bin 46
dydt(46)=p.syn_divJ -(p.mu+p.deg_divJ)*y(46) ...
    -p.fb_DivJ*(1+p.alpha_DivJSpmX*(y(46-32)+y(46-28)))*y(46)+p.bf_DivJ*y(46+4) ...%SpmXp 17-20
    + 4*p.D_DivJ*(y(45)-y(46))/((y(73)+y(74))^2)+p.D_DivJ*(y(47)-y(46))/(y(74)^2);
%Bin 47
dydt(47)=p.syn_divJ -(p.mu+p.deg_divJ)*y(47) ...
    -p.fb_DivJ*(1+p.alpha_DivJSpmX*(y(47-32)+y(47-28)))*y(47)+p.bf_DivJ*y(47+4) ...%SpmXp 17-20
    + 4*p.D_DivJ*(y(48)-y(47))/((y(73)+y(74))^2)+p.D_DivJ*(y(46)-y(47))/(y(74)^2);
%Bin 48
dydt(48)=0*p.syn_divJ -(p.mu+p.deg_divJ)*y(48) ...
    -p.fb_DivJ*(1+p.alpha_DivJSpmX*(y(48-32)+y(48-28)))*y(48)+p.bf_DivJ*y(48+4) ...%SpmXp 17-20
    + 4*p.D_DivJ*(y(47)-y(48))/((y(73)+y(74))^2);

%% DivJb
for i=49:52
    dydt(i)=-(p.mu+p.deg_divJ)*y(i) ...
        +p.fb_DivJ*(1+p.alpha_DivJSpmX*(y(i-36)+y(i-32)))*y(i-4)-p.bf_DivJ*y(i) ;%SpmXp 17-20
end
%% DivK
%Bin 53
dydt(53)=0*p.syn_divK1+0*p.syn_divK2*y(53-20)^2/(y(53-20)^2+p.Ja_DivKCtrA^2) ...%CtrAP 33-36
    -(p.mu+p.deg_divK)*y(53) ...
    -p.phoDivK*(p.alpha_DivKDivJ*y(53-4)+1)*y(53)+p.dephoDivK*(y(53-16)+y(53-12))*(y(53+4)+y(53+8)) ...%DivJ 45-48&49-52; PleC 37-40&41-44
+4*p.D_DivK*(y(54)-y(53))/((y(73)+y(74))^2);
%Bin 54
dydt(54)=p.syn_divK1+p.syn_divK2*y(54-20)^2/(y(54-20)^2+p.Ja_DivKCtrA^2) ...%CtrAP 33-36
    -(p.mu+p.deg_divK)*y(54) ...
    -p.phoDivK*(p.alpha_DivKDivJ*y(54-4)+1)*y(54)+p.dephoDivK*(y(54-16)+y(54-12))*(y(54+4)+y(54+8)) ...%DivJ 45-48&49-52; PleC 37-40&41-44
+4*p.D_DivK*(y(53)-y(54))/((y(73)+y(74))^2)+p.D_DivK*(y(55)-y(54))/(y(74)^2);
%Bin 55
dydt(55)=p.syn_divK1+p.syn_divK2*y(55-20)^2/(y(55-20)^2+p.Ja_DivKCtrA^2) ...%CtrAP 33-36
    -(p.mu+p.deg_divK)*y(55) ...
    -p.phoDivK*(p.alpha_DivKDivJ*y(55-4)+1)*y(55)+p.dephoDivK*(y(55-16)+y(55-12))*(y(55+4)+y(55+8)) ...%DivJ (45-48&)DivJb49-52; PleC 37-40&41-44
+4*p.D_DivK*(y(56)-y(55))/((y(73)+y(74))^2)+p.D_DivK*(y(54)-y(55))/(y(74)^2);
%Bin 56
dydt(56)=0*p.syn_divK1+0*p.syn_divK2*y(56-20)^2/(y(56-20)^2+p.Ja_DivKCtrA^2) ...%CtrAP 33-36
    -(p.mu+p.deg_divK)*y(56) ...
    -p.phoDivK*(p.alpha_DivKDivJ*y(56-4)+1)*y(56)+p.dephoDivK*(y(56-16)+y(56-12))*(y(56+4)+y(56+8)) ...%DivJ 45-48&49-52; PleC 37-40&41-44
+4*p.D_DivK*(y(55)-y(56))/((y(73)+y(74))^2);
%% DivKPf
%Bin 57
dydt(57)= -(p.mu+p.deg_divK)*y(57) ...
    +p.phoDivK*(p.alpha_DivKDivJ*y(57-8)+1)*y(57-4)-p.dephoDivK*(y(57-20)+y(57-16))*y(57) ...%DivJ 45-48&49-52; PleC 37-40&41-44
-p.fb_DivKP*y(57)+p.bf_DivKP*y(57+4)+4*p.D_DivKP*(y(58)-y(57))/((y(73)+y(74))^2);
%Bin 58
dydt(58)=(p.mu+p.deg_divK)*y(58) ...
    +p.phoDivK*(p.alpha_DivKDivJ*y(58-8)+1)*y(58-4)-p.dephoDivK*(y(58-20)+y(58-16))*y(58) ...%DivJ 45-48&49-52; PleC 37-40&41-44
-p.fb_DivKP*y(58)+p.bf_DivKP*y(58+4)+4*p.D_DivKP*(y(57)-y(58))/((y(73)+y(74))^2)+p.D_DivKP*(y(59)-y(58))/(y(74)^2);
%Bin 59
dydt(59)=(p.mu+p.deg_divK)*y(59) ...
    +p.phoDivK*(p.alpha_DivKDivJ*y(59-8)+1)*y(59-4)-p.dephoDivK*(y(59-20)+y(59-16))*y(59) ...%DivJ 45-48&49-52; PleC 37-40&41-44
-p.fb_DivKP*y(59)+p.bf_DivKP*y(59+4)+4*p.D_DivKP*(y(60)-y(59))/((y(73)+y(74))^2)+p.D_DivKP*(y(58)-y(59))/(y(74)^2);
%Bin 60
dydt(60)= -(p.mu+p.deg_divK)*y(60) ...
    +p.phoDivK*(p.alpha_DivKDivJ*y(60-8)+1)*y(60-4)-p.dephoDivK*(y(60-20)+y(60-16))*y(60) ...%DivJ 45-48&49-52; PleC 37-40&41-44
-p.fb_DivKP*y(60)+p.bf_DivKP*y(60+4)+4*p.D_DivKP*(y(59)-y(60))/((y(73)+y(74))^2);
%% DivKPb
for i=61:64
% dydt(i)=0;
dydt(i)=-(p.mu+p.deg_divK)*y(i)+p.fb_DivKP*y(i-4)-p.bf_DivKP*y(i);
end%S - PodJ, CtrA, PleC, PerP
%% PerP
m_perP=((1-m)*y(72)+m);
%Bin 65
dydt(65)=0*p.syn_perP*m_perP*y(65-32)^2/(y(65-32)^2+p.Ja_PerPCtrA^2) ...%CtrAP 33-36
    -(p.mu+p.deg_perP)*y(65) ...
    +4*p.D_PerP*(y(66)-y(65))/((y(73)+y(74))^2);
%Bin 66
dydt(66)=p.syn_perP*m_perP*y(66-32)^2/(y(66-32)^2+p.Ja_PerPCtrA^2) ...%CtrAP 33-36
    -(p.mu+p.deg_perP)*y(66) ...
    +4*p.D_PerP*(y(65)-y(66))/((y(73)+y(74))^2)+p.D_PerP*(y(67)-y(66))/(y(74)^2);
%Bin 67
dydt(67)=p.syn_perP*m_perP*y(67-32)^2/(y(67-32)^2+p.Ja_PerPCtrA^2) ...%CtrAP 33-36
    -(p.mu+p.deg_perP)*y(67) ...
    +4*p.D_PerP*(y(68)-y(67))/((y(73)+y(74))^2)+p.D_PerP*(y(66)-y(67))/(y(74)^2);
%Bin 68
dydt(68)=0*p.syn_perP*m_perP*y(68-32)^2/(y(68-32)^2+p.Ja_PerPCtrA^2) ...%CtrAP 33-36
    -(p.mu+p.deg_perP)*y(68) ...
    +4*p.D_PerP*(y(67)-y(68))/((y(73)+y(74))^2);
%% S
for i=69:72
% dydt(i)=0;
dydt(i)=0;
end%S - PodJ, CtrA, PleC, PerP
%% cell growth equation
dydt(73)=p.mu*y(73);
dydt(74) = p.mu*y(74);








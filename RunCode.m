%%runcode
%%1/16/2021

clear all
close all
%% first cycle
load('T_6.mat');%load parameters
%% version of params
ver=6;   

 celltype='SW';
 mutant='WT';%'deltaPodJ';'PodJ+'
 mutant='PodJ+';%'deltaPodJ';
%initial values
y0=zeros(74,1);%SW IC - first cell cycle
    y0(5:8)=10e-6;%PodJp
y0(9:11)=0.001; y0(12)=0.1;%PodJS
y0(25:27)=0.5; y0(28)=2;%25;%PopZp
y0(29:32)=0.2; y0(33:36)=0.5;%CtrA and CtrAP
y0(37:40)=0.05; %PleCf
y0(44)=0.05;%PleCb
y0(53:56)=0.2/2;%DivK
y0(73)=0.02*20; y0(74)=0.02*30;%length of polar and central compartment
if strcmp(mutant,'deltaPodJ')
y0(1:12)=0;
end
if strcmp(celltype,'ST')
     output=main_SW(T,y0,34,ver,mutant);
     yout=output.yout; yout=yout';
     y0=yout(:,end);
end

 CycleNum=5;%# of simulated cycles
 for i=1:CycleNum
   TITLE= [num2str(i) 'cellcyle'];
[Y, time, y0_,TE,IE]=main1(T,y0,celltype,ver,mutant);%simulation
graphcellcycle(Y,time,celltype,mutant,TITLE,1,0,1)%total graph
y0=IniValue(Y,celltype);%update y0 of next cycle
 end

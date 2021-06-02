function y0=IniValue(yout,celltype)
Y0=yout(:,end);
L=Y0(73)+Y0(74);%length of half-cell at the end time point
a=Y0(73)-0.2*L; b=Y0(74)-0.5*L;
if strcmp(celltype,'SW')
for i=1:17
y0_(i*4-4+1)=Y0(i*4-4+2);
y0_(i*4-4+2)=Y0(i*4-4+2);
y0_(i*4-4+4)=Y0(i*4-4+1);
y0_(i*4-4+3)=(Y0(i*4-4+1)*a+Y0(i*4-4+2)*b)/(0.3*L);
end
y0_(73)=0.02*20; y0_(74)=0.02*30;
elseif strcmp(celltype,'ST')
    for i=1:17
    y0_(i*4-4+1)=Y0(i*4-4+3);
y0_(i*4-4+2)=Y0(i*4-4+3);
y0_(i*4-4+4)=Y0(i*4-4+4);
y0_(i*4-4+3)=(Y0(i*4-4+4)*a+Y0(i*4-4+3)*b)/(0.3*L);
    end
y0_(73)=0.024*20; y0_(74)=0.024*30;
end
y0=y0_';
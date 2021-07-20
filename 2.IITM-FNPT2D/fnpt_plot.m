clc
%clear all
close all

%data=load('Output PHIT.dat');

nx=101;
ny=21;
nt=1800;
dt=0.01;

cx=zeros(nx,ny);
cy=zeros(nx,ny);
p=zeros(nx,ny);

l=0;
for i=0:nt
    tval=i*dt;
   
    for j=1:nx
        for k=1:ny
            l=l+1;
            cx(j,k)=data(l,2);
            cy(j,k)=data(l,3);
            p(j,k)=data(l,4);                        
        end
    end
    if(mod(i,10)==0)        
        clf;
        contourf(cx,cy,p)            
        axis([24.5 29.5 -0.1 0.1]);
        title(strcat('Time=',num2str(tval)));               
        pause(0.05)
    end
end
% Simple visualization of deformation gradient on a rectangle
% simple shear
% pure shear
% ref: https://www.continuummechanics.org/deformationgradient.html

clc
close all
clear all
npts = 11;
X = linspace(0,1,npts);
Y = linspace(0,1,npts);
[Xv,Yv] = meshgrid(X,Y);
figure(1)
plot(Xv,Yv,'.',MarkerSize=20,color='k');hold on;

%Def_grad_simple = [0,1;-1,0];
Def_grad_simple = [1,0;0.5,1]
Def_grad_pure = [1,.5;.5,1]
%Def_grad_pure = [0,.5;.5,0]

phi = 60;
Def_grad_rotate = [cosd(phi),-sind(phi);sind(phi),cosd(phi)]
%xy = Def_grad_simple * [X,Y;X,Y]
%x = Def_grad_simple*[Xv,Yv;Xv,Yv]
x = Def_grad_simple(1,1)*Xv + Def_grad_simple(1,2)*Yv;
y = Def_grad_simple(2,1)*Xv + Def_grad_simple(2,2)*Yv;

xp = Def_grad_pure(1,1)*Xv + Def_grad_pure(1,2)*Yv;
yp = Def_grad_pure(2,1)*Xv + Def_grad_pure(2,2)*Yv;

xr = Def_grad_rotate(1,1)*Xv + Def_grad_rotate(1,2)*Yv;
yr = Def_grad_rotate(2,1)*Xv + Def_grad_rotate(2,2)*Yv;

%figure(2)
% x=xy(1,:)
% y=xy(2,:)
% [xv,yv] = meshgrid(x,y)
h0=plot(Xv(1),Yv(1),'.',MarkerSize=20,color='k',DisplayName='undeformed');hold on;
h1=plot(x(1),y(1),'s',MarkerSize=20,color='r',DisplayName='simple'); hold on;
h2=plot(xp(1),yp(1),'h',MarkerSize=20,color='b',DisplayName='pure'); hold on;
h3=plot(xr(1),yr(1),'s',MarkerSize=20,color='k',DisplayName='rigid rotate'); hold on;

h4=plot(x,y,'s',MarkerSize=20,color='r'); hold on;
%legend('simple')
h5=plot(xp,yp,'h',MarkerSize=20,color='b'); hold on;
%legend('pure')

h6=plot(xr,yr,'s',MarkerSize=20,color='k'); hold on;

legend([h0,h1,h2,h3],FontSize=12);



grid on;


x0=100;
y0=100;
width=1000;
height=1000;
set(gcf,'position',[x0,y0,width,height]);
axis square ;
axis([-1 1.5 0 2.5]);

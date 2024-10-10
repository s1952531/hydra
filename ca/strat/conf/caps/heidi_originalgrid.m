%function that creates linspaces over the rectangle domain and 
% generates a grid in the original polygon. 
% date: 7/01/2019

function [X,Xxt,Yxt,xt,yt,vert]=heidi_originalgrid(nx,ny,prevert,ymin,ntheta,nbot)
%nx: the number of x grid lines
%ny: the number of y grid lines
%vert: vector of vertices in the deformed domain
%corn: corners to be chosen when constructing the rectangle in conformal
%map

% Process Inputs and Initialize Defaults:
nargs = 6;
for k = nargin:nargs-1
    switch k
        case 0
            nx = 400;
        case 1
            ny = 200;
        case 2
            prevert = [0.95+i i 0 2-0.5i 2+1.5i 1+1.5i];
        case 3
            ymin=0.2;
        case 4 
            ntheta=36;
        case 5
            nbot=10;
        otherwise
    end
end

%calculate lower boundary:
x3=real(prevert(3));
x4=real(prevert(4));
y3=imag(prevert(3));
y4=imag(prevert(4));
xbot=linspace(x3, x4,nbot+1);
hbot=y4-y3;
wbot=x4-x3;
ybot=y3+hbot*((xbot-x3)/wbot).^2;
zbot=xbot(2:nbot)+i*ybot(2:nbot);

newvert=[prevert(1:3) zbot prevert(4:6)];

%generate polygon vertices:
%calculate radius:
d=real(prevert(6))-real(prevert(1));
r=d/2;

%create an array of theta values:
theta=linspace(0,-pi,ntheta+1);
ctheta=cos(theta); stheta=sin(theta);

%define centre of semi-circle:
x_0=(real(prevert(1))+real(prevert(end)))/2;
y_0=ymin+r;

%define arrays for real and imaginary parts of semi-circle domain:
xh=x_0+r*ctheta;
yh=y_0+r*stheta;
zh=xh+i*yh;

%define vertices:
vert=[newvert zh];

%set corners:
corn=[2 nbot+4 nbot+5 1];

%define a polygon and generate a conformal map of the defomed domain:
poly=polygon(vert);
fun=rectmap(poly,corn);

kx=nx+1;
ky=ny+1;

%Lxt=10.49147;
%Lyt=pi/2.0;
small=1.0e-8;
Lxt=max(imag(fun.prevertex));
Lyt=max(real(fun.prevertex));
xt=linspace(small,Lxt-small,kx);
yt=linspace(Lyt-small,-Lyt+small,ky);

X=zeros(kx,ky);
Xxt=zeros(kx,ky);
Yxt=zeros(kx,ky);

%now calculate derivates:
for n=1:ky
   zt=yt(n)+1i*xt;
   X(:,n)=eval(fun,zt);
   Xdf=evaldiff(fun,zt);
   Xxt(:,n)=-imag(Xdf);
   Yxt(:,n)=real(Xdf);
   %figure(1)
   %plot(X(:,n))
end

%Calculate the derivatives along the edges:
% dx=xt(2)-xt(1);
% dy=yt(2)-yt(1);
% hdxi=1.0/(2.0*dx);
% hdyi=1.0/(2.0*dy);

% Xxt(1,:)=hdxi*real(4.0*X(2,:)-3.0*X(1,:)-X(3,:));
% Yxt(1,:)=hdxi*imag(4.0*X(2,:)-3.0*X(1,:)-X(3,:));
% Xxt(kx,:)=hdxi*real(3.0*X(kx,:)+X(kx-2,:)-4.0*X(kx-1,:));
% Yxt(kx,:)=hdxi*imag(3.0*X(kx,:)+X(kx-2,:)-4.0*X(kx-1,:));
% 
% Xxt(:,1)=-hdyi*imag(4.0*X(:,2)-3.0*X(:,1)-X(:,3));
% Yxt(:,1)=hdyi*real(4.0*X(:,2)-3.0*X(:,1)-X(:,3));
% Xxt(:,ky)=-hdyi*imag(3.0*X(:,ky)+X(:,ky-2)-4.0*X(:,ky-1));
% Yxt(:,ky)=hdyi*real(3.0*X(:,ky)+X(:,ky-2)-4.0*X(:,ky-1));

%Calculate derivatives at end points:
Xxt(1,1)=2.0*Xxt(2,1)-Xxt(3,1);
Xxt(kx,1)=2.0*Xxt(kx-1,1)-Xxt(kx-2,1);
Xxt(1,ky)=2.0*Xxt(2,ky)-Xxt(3,ky);
Xxt(kx,ky)=2.0*Xxt(kx-1,ky)-Xxt(kx-2,ky);

Yxt(1,1)=2.0*Yxt(2,1)-Yxt(3,1);
Yxt(kx,1)=2.0*Yxt(kx-1,1)-Yxt(kx-2,1);
Yxt(1,ky)=2.0*Yxt(2,ky)-Yxt(3,ky);
Yxt(kx,ky)=2.0*Yxt(kx-1,ky)-Yxt(kx-2,ky);























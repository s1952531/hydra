%function that creates linspaces over the rectangle domain and 
% generates a grid in the original polygon. 
% date: 7/01/2019

function [X,Y,Xxt,Yxt,xt,yt]=heidi_originalgrid(nx,ny,vert,corn)
%nx: the number of x grid lines
%ny: the number of y grid lines
%vert: vector of vertices in the deformed domain
%corn: corners to be chosen when constructing the rectangle in conformal
%map

% Process Inputs and Initialize Defaults:
nargs = 4;
for k = nargin:nargs-1
    switch k
        case 0
            nx = 200;
        case 1
            ny = 100;
        case 2
            vert = [0 1 2-0.5i 2+1.5i 1+1.5i 1+0.2i 0.95+0.2i 0.95+i i];
        case 3
            corn = [1 3 4 9];
        otherwise
    end
end

%define a polygon and generate a conformal map of the defomed domain:
poly=polygon(vert)
fun=rectmap(poly,corn)

kx=nx+1;
ky=ny+1;

%Lxt=10.49147;
%Lyt=pi/2.0;
Lxt=max(imag(fun.prevertex))-1e-8
Lyt=max(real(fun.prevertex))
xt=linspace(0,Lxt,kx);
yt=linspace(Lyt,-Lyt,ky);

X=zeros(kx,ky);
Y=zeros(ky,kx);
Xxt=zeros(kx,ky);
Yxt=zeros(kx,ky);

%now create zpx+1 vectors:
hold on 
for n=1:ky
   zt=yt(n)+1i*xt
   X(:,n)=eval(fun,zt);
   Xdf=evaldiff(fun,zt);
   Xxt(:,n)=-imag(Xdf);
   Yxt(:,n)=real(Xdf);
   figure(1)
   plot(X(:,n))
end

for m=1:kx
    zt=1i*xt(m)+yt
    Y(:,m)=eval(fun,zt);
    figure(1)
    plot(Y(:,m))
end
hold off






















%% Separability

clear; clc; close all;

n = 10;
f1 = @(x) exp(-2*x); f2 = @(y) sin(2*pi*y); f = @(x,y) f1(x).*f2(y);
x = linspace(0,1,n); y = linspace(0,1,n);
[X,Y] = meshgrid(x,y);
A = f(X.',Y.'); a1 = f1(x); a2 = f2(y);
norm(A-a1.'*a2)

%% Case 2

n = 10;
f = @(x,y) sin(pi*(x+y));
x = linspace(0,1,n); y = linspace(0,1,n);
[X,Y] = meshgrid(x,y);
A = f(X.',Y.'); 
sv = svd(A);
figure(1)
semilogy(1:n,sv,'o')

%% Case 3
cheb.xy
ff=@(x,y)exp(-40*(x.^2-x.*y+2*y.^2-1/2).^2);
f=chebfun2(ff);
levels = 0.1:0.1:0.9;
contour(f,levels); 
axis([-1 1 -1 1]);
axis square
colormap("hot")
set(gcf,"Color","none")
export_fig chebfun1.pdf

%and its decomposition
levels = 0.1:0.1:0.9;
for k = 1:9
    figure(k+1)
    contour(chebfun2(ff,k),levels)
    xlim([-1 1]), axis equal
    colormap("hot")
    set(gcf,"Color","none")
    export_fig(sprintf("chebfun%d.pdf",k+1));
end

%% the SVD
n = 30;
x = linspace(0,1,n); y = linspace(0,1,n);
[X,Y] = meshgrid(x,y);
A = ff(X.',Y.'); 
sv = svd(A);
figure()
semilogy(1:n,sv,'o')
hline(eps,'k--')
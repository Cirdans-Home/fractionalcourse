close all;

N = 7000;
%% Brownian motion
x = cumsum(randn(N,1));
y = cumsum(randn(N,1));

%% Levy distribution
pd_levy = makedist('Stable','alpha',1.5,'beta',0,'gam',1,'delta',0);
xl = cumsum(random(pd_levy,N,1));
yl = cumsum(random(pd_levy,N,1));

figure(1)
plot(x(1),y(1),'ro',x,y,'b-')
axis square
set(gcf,'Color','None');
export_fig brown.eps

figure(2)
plot(xl(1),yl(1),'ro',xl,yl,'b-')
axis square
set(gcf,'Color','None');
export_fig levy.eps
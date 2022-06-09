%% ML Phase plot

x = linspace(-5,5,100);
y = linspace(-5,5,100);

[X,Y] = meshgrid(x,y);
C = X+1i*Y;

Z = arrayfun( @(x) ml(x,0.5), C);

figure(1)
contourf(X,Y,angle(Z))
title('arg(E_{0.5}(z))')
axis square


figure(2)
contourf(X,Y,log10(abs(Z)))
title('abs(E_{0.5}(z))')
axis square
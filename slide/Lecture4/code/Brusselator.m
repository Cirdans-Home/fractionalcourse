%% Ordinary Brusselator

a = 1 ; mu = 4 ;
param = [ a , mu ] ;
f_fun = @(t,y,par) [ ...
    par(1) - (par(2)+1)*y(1) + y(1)^2*y(2) ; ...
    par(2)*y(1) - y(1)^2*y(2) ] ;
t0 = 0 ; T = 100 ;
y0 = [ 1 ; 1] ;

[T,Y] = ode45(@(t,y) f_fun(t,y,param),[t0,T],y0);

figure(1)
plot(Y(:,1),Y(:,2),'k-',y(1,1),y(2,1),'ro')


%% Fractional Brusselator

clear; clc; close all;


a = 1 ; mu = 4 ;
param = [ a , mu ] ;
f_fun = @(t,y,par) [ ...
    par(1) - (par(2)+1)*y(1) + y(1)^2*y(2) ; ...
    par(2)*y(1) - y(1)^2*y(2) ] ;
J_fun = @(t,y,par) [ ...
    -(par(2)+1) + 2*y(1)*y(2) , y(1)^2 ; ...
    par(2) - 2*y(1)*y(2) , -y(1)^2 ] ;
t0 = 0 ; T = 100 ;
y0 = [ 1 ; 1] ;

alpha = [0.51,0.46] ;
h = 1e-2;
[t, y] = fde_pi1_ex(alpha,f_fun,t0,T,y0,h,param) ;

%% Comparison Brusselator

f_fun = @(t,y,par) [ ...
    par(1) - (par(2)+1)*y(1) + y(1)^2*y(2) ; ...
    par(2)*y(1) - y(1)^2*y(2) ] ;
t0 = 0 ; T = 100 ;
y0 = [ 1 ; 1] ;

[T,Y] = ode45(@(t,y) f_fun(t,y,param),[t0,T],y0);

figure(1)
plot(y(1,:),y(2,:),'k-',...
    Y(:,1),Y(:,2),'r--',...
    y(1,1),y(2,1),'o')

legend('Fractional Brusselator','Ordinary Brusselator','y_0')
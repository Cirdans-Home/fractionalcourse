%% Zolotarev 4th problem

cd('rktoolbox'); addpath(fullfile(cd)); cd("..");
clear; clc; close all;

b = 3; % E = [-b,-1] and F = [1,b]
k = 8; % Degree of rational approximant to sign.
% Solution to Z's fourth problem:
r = rkfun.gallery('sign', k/2, b);  

% Plot the computed rational function:
x = linspace(-5, 5, 1000);
y1 = linspace(-3, -1, 1000); y2 = linspace(1, 3, 1000);
fill([-b -1 -1 -b -b], 1.5*[-1 -1 1 1 -1], .9*[1 1 1] ), hold on
fill([b 1 1 b b], 1.5*[-1 -1 1 1 -1], .9*[1 1 1] )
[~,l1,l2] = plotyy(x, r(x), [y1 0 y2], [(1-abs(r(y1))) NaN (1-abs(r(y2)))]);
l1.LineWidth = 2; l2.LineWidth = 2;
text(-2.1,-1.4,'E','fontsize',16)
text(2,-1.4,'F','fontsize',16)
title(sprintf('best R_{%d%d} approximant to sign',k,k))
xlabel('x'), hold off

%% Zolotarev fourth problem
% Extrema for [-1,-1/b]\cup [1/b,1]:
K = ellipke(1-1/b^2);
[sn, cn, dn] = ellipj((0:k)*K/k, 1-1/b^2);
extrema = b*dn;   % Transplant to [-b,-1]\cup [1,b]
vals = 1-r(extrema);
c = mean( vals(1:2:end) );
e = eig( [ 2-4/c^2 1 ; 1 0 ] );
Zk = min(abs(e))

figure(2)
x = linspace(1, b, 1000);
plot(x, r(x), 'linewidth', 2), hold on,
plot(extrema, r(extrema), '.r', 'markersize', 30)
title('extrema of sign approximation error')
xlabel('x'), hold off

%% Zolotarev extrema
figure(3)
R = @(x) (1 + (1+Zk)/(1-Zk)*r(x))./(1 - (1+Zk)/(1-Zk)*r(x));
x = linspace(-b, b, 5000);
fill([-b -1 -1 -b -b], 1e4*[-1 -1 1 1 -1], .9*[1 1 1] ), hold on
fill([b 1 1 b b], 1e4*[-1 -1 1 1 -1], .9*[1 1 1] )
plot(x, R(x), 'linewidth', 2), ylim([-1e4,1e4])
xlabel('x')
text(-2.1,-1.4,'E','fontsize',16)
text(2,-1.4,'F','fontsize',16)
title('solution to Zolotarev''s third problem'), hold off
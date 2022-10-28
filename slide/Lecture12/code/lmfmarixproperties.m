%% LMF Matrix Properties

clear; clc; close all;

m = 100;
T = linspace(0,10,m); h = T(2)-T(1);
r = zeros(m-1,1);
c = zeros(m-1,1);
r(1:7) = [147/60,-6,15/2,-20/3,15/4,-6/5,1/6];
c(1) = 147/60;
A = toeplitz(r,c);
% Fix boundary conditions
A(1,1) = 1;
A(2,1) = -2; A(2,2) = 3/2;
A(3,1) = 3/2; A(3,2) = -3; A(3,3) = 11/6;
A(4,1) = -4/3; A(4,2) = 3; A(4,3) = -4; A(4,4) = 25/12;
A(5,1) = 5/4; A(5,2) = -10/3; A(5,3) = 5; A(5,4) = -5; A(5,5) = 137/60;

f = fov(A);
[V,ev] = eig(full(A),'vector');
fprintf('k2(X) = %1.2e\n',cond(V,2));

%%
figure(1), clf
hold on
plot(real(ev),imag(ev),'ro','DisplayName','Eigenvalues')
plot(f,'k-','LineWidth',2,'DisplayName','FoV')
hold off
axis tight square
legend('Location','northoutside')
set(gcf,'Color','white')
export_fig(sprintf("fov_m%d.pdf",m))
close all
clear
clc

[ro1,si1] = rosi_bdf(2,fix(2/2)+1);
[ro2,si2] = rosi_bdf(3,fix(3/2)+1);
[ro3,si3] = rosi_bdf(4,fix(4/2)+1);
[ro4,si4] = rosi_bdf(5,fix(5/2)+1);

g1 = @(theta) polyval(ro1,exp(1i*theta))./polyval(si1,exp(1i*theta));
g2 = @(theta) polyval(ro2,exp(1i*theta))./polyval(si2,exp(1i*theta));
g3 = @(theta) polyval(ro3,exp(1i*theta))./polyval(si3,exp(1i*theta));
g4 = @(theta) polyval(ro4,exp(1i*theta))./polyval(si4,exp(1i*theta));

t = 0:1e-2:2*pi;
figure(1)
subplot(1,2,1);
plot(real(g1(t)),imag(g1(t)),'b-',real(g2(t)),imag(g2(t)),'r--',...
    real(g3(t)),imag(g3(t)),'k.-',real(g4(t)),imag(g4(t)),'c:','LineWidth',2)
title('Boundary loci for GBDF formulas','Interpreter','LateX','FontSize',14);
legend('k=2','k=3','k=4','k=5');
axis tight
axis square

[ro1,si1] = rosi_gam(2,1);
[ro2,si2] = rosi_gam(4,2);
[ro3,si3] = rosi_gam(6,3);
[ro4,si4] = rosi_gam(8,4);

g1 = @(theta) polyval(ro1,exp(1i*theta))./polyval(si1,exp(1i*theta));
g2 = @(theta) polyval(ro2,exp(1i*theta))./polyval(si2,exp(1i*theta));
g3 = @(theta) polyval(ro3,exp(1i*theta))./polyval(si3,exp(1i*theta));
g4 = @(theta) polyval(ro4,exp(1i*theta))./polyval(si4,exp(1i*theta));

t = 0:1e-2:2*pi;
figure(1)
subplot(1,2,2);
plot(real(g1(t)),imag(g1(t)),'b-',real(g2(t)),imag(g2(t)),'r--',...
    real(g3(t)),imag(g3(t)),'k.-',real(g4(t)),imag(g4(t)),'c:','LineWidth',2)
title('Boundary loci for GAM formulas','Interpreter','LateX','FontSize',14);
legend('k=2','k=4','k=6','k=8');
axis tight
axis square
%% Tester for fractional integral

clear; clc; close all;

f     = @(t,mu) t.^mu;
Itrue = @(alpha,mu,t) gamma(mu+1)*t.^(alpha+mu)/ gamma(alpha+mu+1);
mu    = 1;
alpha = 1.5;

N = 100;
lambda = 1;
I = constfracint(@(t) f(t,mu),0,1,alpha,N,1);
fprintf('Relative error is: %e\n',abs(I-Itrue(alpha,mu,1))./abs(Itrue(alpha,mu,1)));

err1 = [];
Nval = [10,100,1000,10000,100000];
for N = Nval
   I = constfracint(@(t) f(t,mu),0,1,alpha,N,1);
   err1 =  [err1,abs(I-Itrue(alpha,mu,1))./abs(Itrue(alpha,mu,1))];
end

err0 = [];
Nval = [10,100,1000,10000,100000];
for N = Nval
   I = constfracint(@(t) f(t,mu),0,1,alpha,N,0);
   err0 =  [err0,abs(I-Itrue(alpha,mu,1))./abs(Itrue(alpha,mu,1))];
end

errhalf = [];
Nval = [10,100,1000,10000,100000];
for N = Nval
   I = constfracint(@(t) f(t,mu),0,1,alpha,N,0.5);
   errhalf =  [errhalf,abs(I-Itrue(alpha,mu,1))./abs(Itrue(alpha,mu,1))];
end

figure(1)
loglog(Nval,err1,'o-',Nval,err0,'x--',Nval,errhalf,'s-.','LineWidth',2);
xlabel('N');
ylabel('Relative Error');
legend({'\lambda = 1','\lambda = 0','\lambda = 0.5'},'Location','southwest');


for i=1:length(Nval)
   fprintf("%1.0e & %1.1e & %1.1e & %1.1e \\\\\n",Nval(i),err0(i),err1(i),errhalf(i));
end

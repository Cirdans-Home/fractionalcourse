%% TEST OF GLT PRECONDITIONER

clear; close all; clc;

index = 1;
cval = [];
k = 1;
alpha = 1.9;
for N = 2.^(5:10)
    hN = 1/(N-1);  x = 0:hN:1; dt = hN;
    dplus=@(x)gamma(3-alpha).*x.^alpha;
    dminus=@(x)gamma(3-alpha).*(1-x).^alpha;
    G = glmatrix(N,alpha); % Discretize
    Gr = G; Grt = G.'; I = eye(N,N);
    Dplus = diag(dplus(x));
    Dminus = diag(dminus(x));
    nu = hN^alpha/dt;
    A = nu*I -(Dplus*Gr + Dminus*Grt);
    e = ones(N,1);
    if k == 2
        T = spdiags([-e,2*e,-e],-1:1,N,N);
        P = nu*I + (Dplus*T + Dminus*T.');
    elseif k == 1
        T = spdiags([e,-e],0:1,N,N);
        P = nu*I + (Dplus*T + Dminus*T.');
    end
    evA = svd(A);
    X = P\A;
    evprec = svd(X);
    cval = [cval,cond(X)];
    figure(1)
    hold on
    semilogx(real(evA),index*e,'k+',...
        real(evprec),(index+1)*e,'kx');
    hold off
    index = index+2;
end
yticks(1:2:index)
yticklabels(2.^(5:10))
set(gca,'XScale','log')
ylabel('N')
xlabel('Singular Values')
set(gcf,'Color','none')
export_fig(sprintf("spec_%d_k%d.pdf",10*alpha,k));

figure(2)
N = 2.^(5:10);
fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,0],...
    'Upper',[Inf,Inf],...
    'StartPoint',[1 abs(alpha-k)]);
ft = fittype('a*x^b','options',fo);
[curve,gof] = fit(N.',cval.',ft);
plot(N,N.^abs(alpha-k),'--',N,cval,'o',N,curve(N),':','LineWidth',2)
axis tight
legend({'N^{|\alpha-k|}','K_2(P_{k,N}^{-1}A_N)','Fit'},...
    'Location','best','FontSize',14)
set(gcf,'Color','none')
export_fig(sprintf("cond_%d_k%d.pdf",10*alpha,k));
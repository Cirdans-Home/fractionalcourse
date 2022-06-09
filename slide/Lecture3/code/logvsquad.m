%% N^2 vs N log2(N)^2

k = 2:14;
N = 2.^k;

figure(1)
loglog(N,N.^2,'o-',N,N.*log2(N).^2,'x--','LineWidth',2)
xlabel('N');
grid on
axis square
legend({'N^2','N log_2(N)^2'},'Location','northwest');
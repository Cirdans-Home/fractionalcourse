%% Zeros of the generating function

t = linspace(-pi,pi,1000);
f = @(alpha) -exp(-1i*t).*(1-exp(1i*t)).^alpha;
p = @(alpha) f(alpha) + conj(f(alpha));

figure(1)
plot(t,p(1.2)./max(p(1.2)),t,p(1.5)./max(p(1.5)),...
    t,p(1.8)./max(p(1.8)),t,p(2)./max(p(2)),'LineWidth',2);
legend({'\alpha=1.2','\alpha=1.5','\alpha=1.8','\alpha=2'},...
    'Location','north');
xticks([-pi/10,0,pi/10])
xticklabels({'-\pi/10','0','\pi/10'});
xlim([-pi/10 pi/10])
%axis tight
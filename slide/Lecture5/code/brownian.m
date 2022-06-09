%% The Brownian motion

x = linspace(0,1,20);
dx = x(2)-x(1);

pos(1) = x(10);

for i=1:10
   p = rand(1);
   pos(i+1) = (pos(i)+dx)*(p < 0.5) + (pos(i)-dx)*(p > 0.5); 
end

figure(1)
plot(real(x),imag(x),'ko-','LineWidth',2)
for i=1:length(pos)
    plot(real(x),imag(x),'ko-',real(pos(i)),imag(pos(i)),'rx','LineWidth',2)
    text(real(pos(i)),-0.6,sprintf('W_{%d}',i));
    set(gcf,'Color','None')
    export_fig(sprintf('randomwalk1_%02d.png',i));
end
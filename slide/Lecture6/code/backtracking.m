%% BACKTRACKING

clear; clc; close all;

load('mycielskian4.mat');

A = Problem.A > 0;
D = spdiags(diag(A*A),0,size(A,1),size(A,2));
I = speye(size(A,1),size(A,2));
p2A = A*A - D;
p3A = A*p2A - (D - I)*A;

G = graph(Problem.A,'omitselfloops');

figure(3)
h = plot(G);

fid = fopen('mycielskian4.csv','w+');
fprintf(fid,'id, x, y ,size,color ,opacity,label,IdAsLabel,NoLabel \n');
for i=1:G.numnodes
    fprintf(fid,'%d,%f,%f,.3,orange,1,%d,true,true\n',i,h.XData(i),h.YData(i),i);
end
fclose(fid);

fid = fopen('mycielskian4edge.csv','w+');
fprintf(fid,'u,v,Direct\n');
for i=1:G.numedges
    fprintf(fid,'%d,%d,false\n',G.Edges.EndNodes(i,1),G.Edges.EndNodes(i,2));
end
fclose(fid);

figure(1)
heatmap(A*A);
set(gcf,'Color','None');
export_fig('back2.eps');

figure(2)
heatmap(p2A);
set(gcf,'Color','None');
export_fig('nonback2.eps');

figure(3)
heatmap(A*A*A);
set(gcf,'Color','None');
export_fig('back3.eps');

figure(4)
heatmap(p3A);
set(gcf,'Color','None');
export_fig('nonback3.eps');
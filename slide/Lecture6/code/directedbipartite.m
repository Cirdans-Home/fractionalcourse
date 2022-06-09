%% From Directed to Bipartite

clear; clc; close all;

load('GD01_b.mat');

G = digraph(Problem.A);

figure(1)
h = plot(G)%,'XData',Problem.aux.coord(:,1),'YData',Problem.aux.coord(:,2));

fid = fopen('GD01_b.csv','w+');
fprintf(fid,'id, x, y ,size,color ,opacity,label,IdAsLabel,NoLabel \n');
for i=1:G.numnodes
    fprintf(fid,'%d,%f,%f,.15,orange,1,%d,true,true\n',i,h.XData(i),h.YData(i),i);
end
fclose(fid);

fid = fopen('GD01_bedge.csv','w+');
fprintf(fid,'u,v,Direct\n');
for i=1:G.numedges
    fprintf(fid,'%d,%d,true\n',G.Edges.EndNodes(i,1),G.Edges.EndNodes(i,2));
end
fclose(fid);

O = sparse(size(Problem.A,1),size(Problem.A,2));

Acal = [O,Problem.A;Problem.A',O];
G = graph(Acal);

figure(2)
h = plot(G);

fid = fopen('GD01_b_bip.csv','w+');
fprintf(fid,'id, x, y ,size,color ,opacity,label,IdAsLabel,NoLabel \n');
for i=1:G.numnodes
    fprintf(fid,'%d,%f,%f,.15,orange,1,%d,true,true\n',i,h.XData(i),h.YData(i),i);
end
fclose(fid);

fid = fopen('GD01_bedge_bip.csv','w+');
fprintf(fid,'u,v,Direct\n');
for i=1:G.numedges
    fprintf(fid,'%d,%d,false\n',G.Edges.EndNodes(i,1),G.Edges.EndNodes(i,2));
end
fclose(fid);
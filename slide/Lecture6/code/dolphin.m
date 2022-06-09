clear; clc; close all;

load('dolphins.mat');

G = graph(Problem.A,'omitselfloops');
v = eig(Problem.A);

figure(1)
h = plot(G);

fid = fopen('dolphinssimplevertex.csv','w+');
fprintf(fid,'id, x, y ,size,color ,opacity,label,IdAsLabel,NoLabel \n');
for i=1:G.numnodes
    fprintf(fid,'%d,%f,%f,.3,orange,1,%d,true,true\n',i,h.XData(i),h.YData(i),i);
end
fclose(fid);

fid = fopen('dolphinssimpleedges.csv','w+');
fprintf(fid,'u,v,Direct\n');
for i=1:G.numedges
    fprintf(fid,'%d,%d,false\n',G.Edges.EndNodes(i,1),G.Edges.EndNodes(i,2));
end
fclose(fid);

figure(2)
plot(1:G.numnodes,v,'o')
xlabel('i')
ylabel('\lambda_i')
matlab2tikz('dolphinsvalue.tikz')

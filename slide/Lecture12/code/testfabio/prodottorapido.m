%% TEST PER CALCOLARE PRODOTTO RAPIDO Mv con M somma di kronecker

clear;
clc;

s = 10;
m = 20;

metodo = 0;
passi = 5;
[ALPHA,BETA]=mab(metodo, passi, s);
A1=ALPHA(:,2:s+1); A = [1,zeros(1,s);A1];
B1=BETA(:,2:s+1); B = [1,zeros(1,s);B1];

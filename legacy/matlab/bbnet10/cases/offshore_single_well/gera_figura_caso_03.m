%gera topologia do caso 03
clear all
close all
clc

qgas = 1;
qoleo= 1;
case_03

figure
x=get(gcf,'Position');
set(gcf,'Position',[x(1) x(2) 1.5.*x(3) x(4)])
subplot(1,2,1)
aux=topology(parameters,[1 2 3 4],[ 5 6 7],[1 2 3 4],[5 6]);

subplot(1,2,2)
topology_gas(parameters,[6 8 9 10],[7 8 9],-1500,aux)


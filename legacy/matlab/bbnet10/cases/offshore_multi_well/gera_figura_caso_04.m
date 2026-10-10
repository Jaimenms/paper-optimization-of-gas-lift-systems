%gera topologia do caso 03
clear all
close all
clc

qgas1 = 1;
qgas2 = 1;
qgas3 = 1;
qgas4 = 1;
qoleo1= 1;
qoleo2= 1;
qoleo3= 1;
qoleo4= 1;
case_04

figure
x=get(gcf,'Position');
set(gcf,'Position',[x(1) x(2) 1.5.*x(3) x(4)])
subplot(1,2,1)
aux=topology(parameters,[1 2 3 4],[5 6 7 11 12 13 17 18 19 23 24 25],[1 2 3 4 10 16 22],[5 6 11 12 17 18 23 24]);

subplot(1,2,2)
topology_gas(parameters,[6 8 9 10 12 14 15 16 18 20 21 22 24 26 27 28],[7 8 9 13 14 15 19 20 21 25 26 27],-1500,aux)

